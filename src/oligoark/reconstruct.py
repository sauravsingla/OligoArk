"""Explicit read-similarity graphs and deterministic consensus reconstruction."""

from __future__ import annotations

import math
import time
from collections import Counter
from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable


@dataclass(frozen=True)
class GraphEdge:
    left: int
    right: int
    weight: float


@dataclass(frozen=True)
class SimilarityGraph:
    nodes: tuple[str, ...]
    edges: tuple[GraphEdge, ...]
    candidate_pairs: int = 0

    def connected_components(self) -> tuple[tuple[int, ...], ...]:
        adjacency: list[list[int]] = [[] for _ in self.nodes]
        for edge in self.edges:
            adjacency[edge.left].append(edge.right)
            adjacency[edge.right].append(edge.left)
        visited: set[int] = set()
        components: list[tuple[int, ...]] = []
        for start in range(len(self.nodes)):
            if start in visited:
                continue
            stack = [start]
            component: list[int] = []
            visited.add(start)
            while stack:
                current = stack.pop()
                component.append(current)
                for neighbor in adjacency[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)
            components.append(tuple(sorted(component)))
        return tuple(components)


@dataclass(frozen=True)
class ReconstructionResult:
    consensus_reads: list[str]
    cluster_sizes: list[int]
    edge_count: int = 0
    component_count: int = 0
    node_count: int = 0
    candidate_pairs: int = 0
    consensus_lengths: tuple[int, ...] = ()
    runtime_seconds: float = 0.0


@runtime_checkable
class EdgeScorer(Protocol):
    """Score how likely two reads are to represent the same underlying strand."""

    def score(self, left: str, right: str) -> float:
        """Return a similarity score, conventionally in [0, 1]."""


@runtime_checkable
class ReadReconstructor(Protocol):
    """Convert noisy/duplicated reads into candidate consensus reads."""

    def reconstruct(self, reads: list[str]) -> ReconstructionResult:
        """Return candidate consensus reads and diagnostics."""


@dataclass(frozen=True)
class LevenshteinEdgeScorer:
    """Deterministic normalized-Levenshtein baseline edge scorer."""

    def score(self, left: str, right: str) -> float:
        return normalized_similarity(left, right)


def edit_distance(a: str, b: str) -> int:
    """Compute Levenshtein distance with O(min(len(a), len(b))) memory."""
    if len(a) < len(b):
        a, b = b, a
    previous = list(range(len(b) + 1))
    for i, char_a in enumerate(a, 1):
        current = [i]
        for j, char_b in enumerate(b, 1):
            current.append(
                min(
                    current[-1] + 1,
                    previous[j] + 1,
                    previous[j - 1] + (char_a != char_b),
                )
            )
        previous = current
    return previous[-1]


def normalized_similarity(a: str, b: str) -> float:
    return 1.0 - edit_distance(a, b) / max(1, len(a), len(b))


def global_align(
    reference: str,
    query: str,
    *,
    match_score: int = 2,
    mismatch_score: int = -1,
    gap_score: int = -2,
) -> tuple[str, str]:
    """Needleman-Wunsch global alignment with deterministic tie breaking."""
    rows = len(reference) + 1
    cols = len(query) + 1
    scores = [[0] * cols for _ in range(rows)]
    trace = [[""] * cols for _ in range(rows)]
    for i in range(1, rows):
        scores[i][0] = i * gap_score
        trace[i][0] = "U"
    for j in range(1, cols):
        scores[0][j] = j * gap_score
        trace[0][j] = "L"

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = scores[i - 1][j - 1] + (
                match_score if reference[i - 1] == query[j - 1] else mismatch_score
            )
            up = scores[i - 1][j] + gap_score
            left = scores[i][j - 1] + gap_score
            best = max(diagonal, up, left)
            scores[i][j] = best
            trace[i][j] = "D" if diagonal == best else ("U" if up == best else "L")

    aligned_reference: list[str] = []
    aligned_query: list[str] = []
    i, j = len(reference), len(query)
    while i > 0 or j > 0:
        direction = trace[i][j]
        if i > 0 and j > 0 and direction == "D":
            aligned_reference.append(reference[i - 1])
            aligned_query.append(query[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or direction == "U"):
            aligned_reference.append(reference[i - 1])
            aligned_query.append("-")
            i -= 1
        else:
            aligned_reference.append("-")
            aligned_query.append(query[j - 1])
            j -= 1
    return "".join(reversed(aligned_reference)), "".join(reversed(aligned_query))


def _qgrams(sequence: str, width: int = 5) -> frozenset[str]:
    if width < 1:
        raise ValueError("qgram width must be positive")
    if len(sequence) < width:
        return frozenset({sequence})
    return frozenset(
        sequence[index : index + width]
        for index in range(len(sequence) - width + 1)
    )


def _qgram_counts(sequence: str, width: int) -> Counter[str]:
    if width < 1:
        raise ValueError("qgram width must be positive")
    if len(sequence) < width:
        return Counter({sequence: 1})
    return Counter(
        sequence[index : index + width]
        for index in range(len(sequence) - width + 1)
    )


def _jaccard(left: frozenset[str], right: frozenset[str]) -> float:
    union = left | right
    if not union:
        return 1.0
    return len(left & right) / len(union)


def build_similarity_graph(
    reads: list[str],
    *,
    threshold: float = 0.90,
    scorer: EdgeScorer | None = None,
    qgram_width: int = 5,
    use_qgram_prefilter: bool = True,
) -> SimilarityGraph:
    """Build explicit weighted read graph using all qualifying candidate pairs."""
    if not 0 <= threshold <= 1:
        raise ValueError("threshold must be between 0 and 1")
    if qgram_width < 1:
        raise ValueError("qgram_width must be positive")
    resolved_scorer = scorer or LevenshteinEdgeScorer()
    nodes = tuple(reads)
    edges: list[GraphEdge] = []

    all_pairs = [
        (left, right)
        for left in range(len(nodes))
        for right in range(left + 1, len(nodes))
    ]
    if use_qgram_prefilter:
        signatures = [_qgram_counts(read, qgram_width) for read in nodes]
        postings: dict[str, list[tuple[int, int]]] = {}
        for index, signature in enumerate(signatures):
            for qgram, count in signature.items():
                postings.setdefault(qgram, []).append((index, count))

        shared_counts: dict[tuple[int, int], int] = {}
        for entries in postings.values():
            for offset, (left, left_count) in enumerate(entries):
                for right, right_count in entries[offset + 1 :]:
                    pair = (left, right)
                    shared_counts[pair] = (
                        shared_counts.get(pair, 0) + min(left_count, right_count)
                    )

        candidate_pair_indexes: list[tuple[int, int]] = []
        for left, right in all_pairs:
            left_length = len(nodes[left])
            right_length = len(nodes[right])
            max_length = max(1, left_length, right_length)
            length_similarity = 1.0 - abs(left_length - right_length) / max_length
            if length_similarity < threshold:
                continue
            if min(left_length, right_length) < qgram_width:
                candidate_pair_indexes.append((left, right))
                continue

            max_edit_distance = math.floor(
                (1.0 - threshold) * max_length + 1e-12
            )
            left_qgrams = left_length - qgram_width + 1
            right_qgrams = right_length - qgram_width + 1
            required_shared = math.ceil(
                (
                    left_qgrams
                    + right_qgrams
                    - 2 * qgram_width * max_edit_distance
                )
                / 2
            )
            if required_shared <= 0 or shared_counts.get((left, right), 0) >= required_shared:
                candidate_pair_indexes.append((left, right))
    else:
        candidate_pair_indexes = all_pairs

    for left, right in candidate_pair_indexes:
        max_length = max(1, len(nodes[left]), len(nodes[right]))
        length_similarity = 1.0 - abs(len(nodes[left]) - len(nodes[right])) / max_length
        if length_similarity < threshold:
            continue
        weight = resolved_scorer.score(nodes[left], nodes[right])
        if weight >= threshold:
            edges.append(GraphEdge(left, right, weight))
    return SimilarityGraph(
        nodes=nodes,
        edges=tuple(edges),
        candidate_pairs=len(candidate_pair_indexes),
    )


def _winner(counter: Counter[str], preferred: str = "") -> str:
    if not counter:
        return preferred
    best_count = max(counter.values())
    tied = sorted(key for key, value in counter.items() if value == best_count)
    if preferred in tied:
        return preferred
    return tied[0]


def _medoid(cluster: list[str]) -> str:
    return min(
        cluster,
        key=lambda candidate: (
            sum(edit_distance(candidate, other) for other in cluster),
            len(candidate),
            candidate,
        ),
    )


def medoid_consensus(cluster: list[str]) -> str:
    """Simple non-alignment baseline retained for ablation experiments."""
    if not cluster:
        raise ValueError("cluster must not be empty")
    if len(cluster) == 1:
        return cluster[0]
    medoid = _medoid(cluster)
    same_length = [read for read in cluster if len(read) == len(medoid)]
    if len(same_length) < 2:
        return medoid
    chars: list[str] = []
    for position, preferred in enumerate(medoid):
        vote = Counter(read[position] for read in same_length)
        chars.append(_winner(vote, preferred))
    return "".join(chars)


def _alignment_consensus_with_reference(cluster: list[str], reference: str) -> str:
    """Align every trace to a supplied reference and vote over bases and insertion slots."""
    if not cluster:
        raise ValueError("cluster must not be empty")
    if not reference:
        raise ValueError("reference must not be empty")

    base_votes = [Counter[str]() for _ in reference]
    insertion_votes = [Counter[str]() for _ in range(len(reference) + 1)]
    for read in sorted(cluster):
        aligned_reference, aligned_read = global_align(reference, read)
        per_read_insertions = [""] * (len(reference) + 1)
        reference_position = 0
        for reference_base, read_base in zip(
            aligned_reference,
            aligned_read,
            strict=True,
        ):
            if reference_base == "-":
                if read_base != "-":
                    per_read_insertions[reference_position] += read_base
                continue
            base_votes[reference_position][read_base] += 1
            reference_position += 1
        for slot, insertion in enumerate(per_read_insertions):
            insertion_votes[slot][insertion] += 1

    output: list[str] = []
    for position, preferred in enumerate(reference):
        insertion = _winner(insertion_votes[position], "")
        if insertion:
            output.append(insertion)
        base = _winner(base_votes[position], preferred)
        if base != "-":
            output.append(base)
    tail = _winner(insertion_votes[len(reference)], "")
    if tail:
        output.append(tail)
    return "".join(output)


def alignment_consensus(cluster: list[str]) -> str:
    """Single-pass medoid-anchored alignment consensus retained as a v0.5 baseline."""
    if not cluster:
        raise ValueError("cluster must not be empty")
    if len(cluster) == 1:
        return cluster[0]
    return _alignment_consensus_with_reference(cluster, _medoid(cluster))


def iterative_trace_consensus(cluster: list[str], *, rounds: int = 3) -> str:
    """Iteratively refine a trace consensus using the previous consensus as reference.

    This deterministic baseline is intended for modest insertion/deletion noise. It is
    not a probabilistic trace-reconstruction algorithm and is not claimed state of art.
    """
    if not cluster:
        raise ValueError("cluster must not be empty")
    if rounds < 1:
        raise ValueError("rounds must be positive")
    if len(cluster) == 1:
        return cluster[0]

    reference = _medoid(cluster)
    for _ in range(rounds):
        updated = _alignment_consensus_with_reference(cluster, reference)
        if updated == reference:
            break
        reference = updated
    return reference


@dataclass(frozen=True)
class GraphConsensusReconstructor:
    """Explicit similarity graph + connected components + deterministic consensus."""

    threshold: float = 0.90
    scorer: EdgeScorer = field(default_factory=LevenshteinEdgeScorer)
    qgram_width: int = 5
    use_qgram_prefilter: bool = True
    consensus_mode: str = "alignment"

    def _validate(self) -> None:
        if not 0 <= self.threshold <= 1:
            raise ValueError("threshold must be between 0 and 1")
        if self.qgram_width < 1:
            raise ValueError("qgram_width must be positive")
        if self.consensus_mode not in {"medoid", "alignment"}:
            raise ValueError("consensus_mode must be 'medoid' or 'alignment'")

    def reconstruct(self, reads: list[str]) -> ReconstructionResult:
        self._validate()
        started = time.perf_counter()
        graph = build_similarity_graph(
            reads,
            threshold=self.threshold,
            scorer=self.scorer,
            qgram_width=self.qgram_width,
            use_qgram_prefilter=self.use_qgram_prefilter,
        )
        components = graph.connected_components()
        consensus_function = (
            alignment_consensus if self.consensus_mode == "alignment" else medoid_consensus
        )
        clusters = [[graph.nodes[index] for index in component] for component in components]
        consensus_reads = [consensus_function(cluster) for cluster in clusters]
        return ReconstructionResult(
            consensus_reads=consensus_reads,
            cluster_sizes=[len(cluster) for cluster in clusters],
            edge_count=len(graph.edges),
            component_count=len(components),
            node_count=len(graph.nodes),
            candidate_pairs=graph.candidate_pairs,
            consensus_lengths=tuple(len(read) for read in consensus_reads),
            runtime_seconds=round(time.perf_counter() - started, 8),
        )


@dataclass(frozen=True)
class TraceConsensusReconstructor:
    """Multi-threshold explicit-graph reconstruction with iterative trace consensus.

    Multiple thresholds intentionally generate multiple candidate consensuses. Downstream
    frame CRC, ECC, and archive SHA-256 verification remain the acceptance gate.
    """

    thresholds: tuple[float, ...] = (0.94, 0.90, 0.86)
    scorer: EdgeScorer = field(default_factory=LevenshteinEdgeScorer)
    qgram_width: int = 5
    use_qgram_prefilter: bool = True
    rounds: int = 3
    minimum_component_size: int = 2

    def _validate(self) -> None:
        if not self.thresholds:
            raise ValueError("thresholds must not be empty")
        if any(not 0 <= threshold <= 1 for threshold in self.thresholds):
            raise ValueError("thresholds must be between 0 and 1")
        if self.qgram_width < 1:
            raise ValueError("qgram_width must be positive")
        if self.rounds < 1:
            raise ValueError("rounds must be positive")
        if self.minimum_component_size < 2:
            raise ValueError("minimum_component_size must be at least 2")

    def reconstruct(self, reads: list[str]) -> ReconstructionResult:
        self._validate()
        started = time.perf_counter()
        candidates: list[str] = []
        seen_candidates: set[str] = set()
        processed_components: set[tuple[int, ...]] = set()
        cluster_sizes: list[int] = []
        component_count = 0
        thresholds = sorted(set(self.thresholds), reverse=True)

        base_graph = build_similarity_graph(
            reads,
            threshold=min(thresholds),
            scorer=self.scorer,
            qgram_width=self.qgram_width,
            use_qgram_prefilter=self.use_qgram_prefilter,
        )
        for threshold in thresholds:
            graph = SimilarityGraph(
                nodes=base_graph.nodes,
                edges=tuple(
                    edge for edge in base_graph.edges if edge.weight >= threshold
                ),
                candidate_pairs=base_graph.candidate_pairs,
            )
            components = graph.connected_components()
            component_count += len(components)

            for component in components:
                if (
                    len(component) < self.minimum_component_size
                    or component in processed_components
                ):
                    continue
                processed_components.add(component)
                cluster = [graph.nodes[index] for index in component]
                cluster_sizes.append(len(cluster))
                generated = (
                    iterative_trace_consensus(cluster, rounds=self.rounds),
                    alignment_consensus(cluster),
                    medoid_consensus(cluster),
                )
                for candidate in generated:
                    if candidate not in seen_candidates:
                        seen_candidates.add(candidate)
                        candidates.append(candidate)

        return ReconstructionResult(
            consensus_reads=candidates,
            cluster_sizes=cluster_sizes,
            edge_count=len(base_graph.edges),
            component_count=component_count,
            node_count=len(reads),
            candidate_pairs=base_graph.candidate_pairs,
            consensus_lengths=tuple(len(candidate) for candidate in candidates),
            runtime_seconds=round(time.perf_counter() - started, 8),
        )


def graph_cluster_consensus(
    reads: list[str],
    threshold: float = 0.90,
    *,
    scorer: EdgeScorer | None = None,
    qgram_width: int = 5,
    use_qgram_prefilter: bool = True,
    consensus_mode: str = "alignment",
) -> ReconstructionResult:
    """Compatibility wrapper around GraphConsensusReconstructor."""
    reconstructor = GraphConsensusReconstructor(
        threshold=threshold,
        scorer=scorer or LevenshteinEdgeScorer(),
        qgram_width=qgram_width,
        use_qgram_prefilter=use_qgram_prefilter,
        consensus_mode=consensus_mode,
    )
    return reconstructor.reconstruct(reads)

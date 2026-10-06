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


def _trace_anchor_order(
    cluster: list[str],
    *,
    target_length: int | None,
) -> list[str]:
    """Rank unique observed traces without using the unknown reference sequence."""
    unique = sorted(set(cluster))
    return sorted(
        unique,
        key=lambda candidate: (
            sum(edit_distance(candidate, other) for other in cluster),
            abs(len(candidate) - target_length) if target_length is not None else 0,
            len(candidate),
            candidate,
        ),
    )


def _refine_trace_from_anchor(
    cluster: list[str],
    reference: str,
    *,
    rounds: int,
) -> str:
    current = reference
    for _ in range(rounds):
        updated = _alignment_consensus_with_reference(cluster, current)
        if updated == current:
            break
        current = updated
    return current


def _trace_candidate_score(
    candidate: str,
    cluster: list[str],
    *,
    target_length: int | None,
    length_penalty: float,
) -> tuple[float, int, int, str]:
    """Score a candidate only from observed traces plus an optional known oligo length."""
    length_delta = (
        abs(len(candidate) - target_length) if target_length is not None else 0
    )
    trace_distance = sum(edit_distance(candidate, read) for read in cluster)
    return (
        trace_distance + length_penalty * length_delta * len(cluster),
        length_delta,
        len(candidate),
        candidate,
    )


def multistart_trace_consensus(
    cluster: list[str],
    *,
    target_length: int | None = None,
    anchors: int = 2,
    rounds: int = 2,
    bidirectional: bool = True,
    length_penalty: float = 1.0,
) -> str:
    """Low-compute multi-start trace consensus with optional bidirectional refinement.

    Candidate anchors are ranked only by their agreement with the observed reads. Each
    selected anchor is refined by the existing deterministic alignment consensus. When
    bidirectional mode is enabled, the same procedure is repeated on reversed reads so
    deterministic alignment tie-breaking does not always favor the same end of a strand.
    The final candidate minimizes total edit distance to the observed reads, with an
    optional penalty for deviating from a caller-supplied known oligo length.

    This is a lightweight deterministic research baseline; it is not a copy of BBS and
    does not use the unknown reference sequence during reconstruction.
    """
    if not cluster:
        raise ValueError("cluster must not be empty")
    if anchors < 1:
        raise ValueError("anchors must be positive")
    if rounds < 1:
        raise ValueError("rounds must be positive")
    if target_length is not None and target_length < 1:
        raise ValueError("target_length must be positive")
    if length_penalty < 0:
        raise ValueError("length_penalty must be non-negative")
    if len(cluster) == 1:
        return cluster[0]

    candidates: list[str] = []
    seen: set[str] = set()

    def add_candidate(candidate: str) -> None:
        if candidate not in seen:
            seen.add(candidate)
            candidates.append(candidate)

    forward_anchors = _trace_anchor_order(
        cluster,
        target_length=target_length,
    )[:anchors]
    for anchor in forward_anchors:
        add_candidate(anchor)
        add_candidate(
            _refine_trace_from_anchor(
                cluster,
                anchor,
                rounds=rounds,
            )
        )

    if bidirectional:
        reversed_cluster = [read[::-1] for read in cluster]
        reverse_anchors = _trace_anchor_order(
            reversed_cluster,
            target_length=target_length,
        )[:anchors]
        for anchor in reverse_anchors:
            reversed_candidate = _refine_trace_from_anchor(
                reversed_cluster,
                anchor,
                rounds=rounds,
            )
            add_candidate(reversed_candidate[::-1])

    return min(
        candidates,
        key=lambda candidate: _trace_candidate_score(
            candidate,
            cluster,
            target_length=target_length,
            length_penalty=length_penalty,
        ),
    )


def _alignment_vote_evidence(
    reference: str,
    cluster: list[str],
) -> tuple[list[Counter[str]], list[Counter[str]]]:
    """Collect base and insertion votes after aligning reads to a candidate."""
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
    return base_votes, insertion_votes


def _homopolymer_run_after_insertion(
    candidate: str,
    slot: int,
    base: str,
) -> int:
    left = 0
    index = slot - 1
    while index >= 0 and candidate[index] == base:
        left += 1
        index -= 1
    right = 0
    index = slot
    while index < len(candidate) and candidate[index] == base:
        right += 1
        index += 1
    return left + 1 + right


def _targeted_insertion_operations(
    candidate: str,
    insertion_votes: list[Counter[str]],
    *,
    min_homopolymer_run: int,
) -> set[tuple[int, str]]:
    operations: set[tuple[int, str]] = set()
    dna_bases = set("ACGT")

    for slot, votes in enumerate(insertion_votes):
        nonempty = [sequence for sequence in votes if sequence]
        if not nonempty:
            continue
        for sequence in nonempty:
            operations.update((slot, base) for base in sequence if base in dna_bases)
        if slot > 0:
            operations.add((slot, candidate[slot - 1]))
        if slot < len(candidate):
            operations.add((slot, candidate[slot]))

    for slot in range(len(candidate) + 1):
        if slot < len(candidate):
            base = candidate[slot]
            run = 0
            index = slot
            while index < len(candidate) and candidate[index] == base:
                run += 1
                index += 1
            if run >= min_homopolymer_run:
                operations.add((slot, base))
        if slot > 0:
            base = candidate[slot - 1]
            run = 0
            index = slot - 1
            while index >= 0 and candidate[index] == base:
                run += 1
                index -= 1
            if run >= min_homopolymer_run:
                operations.add((slot, base))

    if not operations:
        for slot in range(len(candidate) + 1):
            for base in "ACGT":
                operations.add((slot, base))
    return operations


def targeted_trace_consensus(
    cluster: list[str],
    *,
    target_length: int,
    anchors: int = 3,
    rounds: int = 1,
    bidirectional: bool = True,
    length_penalty: float = 1.0,
    homopolymer_weight: float = 1.0,
    min_homopolymer_run: int = 2,
    substitution_min_gain: float = 0.0,
    substitution_homopolymer_weight: float = 0.0,
    substitution_homopolymer_min_reads: int = 10,
) -> str:
    """Bounded confidence-guided one-edit repair over multi-start consensus.

    The method starts from the lightweight multi-start consensus. A one-base-short
    candidate receives a bounded insertion repair using observed insertion evidence
    and homopolymer context. A target-length candidate tests substitutions only at
    its single most ambiguous aligned column. Selection uses observed reads and the
    known oligo length; the unknown reference is never used.
    """
    if not cluster:
        raise ValueError("cluster must not be empty")
    if target_length < 1:
        raise ValueError("target_length must be positive")
    if anchors < 1:
        raise ValueError("anchors must be positive")
    if rounds < 1:
        raise ValueError("rounds must be positive")
    if length_penalty < 0:
        raise ValueError("length_penalty must be non-negative")
    if homopolymer_weight < 0:
        raise ValueError("homopolymer_weight must be non-negative")
    if min_homopolymer_run < 1:
        raise ValueError("min_homopolymer_run must be positive")
    if substitution_min_gain < 0:
        raise ValueError("substitution_min_gain must be non-negative")
    if substitution_homopolymer_weight < 0:
        raise ValueError("substitution_homopolymer_weight must be non-negative")
    if substitution_homopolymer_min_reads < 1:
        raise ValueError("substitution_homopolymer_min_reads must be positive")
    if len(cluster) == 1:
        return cluster[0]

    selected = multistart_trace_consensus(
        cluster,
        target_length=target_length,
        anchors=anchors,
        rounds=rounds,
        bidirectional=bidirectional,
        length_penalty=length_penalty,
    )
    distance_cache: dict[str, int] = {}

    def read_distance_sum(candidate: str) -> int:
        cached = distance_cache.get(candidate)
        if cached is not None:
            return cached
        value = sum(edit_distance(candidate, read) for read in cluster)
        distance_cache[candidate] = value
        return value

    base_votes, insertion_votes = _alignment_vote_evidence(selected, cluster)

    if len(selected) == target_length - 1:
        operations = _targeted_insertion_operations(
            selected,
            insertion_votes,
            min_homopolymer_run=min_homopolymer_run,
        )
        scored: list[tuple[float, int, str]] = []
        for slot, base in operations:
            candidate = selected[:slot] + base + selected[slot:]
            read_score = read_distance_sum(candidate)
            insertion_run_length = _homopolymer_run_after_insertion(
                selected,
                slot,
                base,
            )
            scored.append(
                (
                    read_score - homopolymer_weight * insertion_run_length,
                    read_score,
                    candidate,
                )
            )
        return min(scored)[2]

    if len(selected) != target_length or not selected:
        return selected

    uncertainty: list[tuple[int, int, int, int]] = []
    for position, votes in enumerate(base_votes):
        selected_count = votes[selected[position]]
        alternative_count = max(
            (
                count
                for base, count in votes.items()
                if base != selected[position]
            ),
            default=0,
        )
        uncertainty.append(
            (
                selected_count - alternative_count,
                -alternative_count,
                selected_count,
                position,
            )
        )
    position = min(uncertainty)[3]
    baseline_score = read_distance_sum(selected)

    def run_length_at_position(sequence: str, index: int) -> int:
        base = sequence[index]
        left = index
        while left > 0 and sequence[left - 1] == base:
            left -= 1
        right = index
        while right + 1 < len(sequence) and sequence[right + 1] == base:
            right += 1
        return right - left + 1

    active_homopolymer_weight = (
        substitution_homopolymer_weight
        if len(cluster) >= substitution_homopolymer_min_reads
        else 0.0
    )
    baseline_composite = (
        baseline_score
        - active_homopolymer_weight * run_length_at_position(selected, position)
    )
    alternatives: list[tuple[float, int, str]] = []
    for base in "ACGT":
        if base == selected[position]:
            continue
        candidate = selected[:position] + base + selected[position + 1 :]
        candidate_score = read_distance_sum(candidate)
        candidate_composite = (
            candidate_score
            - active_homopolymer_weight
            * run_length_at_position(candidate, position)
        )
        alternatives.append((candidate_composite, candidate_score, candidate))

    best_composite, _, best_candidate = min(alternatives)
    if baseline_composite - best_composite >= substitution_min_gain:
        return best_candidate
    return selected


def _candidate_qgram_similarity(
    candidate: str,
    cluster: list[str],
    *,
    width: int,
) -> float:
    signature = _qgrams(candidate, width)
    if not cluster:
        return 0.0
    return sum(_jaccard(signature, _qgrams(read, width)) for read in cluster) / len(cluster)


def _bounded_local_candidates(
    candidate: str,
    cluster: list[str],
    *,
    target_length: int,
    top_positions: int,
    max_candidates: int,
) -> list[str]:
    """Generate a tiny deterministic edit neighborhood from low-confidence evidence."""
    if top_positions < 1 or max_candidates < 1:
        return []

    base_votes, insertion_votes = _alignment_vote_evidence(candidate, cluster)
    generated: list[str] = []
    seen: set[str] = {candidate}

    def add(value: str) -> None:
        if value not in seen and len(generated) < max_candidates:
            seen.add(value)
            generated.append(value)

    def ranked_insertions(
        sequence: str,
        votes_by_slot: list[Counter[str]],
    ) -> list[tuple[int, int, int, str]]:
        ranked: list[tuple[int, int, int, str]] = []
        operations = _targeted_insertion_operations(
            sequence,
            votes_by_slot,
            min_homopolymer_run=2,
        )
        for slot, base in operations:
            observed_support = max(
                (
                    count
                    for inserted, count in votes_by_slot[slot].items()
                    if inserted and base in inserted
                ),
                default=0,
            )
            run_length = _homopolymer_run_after_insertion(sequence, slot, base)
            ranked.append((-observed_support, -run_length, slot, base))
        return sorted(set(ranked))

    if len(candidate) == target_length - 1:
        for _, _, slot, base in ranked_insertions(candidate, insertion_votes):
            add(candidate[:slot] + base + candidate[slot:])
        return generated

    if len(candidate) == target_length - 2:
        beam_width = max(2, min(4, top_positions + 1))
        first_steps = ranked_insertions(candidate, insertion_votes)[:beam_width]
        for _, _, slot, base in first_steps:
            first = candidate[:slot] + base + candidate[slot:]
            _, second_insertion_votes = _alignment_vote_evidence(first, cluster)
            second_steps = ranked_insertions(first, second_insertion_votes)[:beam_width]
            for _, _, second_slot, second_base in second_steps:
                add(first[:second_slot] + second_base + first[second_slot:])
        return generated

    if len(candidate) == target_length + 1:
        length_deletions: list[tuple[int, int, int]] = []
        for position, votes in enumerate(base_votes):
            selected_count = votes[candidate[position]]
            gap_count = votes["-"]
            length_deletions.append(
                (selected_count - gap_count, selected_count, position)
            )
        for _, _, position in sorted(length_deletions)[:top_positions]:
            add(candidate[:position] + candidate[position + 1 :])
        return generated

    if len(candidate) != target_length:
        return generated

    uncertainty: list[tuple[int, int, int, int]] = []
    for position, votes in enumerate(base_votes):
        selected_count = votes[candidate[position]]
        alternative_count = max(
            (count for base, count in votes.items() if base != candidate[position]),
            default=0,
        )
        uncertainty.append(
            (
                selected_count - alternative_count,
                -alternative_count,
                selected_count,
                position,
            )
        )

    suspicious_positions = [entry[3] for entry in sorted(uncertainty)[:top_positions]]
    for position in suspicious_positions:
        votes = base_votes[position]
        alternatives = sorted(
            (
                (-count, base)
                for base, count in votes.items()
                if base in "ACGT" and base != candidate[position]
            )
        )
        if not alternatives:
            alternatives = [(0, base) for base in "ACGT" if base != candidate[position]]
        for _, base in alternatives[:1]:
            add(candidate[:position] + base + candidate[position + 1 :])

    # Tiny same-length indel-pair beam. High gap support nominates a deletion;
    # insertion evidence/homopolymer context nominates the compensating insertion.
    ranked_deletions: list[tuple[int, int, int, int]] = []
    for position, votes in enumerate(base_votes):
        selected_count = votes[candidate[position]]
        gap_count = votes["-"]
        ranked_deletions.append(
            (
                selected_count - gap_count,
                -gap_count,
                selected_count,
                position,
            )
        )
    deletion_positions = [
        entry[3] for entry in sorted(ranked_deletions)[:top_positions]
    ]
    insertion_operations = ranked_insertions(candidate, insertion_votes)[:top_positions]
    for delete_position in deletion_positions:
        without = candidate[:delete_position] + candidate[delete_position + 1 :]
        for _, _, slot, base in insertion_operations:
            adjusted_slot = slot - 1 if slot > delete_position else slot
            adjusted_slot = max(0, min(len(without), adjusted_slot))
            add(without[:adjusted_slot] + base + without[adjusted_slot:])

    # Very local shift repair remains as a final cheap fallback.
    for position in suspicious_positions:
        without = candidate[:position] + candidate[position + 1 :]
        for slot in (max(0, position - 1), min(len(without), position + 1)):
            suggested: set[str] = {candidate[position]}
            source_slot = min(slot, len(insertion_votes) - 1)
            suggested.update(
                base
                for sequence in insertion_votes[source_slot]
                if sequence
                for base in sequence
                if base in "ACGT"
            )
            if slot > 0:
                suggested.add(without[slot - 1])
            if slot < len(without):
                suggested.add(without[slot])
            for base in sorted(suggested):
                add(without[:slot] + base + without[slot:])

    return generated


def confidence_fusion_trace_consensus(
    cluster: list[str],
    *,
    target_length: int,
    anchors: int = 3,
    rounds: int = 1,
    top_positions: int = 3,
    max_candidates: int = 24,
    trim_farthest: int = 1,
    qgram_width: int = 4,
    qgram_weight: float = 0.5,
    minimum_score_gain: float = 0.0,
) -> str:
    """Fuse low-cost consensuses and apply a bounded confidence-guided edit search.

    The method never uses the unknown reference. It starts from the frozen targeted
    repair, adds only a tiny edit neighborhood around low-confidence positions, and
    selects a replacement only when observed-read evidence improves by a configured margin.
    """
    if not cluster:
        raise ValueError("cluster must not be empty")
    if target_length < 1:
        raise ValueError("target_length must be positive")
    if anchors < 1 or rounds < 1:
        raise ValueError("anchors and rounds must be positive")
    if top_positions < 1 or max_candidates < 1:
        raise ValueError("top_positions and max_candidates must be positive")
    if trim_farthest < 0 or trim_farthest >= len(cluster):
        raise ValueError("trim_farthest must be between 0 and len(cluster)-1")
    if qgram_width < 1:
        raise ValueError("qgram_width must be positive")
    if qgram_weight < 0 or minimum_score_gain < 0:
        raise ValueError("weights and minimum_score_gain must be non-negative")
    if len(cluster) == 1:
        return cluster[0]

    baseline = targeted_trace_consensus(
        cluster,
        target_length=target_length,
        anchors=anchors,
        rounds=rounds,
        bidirectional=True,
        length_penalty=1.0,
        homopolymer_weight=0.75,
        min_homopolymer_run=2,
        substitution_min_gain=0.0,
        substitution_homopolymer_weight=2.0,
        substitution_homopolymer_min_reads=10,
    )

    # The frozen targeted consensus already includes multi-start and bidirectional
    # refinement. Keep this second-stage search deliberately tiny so the external
    # benchmark remains CPU-friendly.
    candidates: list[str] = [baseline]
    seen: set[str] = {baseline}
    for local in _bounded_local_candidates(
        baseline,
        cluster,
        target_length=target_length,
        top_positions=top_positions,
        max_candidates=max_candidates - 1,
    ):
        if local not in seen:
            seen.add(local)
            candidates.append(local)
        if len(candidates) >= max_candidates:
            break

    distance_cache: dict[str, tuple[int, int]] = {}

    def distance_scores(candidate: str) -> tuple[int, int]:
        cached = distance_cache.get(candidate)
        if cached is not None:
            return cached
        distances = sorted(edit_distance(candidate, read) for read in cluster)
        robust = (
            sum(distances[: len(distances) - trim_farthest])
            if trim_farthest
            else sum(distances)
        )
        total = sum(distances)
        distance_cache[candidate] = (robust, total)
        return robust, total

    def score(candidate: str) -> tuple[float, int, int, str]:
        robust, total = distance_scores(candidate)
        length_delta = abs(len(candidate) - target_length)
        qgram_similarity = _candidate_qgram_similarity(
            candidate,
            cluster,
            width=qgram_width,
        )
        composite = (
            robust
            + 0.05 * total
            + 2.0 * length_delta
            - qgram_weight * qgram_similarity
        )
        return (composite, length_delta, total, candidate)

    baseline_score = score(baseline)[0]
    best = min(candidates, key=score)
    best_score = score(best)[0]
    if baseline_score - best_score >= minimum_score_gain:
        return best
    return baseline


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

"""Extensible read-reconstruction interfaces and deterministic graph baselines."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable


@dataclass(frozen=True)
class ReconstructionResult:
    """Consensus reads plus diagnostic cluster sizes."""

    consensus_reads: list[str]
    cluster_sizes: list[int]


@runtime_checkable
class EdgeScorer(Protocol):
    """Score how likely two reads are to represent the same underlying strand."""

    def score(self, left: str, right: str) -> float:
        """Return a similarity score, conventionally in the inclusive range [0, 1]."""


@runtime_checkable
class ReadReconstructor(Protocol):
    """Convert noisy/duplicated reads into candidate consensus reads."""

    def reconstruct(self, reads: list[str]) -> ReconstructionResult:
        """Return candidate consensus reads and reconstruction diagnostics."""


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
    """Return normalized Levenshtein similarity in [0, 1]."""
    return 1.0 - edit_distance(a, b) / max(1, len(a), len(b))


def _qgrams(sequence: str, width: int = 5) -> frozenset[str]:
    if width < 1:
        raise ValueError("qgram width must be positive")
    if len(sequence) < width:
        return frozenset({sequence})
    return frozenset(sequence[index : index + width] for index in range(len(sequence) - width + 1))


def _jaccard(left: frozenset[str], right: frozenset[str]) -> float:
    union = left | right
    if not union:
        return 1.0
    return len(left & right) / len(union)


def _consensus(cluster: list[str]) -> str:
    """Use a medoid plus same-length majority voting as a deterministic baseline."""
    if len(cluster) == 1:
        return cluster[0]
    medoid = min(
        cluster,
        key=lambda candidate: sum(edit_distance(candidate, other) for other in cluster),
    )
    same_len = [read for read in cluster if len(read) == len(medoid)]
    if len(same_len) < 2:
        return medoid
    chars: list[str] = []
    for position in range(len(medoid)):
        vote = Counter(read[position] for read in same_len)
        chars.append(vote.most_common(1)[0][0])
    return "".join(chars)


@dataclass(frozen=True)
class GraphConsensusReconstructor:
    """Implicit similarity-graph clustering with a pluggable edge scorer.

    The default scorer is deterministic normalized Levenshtein similarity. A future GNN can
    implement EdgeScorer and set use_qgram_prefilter=False if the learned scorer should evaluate
    all candidate pairs.
    """

    threshold: float = 0.90
    scorer: EdgeScorer = field(default_factory=LevenshteinEdgeScorer)
    qgram_width: int = 5
    use_qgram_prefilter: bool = True

    def _validate(self) -> None:
        if not 0 <= self.threshold <= 1:
            raise ValueError("threshold must be between 0 and 1")
        if self.qgram_width < 1:
            raise ValueError("qgram_width must be positive")

    def reconstruct(self, reads: list[str]) -> ReconstructionResult:
        self._validate()
        if not reads:
            return ReconstructionResult([], [])

        clusters: list[list[str]] = []
        representatives: list[str] = []
        signatures: list[frozenset[str]] = []
        prefilter_threshold = max(0.10, self.threshold - 0.35)

        for read in sorted(reads, key=lambda value: (len(value), value)):
            read_signature = _qgrams(read, self.qgram_width)
            best_index = -1
            best_score = self.threshold
            for index, representative in enumerate(representatives):
                max_length = max(1, len(read), len(representative))
                length_similarity = 1.0 - abs(len(read) - len(representative)) / max_length
                if length_similarity < self.threshold:
                    continue
                if (
                    self.use_qgram_prefilter
                    and _jaccard(read_signature, signatures[index]) < prefilter_threshold
                ):
                    continue
                score = self.scorer.score(read, representative)
                if score >= best_score:
                    best_index, best_score = index, score
            if best_index < 0:
                clusters.append([read])
                representatives.append(read)
                signatures.append(read_signature)
            else:
                clusters[best_index].append(read)

        return ReconstructionResult(
            consensus_reads=[_consensus(cluster) for cluster in clusters],
            cluster_sizes=[len(cluster) for cluster in clusters],
        )


def graph_cluster_consensus(
    reads: list[str],
    threshold: float = 0.90,
    *,
    scorer: EdgeScorer | None = None,
    qgram_width: int = 5,
    use_qgram_prefilter: bool = True,
) -> ReconstructionResult:
    """Compatibility wrapper around GraphConsensusReconstructor."""
    reconstructor = GraphConsensusReconstructor(
        threshold=threshold,
        scorer=scorer or LevenshteinEdgeScorer(),
        qgram_width=qgram_width,
        use_qgram_prefilter=use_qgram_prefilter,
    )
    return reconstructor.reconstruct(reads)

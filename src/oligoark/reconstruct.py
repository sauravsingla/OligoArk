"""Graph-inspired read clustering and consensus reconstruction baselines."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class ReconstructionResult:
    consensus_reads: list[str]
    cluster_sizes: list[int]


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


def _qgrams(sequence: str, width: int = 5) -> frozenset[str]:
    if len(sequence) < width:
        return frozenset({sequence})
    return frozenset(sequence[index : index + width] for index in range(len(sequence) - width + 1))


def _jaccard(left: frozenset[str], right: frozenset[str]) -> float:
    union = left | right
    if not union:
        return 1.0
    return len(left & right) / len(union)


def _consensus(cluster: list[str]) -> str:
    """Use a medoid plus same-length majority voting as an indel-aware baseline."""
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


def graph_cluster_consensus(reads: list[str], threshold: float = 0.90) -> ReconstructionResult:
    """Cluster an implicit similarity graph and return deterministic consensus reads.

    A q-gram Jaccard prefilter avoids quadratic edit-distance work for unrelated strands while
    preserving full Levenshtein scoring for plausible neighbors. This remains a deterministic
    non-ML baseline; a future GNN can replace edge scoring without changing the archive codec.
    """
    if not 0 <= threshold <= 1:
        raise ValueError("threshold must be between 0 and 1")

    clusters: list[list[str]] = []
    representatives: list[str] = []
    signatures: list[frozenset[str]] = []
    prefilter_threshold = max(0.10, threshold - 0.35)

    for read in sorted(reads, key=lambda value: (len(value), value)):
        read_signature = _qgrams(read)
        best_index = -1
        best_score = threshold
        for index, representative in enumerate(representatives):
            max_length = max(1, len(read), len(representative))
            length_similarity = 1.0 - abs(len(read) - len(representative)) / max_length
            if length_similarity < threshold:
                continue
            if _jaccard(read_signature, signatures[index]) < prefilter_threshold:
                continue
            score = normalized_similarity(read, representative)
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

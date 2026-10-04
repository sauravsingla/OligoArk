"""Graph-inspired read clustering and consensus reconstruction baselines."""

from __future__ import annotations
from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class ReconstructionResult:
    consensus_reads: list[str]
    cluster_sizes: list[int]


def edit_distance(a: str, b: str) -> int:
    if len(a) < len(b):
        a, b = b, a
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        current = [i]
        for j, cb in enumerate(b, 1):
            current.append(min(current[-1] + 1, previous[j] + 1,
                               previous[j - 1] + (ca != cb)))
        previous = current
    return previous[-1]


def normalized_similarity(a: str, b: str) -> float:
    return 1.0 - edit_distance(a, b) / max(1, len(a), len(b))


def _consensus(cluster: list[str]) -> str:
    if len(cluster) == 1:
        return cluster[0]
    medoid = min(cluster, key=lambda c: sum(edit_distance(c, other) for other in cluster))
    same_len = [read for read in cluster if len(read) == len(medoid)]
    if len(same_len) < 2:
        return medoid
    chars: list[str] = []
    for pos in range(len(medoid)):
        chars.append(Counter(read[pos] for read in same_len).most_common(1)[0][0])
    return "".join(chars)


def graph_cluster_consensus(reads: list[str], threshold: float = 0.90) -> ReconstructionResult:
    """Greedy connected-component baseline over an implicit sequence-similarity graph."""
    clusters: list[list[str]] = []
    for read in sorted(reads, key=lambda r: (len(r), r)):
        best_idx, best_score = -1, threshold
        for idx, cluster in enumerate(clusters):
            score = normalized_similarity(read, cluster[0])
            if score >= best_score:
                best_idx, best_score = idx, score
        if best_idx < 0:
            clusters.append([read])
        else:
            clusters[best_idx].append(read)
    return ReconstructionResult([_consensus(c) for c in clusters], [len(c) for c in clusters])

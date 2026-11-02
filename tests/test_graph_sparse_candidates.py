"""Graph q-gram candidate enumeration must preserve exhaustive graph edges."""

import math
import random
from collections import Counter
from itertools import combinations

import pytest

from oligoark.reconstruct import build_similarity_graph


def reference_prefilter_count(reads: list[str], threshold: float, width: int) -> int:
    """Quadratic reference oracle for small, deterministic regression inputs."""
    total = 0
    for left, right in combinations(reads, 2):
        longest = max(1, len(left), len(right))
        if 1 - abs(len(left) - len(right)) / longest < threshold:
            continue
        if min(len(left), len(right)) >= width:
            edits = math.floor((1 - threshold) * longest + 1e-12)
            needed = math.ceil(
                (len(left) + len(right) - 2 * width + 2 - 2 * width * edits) / 2
            )
            a = Counter(left[i : i + width] for i in range(len(left) - width + 1))
            b = Counter(right[i : i + width] for i in range(len(right) - width + 1))
            if needed > 0 and sum((a & b).values()) < needed:
                continue
        total += 1
    return total


@pytest.mark.parametrize("threshold", [0.0, 0.7, 0.90, 0.99, 1.0])
@pytest.mark.parametrize("mixed_lengths", [False, True])
def test_sparse_pairs_preserve_graph_and_candidate_counts(
    threshold: float, mixed_lengths: bool
) -> None:
    rng = random.Random(20261102)
    lengths = [4, 12, 24, 60] if mixed_lengths else [48, 48, 48, 48]
    reads = [
        "".join(rng.choices("ACGT", k=rng.choice(lengths)))
        for _ in range(27)
    ]
    reads.extend([reads[1], reads[7], reads[9][:4] + "T" + reads[9][5:]])
    indexed = build_similarity_graph(reads, threshold=threshold, use_qgram_prefilter=True)
    exhaustive = build_similarity_graph(reads, threshold=threshold, use_qgram_prefilter=False)
    assert indexed.edges == exhaustive.edges
    assert indexed.candidate_pairs == reference_prefilter_count(reads, threshold, 5)


def test_sparse_pairs_skip_unrelated_sequences_without_scoring() -> None:
    class CountingScorer:
        def __init__(self) -> None:
            self.calls = 0

        def score(self, left: str, right: str) -> float:
            self.calls += 1
            return float(left == right)

    scorer = CountingScorer()
    graph = build_similarity_graph(
        ["A" * 48, "C" * 48, "G" * 48, "T" * 48],
        threshold=0.99,
        scorer=scorer,
    )
    assert graph.edges == ()
    assert graph.candidate_pairs == 0
    assert scorer.calls == 0

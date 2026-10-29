"""Parity tests for streamed graph scoring against the previous materialized path."""

from __future__ import annotations

import math
import random

import pytest

from oligoark.reconstruct import (
    GraphEdge,
    LevenshteinEdgeScorer,
    SimilarityGraph,
    _qgram_counts,
    build_similarity_graph,
)


def _materialized_reference(
    reads: list[str], threshold: float, width: int, prefilter: bool
) -> SimilarityGraph:
    """Independent copy of the old candidate-list algorithm for parity testing."""
    nodes = tuple(reads)
    all_pairs = (
        (left, right)
        for left in range(len(nodes))
        for right in range(left + 1, len(nodes))
    )
    if prefilter:
        signatures = [_qgram_counts(read, width) for read in nodes]
        postings: dict[str, list[tuple[int, int]]] = {}
        for index, signature in enumerate(signatures):
            for gram, count in signature.items():
                postings.setdefault(gram, []).append((index, count))
        shared: dict[tuple[int, int], int] = {}
        for entries in postings.values():
            for offset, (left, left_count) in enumerate(entries):
                for right, right_count in entries[offset + 1 :]:
                    pair = (left, right)
                    shared[pair] = shared.get(pair, 0) + min(left_count, right_count)
        candidates: list[tuple[int, int]] = []
        for left, right in all_pairs:
            left_len, right_len = len(nodes[left]), len(nodes[right])
            maximum = max(1, left_len, right_len)
            if 1.0 - abs(left_len - right_len) / maximum < threshold:
                continue
            if min(left_len, right_len) < width:
                candidates.append((left, right))
                continue
            max_edits = math.floor((1.0 - threshold) * maximum + 1e-12)
            required = math.ceil(
                (left_len - width + 1 + right_len - width + 1 - 2 * width * max_edits)
                / 2
            )
            if required <= 0 or shared.get((left, right), 0) >= required:
                candidates.append((left, right))
    else:
        candidates = list(all_pairs)
    scorer = LevenshteinEdgeScorer()
    edges = []
    for left, right in candidates:
        maximum = max(1, len(nodes[left]), len(nodes[right]))
        if 1.0 - abs(len(nodes[left]) - len(nodes[right])) / maximum < threshold:
            continue
        score = scorer.score(nodes[left], nodes[right])
        if score >= threshold:
            edges.append(GraphEdge(left, right, score))
    return SimilarityGraph(nodes, tuple(edges), candidate_pairs=len(candidates))


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("width", [1, 5])
@pytest.mark.parametrize("threshold", [0.0, 0.85, 1.0])
@pytest.mark.parametrize("prefilter", [False, True])
def test_streamed_graph_exactly_matches_materialized_reference(
    seed: int, width: int, threshold: float, prefilter: bool
) -> None:
    rng = random.Random(seed)
    reads = [
        "".join(rng.choices("ACGT", k=rng.randrange(0, 20)))
        for _ in range(9)
    ]
    reads += [reads[0], reads[1]]
    actual = build_similarity_graph(
        reads, threshold=threshold, qgram_width=width, use_qgram_prefilter=prefilter
    )
    assert actual == _materialized_reference(reads, threshold, width, prefilter)
    assert actual.connected_components() == _materialized_reference(
        reads, threshold, width, prefilter
    ).connected_components()


def test_streamed_exhaustive_mode_counts_pairs_even_when_lengths_reject() -> None:
    reads = ["", "A", "A" * 20, "T" * 20]
    actual = build_similarity_graph(reads, threshold=1.0, use_qgram_prefilter=False)
    assert actual.candidate_pairs == 6
    assert actual == _materialized_reference(reads, 1.0, 5, False)


def test_streamed_custom_scorer_keeps_left_right_order() -> None:
    class RecordingScorer:
        def __init__(self) -> None:
            self.calls: list[tuple[str, str]] = []

        def score(self, left: str, right: str) -> float:
            self.calls.append((left, right))
            return 1.0

    reads = ["ACGTACGT", "ACGTACGA", "ACGTTCGA", "TTTTTTTT"]
    scorer = RecordingScorer()
    result = build_similarity_graph(
        reads, threshold=0.8, qgram_width=3, scorer=scorer, use_qgram_prefilter=True
    )
    assert scorer.calls == [(reads[e.left], reads[e.right]) for e in result.edges]
    assert len(scorer.calls) == result.candidate_pairs

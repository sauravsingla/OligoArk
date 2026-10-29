"""Performance fast-path regressions; preserve decoding and ranking semantics."""

from __future__ import annotations

import random

import pytest

import oligoark.reconstruct as reconstruction
from oligoark.ecc import _generator_poly, rs_decode, rs_encode


def _reference_distance(a: str, b: str) -> int:
    previous = list(range(len(b) + 1))
    for i, left in enumerate(a, 1):
        current = [i]
        for j, right in enumerate(b, 1):
            current.append(
                min(current[-1] + 1, previous[j] + 1, previous[j - 1] + (left != right))
            )
        previous = current
    return previous[-1]


@pytest.mark.parametrize("seed", range(5))
def test_python_distance_exact_on_sparse_and_random_edits(monkeypatch: pytest.MonkeyPatch, seed: int) -> None:
    monkeypatch.setattr(reconstruction, "_native_edit_distance", None)
    rng = random.Random(seed)
    for _ in range(60):
        left = "".join(rng.choices("ACGT", k=rng.randrange(0, 45)))
        right = "".join(rng.choices("ACGT", k=rng.randrange(0, 45)))
        assert reconstruction.edit_distance(left, right) == _reference_distance(left, right)
        assert reconstruction.edit_distance(right, left) == _reference_distance(left, right)
        assert reconstruction.edit_distance(left, left) == 0
    original = "A" * 80 + "CGT" + "T" * 80
    assert reconstruction.edit_distance(original, "A" * 80 + "CAT" + "T" * 80) == 1


def test_optional_native_distance_exact(monkeypatch: pytest.MonkeyPatch) -> None:
    rapidfuzz = pytest.importorskip("rapidfuzz.distance")
    monkeypatch.setattr(reconstruction, "_native_edit_distance", rapidfuzz.Levenshtein.distance)
    rng = random.Random(41)
    for _ in range(120):
        a = "".join(rng.choices("ACGT", k=rng.randrange(0, 80)))
        b = "".join(rng.choices("ACGT", k=rng.randrange(0, 80)))
        assert reconstruction.edit_distance(a, b) == _reference_distance(a, b)


def test_weighted_medoid_and_anchor_order_match_exhaustive(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(reconstruction, "_native_edit_distance", None)
    rng = random.Random(17)
    for _ in range(30):
        unique = ["".join(rng.choices("ACGT", k=rng.randrange(2, 13))) for _ in range(6)]
        cluster = [rng.choice(unique) for _ in range(15)]
        unique_sorted = sorted(set(cluster))
        naive = {item: sum(_reference_distance(item, other) for other in cluster) for item in unique_sorted}
        assert reconstruction._weighted_read_costs(cluster) == naive
        expected_medoid = min(cluster, key=lambda x: (naive[x], len(x), x))
        assert reconstruction._medoid(cluster) == expected_medoid
        expected_anchors = sorted(
            unique_sorted,
            key=lambda x: (naive[x], abs(len(x) - 8), len(x), x),
        )
        assert reconstruction._trace_anchor_order(cluster, target_length=8) == expected_anchors


def test_qgram_reuse_matches_original() -> None:
    cluster = ["ACGTACGT", "ACGTACGA", "ACGTTACGT"]
    signatures = [reconstruction._qgrams(read, 4) for read in cluster]
    for candidate in ["ACGTACGT", "ACGTACGA", "TCGTACGT"]:
        original = reconstruction._candidate_qgram_similarity(candidate, cluster, width=4)
        accelerated = reconstruction._candidate_qgram_similarity(
            candidate, cluster, width=4, signatures=signatures
        )
        assert accelerated == original


def test_cached_reed_solomon_generator_and_recovery() -> None:
    _generator_poly.cache_clear()
    assert isinstance(_generator_poly(8), tuple)
    assert _generator_poly(8) == _generator_poly(8)
    assert _generator_poly.cache_info().hits >= 1
    for parity in (4, 8, 12):
        payload = bytes(range(64))
        packed = rs_encode(payload, parity)
        assert rs_decode(packed, parity) == payload
    _generator_poly.cache_clear()

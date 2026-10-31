"""Single-indel repair enumeration must preserve all prior CRC-gated candidates."""

from __future__ import annotations

import random
from collections.abc import Iterator

import pytest

from oligoark.dna import SequenceConstraints
from oligoark.framing import (
    _single_indel_candidates,
    decode_frame_resilient,
    encode_frame,
)


def _reference_candidates(sequence: str) -> list[str]:
    """Pre-optimization list-and-set implementation, including candidate order."""
    sequence = sequence.strip().upper()
    if not sequence or set(sequence) - {"A", "C", "G", "T"}:
        return []
    candidates: list[str] = []
    seen: set[str] = set()

    def add(candidate: str) -> None:
        if candidate not in seen:
            seen.add(candidate)
            candidates.append(candidate)

    if len(sequence) % 4 == 1:
        for position in range(len(sequence)):
            add(sequence[:position] + sequence[position + 1 :])
    elif len(sequence) % 4 == 3:
        for position in range(len(sequence) + 1):
            for base in "ACGT":
                add(sequence[:position] + base + sequence[position:])
    return candidates


@pytest.mark.parametrize("seed", range(5))
def test_streamed_candidates_match_previous_exact_order(seed: int) -> None:
    rng = random.Random(seed)
    for length in range(80):
        sequence = "".join(rng.choices("ACGT", k=length))
        assert list(_single_indel_candidates(sequence)) == _reference_candidates(sequence)
    for length in (599, 601):
        sequence = "".join(rng.choices("ACGT", k=length))
        assert list(_single_indel_candidates(sequence)) == _reference_candidates(sequence)
    for sequence in ("AAAAA", "AAAAAAA", "ACGNN", "ACGT", "", "acgta"):
        assert list(_single_indel_candidates(sequence)) == _reference_candidates(sequence)


def test_indel_candidates_are_lazy() -> None:
    candidates = _single_indel_candidates("ACGT" * 150 + "A")
    assert isinstance(candidates, Iterator)
    assert next(candidates) == ("ACGT" * 150 + "A")[1:]


@pytest.mark.parametrize("corruption", ("insert", "delete"))
def test_streamed_rescue_retains_exact_crc_gated_recovery(corruption: str) -> None:
    payload = bytes(range(16))
    sequence = encode_frame(
        payload,
        index=7,
        total_data=11,
        is_parity=False,
        rs_nsym=4,
        adaptive_masks=False,
        mask_search_limit=1,
        sequence_constraints=SequenceConstraints(
            min_gc_fraction=0.0, max_gc_fraction=1.0, max_homopolymer=1000
        ),
        compact_framing=True,
        inline_mask_framing=True,
    )
    damaged = (
        sequence[:45] + "A" + sequence[45:]
        if corruption == "insert"
        else sequence[:45] + sequence[46:]
    )
    decoded = decode_frame_resilient(
        damaged,
        rs_nsym=4,
        mask_search_limit=1,
        compact_framing=True,
        inline_mask_framing=True,
        expected_total_data=11,
    )
    assert decoded.index == 7
    assert decoded.payload == payload

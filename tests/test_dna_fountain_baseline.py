from __future__ import annotations

import os

from oligoark.baselines import (
    DnaFountainBaselineConfig,
    decode_dna_fountain_baseline,
    encode_dna_fountain_baseline,
)


def test_clean_room_dna_fountain_baseline_roundtrip() -> None:
    payload = os.urandom(512)
    archive = encode_dna_fountain_baseline(
        payload,
        DnaFountainBaselineConfig(redundancy=2.0),
    )

    recovered = decode_dna_fountain_baseline(archive)

    assert recovered == payload
    assert all(len(sequence) == 152 for sequence in archive.sequences)
    assert archive.logical_bits_per_nucleotide > 0


def test_dna_fountain_baseline_is_sha_verified() -> None:
    payload = b"dna fountain baseline" * 16
    archive = encode_dna_fountain_baseline(
        payload,
        DnaFountainBaselineConfig(redundancy=2.0),
    )

    assert decode_dna_fountain_baseline(archive) == payload
    assert len(archive.sha256) == 64

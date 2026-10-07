from __future__ import annotations

from oligoark.baselines import (
    DnaFountainBaselineConfig,
    RotatingTernaryBaselineConfig,
    decode_dna_fountain_baseline,
    decode_rotating_ternary_baseline,
    encode_dna_fountain_baseline,
    encode_rotating_ternary_baseline,
)


def test_clean_room_dna_fountain_baseline_roundtrip() -> None:
    # Keep CI reproducible: biochemical screening depends on payload bytes, so an
    # os.urandom fixture can occasionally select a rank-deficient droplet set.
    payload = b"dna fountain baseline" * 16
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


def test_rotating_ternary_baseline_roundtrip() -> None:
    payload = bytes(range(256)) * 4
    archive = encode_rotating_ternary_baseline(
        payload,
        RotatingTernaryBaselineConfig(redundancy=0.25),
    )

    recovered = decode_rotating_ternary_baseline(archive)

    assert recovered == payload
    assert max(map(len, archive.sequences)) <= 152
    assert all(
        left != right
        for sequence in archive.sequences
        for left, right in zip(sequence, sequence[1:], strict=False)
    )
    assert archive.logical_bits_per_nucleotide > 0


def test_rotating_ternary_xor_recovers_one_loss_per_group() -> None:
    payload = bytes(range(144))
    archive = encode_rotating_ternary_baseline(
        payload,
        RotatingTernaryBaselineConfig(redundancy=0.25),
    )
    # 12-byte chunks and 25% redundancy imply one parity strand per four data strands.
    kept = [
        sequence
        for index, sequence in enumerate(archive.sequences)
        if index not in {0, 4, 8}
    ]

    assert decode_rotating_ternary_baseline(archive, kept) == payload

from __future__ import annotations

import random

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


def test_dna_fountain_progress_is_diagnostic_only() -> None:
    payload = b"dna fountain baseline" * 16
    config = DnaFountainBaselineConfig(redundancy=2.0)
    reference = encode_dna_fountain_baseline(payload, config)
    observed: list[tuple[int, int, int]] = []
    instrumented = encode_dna_fountain_baseline(
        payload,
        config,
        progress=lambda accepted, attempts, target: observed.append(
            (accepted, attempts, target)
        ),
    )

    assert observed[0] == (0, 0, len(reference.sequences))
    assert observed[-1] == (
        len(instrumented.sequences),
        instrumented.attempts,
        len(instrumented.sequences),
    )
    assert instrumented == reference
    assert decode_dna_fountain_baseline(instrumented) == payload


def test_dna_fountain_seeded_randomized_progress_parity() -> None:
    rng = random.Random(20261101)
    config = DnaFountainBaselineConfig(redundancy=2.0)
    for size in (1, 33, 257, 1024):
        payload = rng.randbytes(size)
        reference = encode_dna_fountain_baseline(payload, config)
        instrumented = encode_dna_fountain_baseline(
            payload,
            config,
            progress=lambda _accepted, _attempts, _target: None,
        )
        assert instrumented == reference
        assert decode_dna_fountain_baseline(instrumented) == payload

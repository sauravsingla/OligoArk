"""Reference parity, error contract and exact DNA conversion integrity."""

from __future__ import annotations

import hashlib
import random

import pytest

from oligoark.dna import BASE_TO_BITS, bytes_to_dna, dna_to_bytes


def _reference_dna_to_bytes(sequence: str) -> bytes:
    """Pre-optimization decoder retained independently for parity tests."""
    sequence = sequence.strip().upper()
    if len(sequence) % 4:
        raise ValueError("DNA sequence length must be divisible by 4")
    try:
        vals = [BASE_TO_BITS[base] for base in sequence]
    except KeyError as exc:
        raise ValueError(f"Invalid DNA base: {exc.args[0]!r}") from exc
    return bytes(
        (vals[i] << 6) | (vals[i + 1] << 4) | (vals[i + 2] << 2) | vals[i + 3]
        for i in range(0, len(vals), 4)
    )


@pytest.mark.parametrize("seed", (0, 1, 7, 31, 20261031))
def test_direct_decode_exact_reference_parity(seed: int) -> None:
    rng = random.Random(seed)
    for _ in range(500):
        size = rng.randrange(0, 256)
        data = rng.randbytes(size)
        sequence = bytes_to_dna(data)
        assert dna_to_bytes(sequence) == _reference_dna_to_bytes(sequence) == data
        assert hashlib.sha256(dna_to_bytes(sequence)).digest() == hashlib.sha256(data).digest()
        assert dna_to_bytes("  " + sequence.lower() + "\n") == data


@pytest.mark.parametrize(
    "sequence",
    ("ACG", "ACGX", "XACG", "ACGTXACG", "ACG😎", "ACGTNACG", " ACGA N  "),
)
def test_direct_decode_preserves_invalid_input_contract(sequence: str) -> None:
    with pytest.raises(ValueError) as reference_error:
        _reference_dna_to_bytes(sequence)
    with pytest.raises(ValueError) as optimized_error:
        dna_to_bytes(sequence)
    assert str(optimized_error.value) == str(reference_error.value)


@pytest.mark.parametrize("size", (0, 1, 38, 150, 256, 1024))
def test_direct_decode_long_round_trip_sha256(size: int) -> None:
    data = random.Random(20261031 + size).randbytes(size)
    recovered = dna_to_bytes(bytes_to_dna(data))
    assert recovered == data
    assert hashlib.sha256(recovered).hexdigest() == hashlib.sha256(data).hexdigest()

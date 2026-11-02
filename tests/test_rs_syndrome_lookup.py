"""Parity and exact-recovery regression tests for cached RS syndrome multipliers."""

from __future__ import annotations

import random

import pytest

from oligoark.ecc import _gf_pow, _poly_eval, _syndromes, rs_decode, rs_encode


def _reference_syndromes(msg: bytes | bytearray, nsym: int) -> list[int]:
    """Previous evaluator, retained independently for byte-exact parity."""
    return [0] + [_poly_eval(msg, _gf_pow(2, i)) for i in range(nsym)]


@pytest.mark.parametrize("nsym", (0, 1, 2, 4, 8, 12, 16, 24, 32))
def test_randomized_syndrome_parity(nsym: int) -> None:
    rng = random.Random(20261102 + nsym)
    for width in (1, 2, 8, 32, 64, 120, 220, 255):
        for _ in range(25):
            raw = rng.randbytes(width)
            for msg in (raw, bytearray(raw)):
                assert _syndromes(msg, nsym) == _reference_syndromes(msg, nsym)


def test_empty_codeword_preserves_error_and_zero_parity_contract() -> None:
    assert _syndromes(b"", 0) == [0]
    with pytest.raises(IndexError):
        _syndromes(b"", 4)


@pytest.mark.parametrize("nsym", (2, 4, 8, 12, 16))
def test_randomized_reed_solomon_exact_recovery(nsym: int) -> None:
    rng = random.Random(120026 + nsym)
    for width in (8, 32, 80, 155):
        for _ in range(12):
            payload = rng.randbytes(width)
            codeword = rs_encode(payload, nsym)
            assert rs_decode(codeword, nsym) == payload
            damaged = bytearray(codeword)
            for position in rng.sample(range(len(damaged)), rng.randrange(nsym // 2 + 1)):
                damaged[position] ^= rng.randrange(1, 256)
            assert rs_decode(bytes(damaged), nsym) == payload

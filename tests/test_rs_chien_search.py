"""Parity and exact-recovery checks for the specialized Reed-Solomon Chien search."""

from __future__ import annotations

import random

import pytest

from oligoark.ecc import (
    ECCDecodeError,
    _find_errors,
    _gf_pow,
    _poly_eval,
    _poly_mul,
    rs_decode,
    rs_encode,
)


def _reference_find_errors(err_loc: list[int], message_len: int) -> list[int]:
    """Pre-optimization generic field-evaluation implementation."""
    errs = len(err_loc) - 1
    positions: list[int] = []
    for i in range(message_len):
        if _poly_eval(err_loc, _gf_pow(2, i)) == 0:
            positions.append(message_len - 1 - i)
    if len(positions) != errs:
        raise ECCDecodeError("could not locate all Reed-Solomon symbol errors")
    return positions


@pytest.mark.parametrize("seed", range(5))
def test_chien_search_matches_reference_for_random_polynomials(seed: int) -> None:
    rng = random.Random(seed)
    for message_len in (1, 16, 80, 255):
        for degree in (0, 1, 2, 4, 8, 16):
            for _ in range(12):
                locator = [rng.randrange(1, 256)]
                locator.extend(rng.randrange(256) for _ in range(degree))
                try:
                    expected = _reference_find_errors(locator, message_len)
                except ECCDecodeError:
                    with pytest.raises(ECCDecodeError):
                        _find_errors(locator, message_len)
                else:
                    assert _find_errors(locator, message_len) == expected


@pytest.mark.parametrize("seed", range(5))
def test_chien_search_locates_known_roots_in_identical_order(seed: int) -> None:
    rng = random.Random(seed)
    for message_len in (32, 80, 255):
        for degree in (1, 2, 4, 8, 16):
            indices = rng.sample(range(message_len), min(degree, message_len))
            locator = [1]
            for index in indices:
                locator = _poly_mul(locator, [1, _gf_pow(2, index)])
            expected = sorted((message_len - 1 - i for i in indices), reverse=True)
            assert _find_errors(locator, message_len) == expected
            assert _find_errors(locator, message_len) == _reference_find_errors(
                locator, message_len
            )


@pytest.mark.parametrize("seed", range(4))
def test_rs_exact_roundtrip_with_correctable_random_symbol_errors(seed: int) -> None:
    rng = random.Random(seed)
    for nsym in (4, 8, 16):
        for size in (1, 16, 80, 200):
            original = rng.randbytes(size)
            codeword = rs_encode(original, nsym)
            for errors in (0, 1, nsym // 2):
                corrupted = bytearray(codeword)
                for position in rng.sample(range(len(corrupted)), errors):
                    corrupted[position] ^= rng.randrange(1, 256)
                assert rs_decode(bytes(corrupted), nsym) == original

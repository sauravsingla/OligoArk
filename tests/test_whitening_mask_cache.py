"""The cached mask path must preserve every original whitening byte and archive digest."""

from __future__ import annotations

import random

import pytest

from oligoark.archive import ArchiveConfig, archive_bytes, recover_bytes
from oligoark.framing import _mask, _mask_pad_int


def _reference_mask(data: bytes, mask_id: int) -> bytes:
    """Frozen pre-optimization implementation; independent parity oracle."""
    if not 0 <= mask_id <= 255:
        raise ValueError("mask_id must fit in one byte")
    legacy = (0, 0x55, 0xAA, 0xFF)
    if mask_id < len(legacy):
        value = legacy[mask_id]
        if value == 0:
            return data
        return bytes(byte ^ value for byte in data)

    state = (0x9E3779B9 ^ (mask_id * 0x45D9F3B)) & 0xFFFFFFFF
    out = bytearray(len(data))
    for index, byte in enumerate(data):
        state ^= (state << 13) & 0xFFFFFFFF
        state ^= state >> 17
        state ^= (state << 5) & 0xFFFFFFFF
        out[index] = byte ^ (state & 0xFF)
    return bytes(out)


@pytest.mark.parametrize("length", (0, 1, 2, 3, 4, 35, 36, 37, 96, 152, 255, 256, 320))
def test_all_mask_ids_match_original_and_round_trip(length: int) -> None:
    rng = random.Random(20261102 + length)
    for mask_id in range(256):
        for _ in range(2):
            data = rng.randbytes(length)
            expected = _reference_mask(data, mask_id)
            assert _mask(data, mask_id) == expected
            assert _mask(expected, mask_id) == data


def test_mask_cache_bounded_and_not_used_for_oversize_inputs() -> None:
    _mask_pad_int.cache_clear()
    for length in (256, 257, 512, 1024):
        data = bytes(range(256)) * (length // 256) + bytes(range(length % 256))
        assert _mask(data, 255) == _reference_mask(data, 255)
    assert _mask_pad_int.cache_info().currsize == 0
    for mask_id in range(1, 256):
        for length in (4, 36, 96):
            assert _mask(b"A" * length, mask_id) == _reference_mask(
                b"A" * length, mask_id
            )
    assert _mask_pad_int.cache_info().currsize <= 1024


def test_mask_rejects_out_of_range_ids_without_changing_contract() -> None:
    for bad_id in (-1, 256, 65536):
        with pytest.raises(ValueError, match="mask_id"):
            _mask(b"ABC", bad_id)


def test_archive_round_trip_after_mask_change() -> None:
    payload = bytes(range(256)) * 6
    config = ArchiveConfig(
        chunk_size=24,
        rs_nsym=4,
        parity_group_size=3,
        compact_framing=True,
        inline_mask_framing=True,
        mask_search_limit=256,
        min_gc_fraction=0.0,
        max_gc_fraction=1.0,
        max_homopolymer=1024,
    )
    archive = archive_bytes(payload, config)
    assert recover_bytes(archive, archive.strands) == payload

"""Regressions for backward-compatible extended compact-inline whitening masks."""

import pytest

from oligoark.dna import SequenceConstraintError, SequenceConstraints, bytes_to_dna
from oligoark.framing import (
    _inline_decode_mask_candidates,
    _inline_mask_from_selector,
    _inline_selector_for_mask,
    compact_inline_frame_hint,
    decode_frame_packed,
    encode_frame_packed,
)
from oligoark.profiles import physical_strand_profile


def test_extended_inline_selector_is_bijective_and_legacy_stable():
    selectors = [_inline_selector_for_mask(mask_id) for mask_id in range(256)]
    assert len(set(selectors)) == 256
    assert selectors[:16] == list(range(0xB0, 0xC0))
    for mask_id, selector in enumerate(selectors):
        assert _inline_mask_from_selector(selector) == mask_id


def test_compact_inline_can_recover_when_all_legacy_masks_are_rejected(monkeypatch):
    # Force a rare mask-exhaustion scenario deterministically; decoding is still
    # guarded by RS and CRC16, not by the mocked sequence-constraint decision.
    def accept_only_extended(self, sequence):
        return sequence.startswith("AAAA")  # mask 16 has selector byte 0x00

    monkeypatch.setattr(SequenceConstraints, "accepts", accept_only_extended)
    kwargs = {
        "index": 43,
        "total_data": 1,
        "is_parity": False,
        "rs_nsym": 2,
        "adaptive_masks": True,
        "compact_framing": True,
        "inline_mask_framing": True,
    }
    payload = bytes(range(30))
    with pytest.raises(SequenceConstraintError):
        encode_frame_packed(payload, mask_search_limit=16, **kwargs)

    packed = encode_frame_packed(payload, mask_search_limit=256, **kwargs)
    assert packed[0] == 0x00
    assert len(packed) * 4 <= 152
    frame = decode_frame_packed(
        packed,
        rs_nsym=2,
        compact_framing=True,
        inline_mask_framing=True,
        expected_total_data=1,
        mask_search_limit=256,
    )
    assert frame.payload == payload
    assert frame.index == 43

    # A prefix-only identity hint must support the extended mask selector too.
    # It is not accepted as recovered data without CRC/RS verification.
    hint = compact_inline_frame_hint(bytes_to_dna(packed))
    assert hint is not None
    assert hint.index == 43
    assert hint.is_parity is False
    assert hint.is_fountain is False
    assert compact_inline_frame_hint(bytes_to_dna(packed)[:-1]) == hint

    # Damage one of the selector's four DNA bases. The bounded fallback must
    # recover the extended mask without scanning all 256 Reed-Solomon variants.
    changed_selector = bytes([packed[0] ^ 0x40]) + packed[1:]
    candidates = _inline_decode_mask_candidates(changed_selector[0], 256)
    assert 16 in candidates
    assert len(candidates) <= 29
    repaired = decode_frame_packed(
        changed_selector,
        rs_nsym=2,
        compact_framing=True,
        inline_mask_framing=True,
        expected_total_data=1,
        mask_search_limit=256,
    )
    assert repaired.payload == payload
    assert repaired.index == 43


def test_v3_profile_can_search_all_inline_selectors():
    profile = physical_strand_profile("oligoark-152-compact-v3")
    assert profile.mask_search_limit == 256
    assert profile.actual_nucleotides == 152

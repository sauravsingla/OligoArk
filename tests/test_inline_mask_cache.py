"""Cache parity and immutability for deterministic inline selector repairs."""

from oligoark.framing import _inline_decode_mask_candidates, _inline_mask_from_selector


def _reference_candidates(selector: int, limit: int) -> tuple[int, ...]:
    result: list[int] = []
    seen: set[int] = set()

    def add(mask: int) -> None:
        if mask < limit and mask not in seen:
            seen.add(mask)
            result.append(mask)

    add(_inline_mask_from_selector(selector))
    for mask in range(min(limit, 16)):
        add(mask)
    for shift in (6, 4, 2, 0):
        base = (selector >> shift) & 3
        for new_base in range(4):
            if new_base != base:
                repaired = (selector & ~(3 << shift)) | (new_base << shift)
                add(_inline_mask_from_selector(repaired))
    return tuple(result)


def test_all_selectors_and_representative_limits_match_reference() -> None:
    for limit in (1, 4, 16, 17, 64, 128, 256):
        for selector in range(256):
            assert _inline_decode_mask_candidates(selector, limit) == _reference_candidates(
                selector, limit
            )


def test_cached_candidates_are_immutable_and_bounded() -> None:
    candidates = _inline_decode_mask_candidates(0xB2, 256)
    assert isinstance(candidates, tuple)
    assert candidates is _inline_decode_mask_candidates(0xB2, 256)
    assert _inline_decode_mask_candidates.cache_info().maxsize == 1024

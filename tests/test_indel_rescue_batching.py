"""Indel-rescue batching must preserve exact archive recovery and bound full peels."""

from __future__ import annotations

import hashlib

import pytest

import oligoark.archive as archive_module
from oligoark.archive import ArchiveConfig, archive_bytes
from oligoark.framing import CompactFrameHint, DecodedFrame


@pytest.mark.parametrize(
    ("total_data", "expected_peels"),
    [(64, 3), (512, 5)],
)
def test_indel_rescue_batches_full_redundancy_passes(
    monkeypatch: pytest.MonkeyPatch,
    total_data: int,
    expected_peels: int,
) -> None:
    """Synthetic valid frames isolate the batching schedule from RS search cost."""
    width = 64
    payload = bytes(range(256)) * (total_data * width // 256)
    config = ArchiveConfig(
        chunk_size=width,
        rs_nsym=0,
        adaptive_masks=False,
        redundancy_scheme="none",
        indel_rescue=True,
        min_gc_fraction=0.0,
        max_gc_fraction=1.0,
        max_homopolymer=1024,
    )
    archive = archive_bytes(payload, config)
    missing = set(range(0, total_data, 4))

    def decoded_frames(*_args: object) -> tuple[object, object, object, object]:
        known = {
            index: payload[index * width : (index + 1) * width]
            for index in range(total_data)
            if index not in missing
        }
        hints = [
            (
                f"candidate-{index}",
                CompactFrameHint(index=index, is_parity=False, is_fountain=False),
            )
            for index in sorted(missing)
        ]
        return known, {}, [], hints

    def resilient(sequence: str, **_kwargs: object) -> DecodedFrame:
        index = int(sequence.split("-")[1])
        return DecodedFrame(
            index=index,
            total_data=total_data,
            payload=payload[index * width : (index + 1) * width],
            is_parity=False,
        )

    original_peel = archive_module._apply_redundancy
    peel_calls = 0

    def counted_peel(*args: object, **kwargs: object) -> dict[int, bytes]:
        nonlocal peel_calls
        peel_calls += 1
        return original_peel(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(archive_module, "_decode_available_frames", decoded_frames)
    monkeypatch.setattr(archive_module, "decode_frame_resilient", resilient)
    monkeypatch.setattr(archive_module, "_apply_redundancy", counted_peel)

    recovered = archive_module.recover_bytes(archive, [])
    assert recovered == payload
    assert hashlib.sha256(recovered).hexdigest() == archive.metadata["sha256"]
    assert peel_calls == expected_peels

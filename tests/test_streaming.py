from __future__ import annotations

import os

import pytest

from oligoark.archive import ArchiveConfig
from oligoark.framing import decode_frame, decode_frame_packed, encode_frame, encode_frame_packed
from oligoark.profiles import physical_strand_profile
from oligoark.streaming import (
    StreamingFaultProfile,
    archive_file_streaming,
    recover_file_streaming,
)


def _scale_config(scheme: str = "xor") -> ArchiveConfig:
    return ArchiveConfig(
        chunk_size=237,
        rs_nsym=0,
        parity_group_size=8,
        adaptive_masks=False,
        redundancy_scheme=scheme,
        fountain_redundancy=0.25,
        min_gc_fraction=0.0,
        max_gc_fraction=1.0,
        max_homopolymer=1024,
        mask_search_limit=1,
    )


def test_packed_frame_is_compatible_with_text_frame() -> None:
    payload = b"packed compatibility"
    config = _scale_config("none")
    packed = encode_frame_packed(
        payload,
        index=3,
        total_data=11,
        is_parity=False,
        rs_nsym=config.rs_nsym,
        adaptive_masks=config.adaptive_masks,
        sequence_constraints=config.sequence_constraints,
        mask_search_limit=config.mask_search_limit,
    )
    sequence = encode_frame(
        payload,
        index=3,
        total_data=11,
        is_parity=False,
        rs_nsym=config.rs_nsym,
        adaptive_masks=config.adaptive_masks,
        sequence_constraints=config.sequence_constraints,
        mask_search_limit=config.mask_search_limit,
    )

    packed_frame = decode_frame_packed(
        packed,
        rs_nsym=config.rs_nsym,
        mask_search_limit=config.mask_search_limit,
    )
    text_frame = decode_frame(
        sequence,
        rs_nsym=config.rs_nsym,
        mask_search_limit=config.mask_search_limit,
    )
    assert packed_frame == text_frame
    assert packed_frame.payload == payload


@pytest.mark.parametrize(
    ("name", "target"),
    [("oligoark-152", 152), ("oligoark-200", 200), ("oligoark-248", 248)],
)
def test_physical_profiles_fit_requested_oligo_length(name: str, target: int) -> None:
    profile = physical_strand_profile(name)
    config = profile.to_archive_config()
    packed = encode_frame_packed(
        b"x" * config.chunk_size,
        index=0,
        total_data=1,
        is_parity=False,
        rs_nsym=config.rs_nsym,
        adaptive_masks=config.adaptive_masks,
        sequence_constraints=config.sequence_constraints,
        mask_search_limit=config.mask_search_limit,
    )
    assert len(packed) * 4 == target
    assert profile.actual_nucleotides == target


@pytest.mark.parametrize("scheme", ["none", "xor", "fountain", "hybrid"])
def test_streaming_roundtrip_all_redundancy_modes(tmp_path, scheme: str) -> None:
    payload = os.urandom(4097)
    source = tmp_path / "source.bin"
    archive = tmp_path / f"{scheme}.oab"
    output = tmp_path / "output.bin"
    source.write_bytes(payload)

    stats = archive_file_streaming(source, archive, _scale_config(scheme))
    report = recover_file_streaming(archive, output)

    assert output.read_bytes() == payload
    assert report.verified_sha256 is True
    assert stats.source_sha256 == report.output_sha256
    assert stats.archive_size_bytes == archive.stat().st_size
    assert stats.encoded_nucleotides > len(payload) * 4


def test_controlled_dropout_is_repaired_without_losing_sha_integrity(tmp_path) -> None:
    payload = os.urandom(16_384)
    source = tmp_path / "source.bin"
    archive = tmp_path / "archive.oab"
    output = tmp_path / "output.bin"
    source.write_bytes(payload)

    archive_file_streaming(source, archive, _scale_config("xor"))
    report = recover_file_streaming(
        archive,
        output,
        fault=StreamingFaultProfile(dropout_rate=0.125, seed=2026),
    )

    assert report.verified_sha256 is True
    assert report.dropped_records > 0
    assert report.xor_recovered_strands == report.dropped_records
    assert output.read_bytes() == payload


def test_unprotected_streaming_archive_reports_controlled_loss(tmp_path) -> None:
    source = tmp_path / "source.bin"
    archive = tmp_path / "archive.oab"
    output = tmp_path / "output.bin"
    source.write_bytes(os.urandom(8192))
    archive_file_streaming(source, archive, _scale_config("none"))

    report = recover_file_streaming(
        archive,
        output,
        fault=StreamingFaultProfile(dropout_rate=0.125, seed=2026),
        strict=False,
    )

    assert report.verified_sha256 is False
    assert report.missing_data_strands > 0

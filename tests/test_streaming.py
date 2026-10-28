from __future__ import annotations

import json
import os

import pytest

from oligoark.archive import ArchiveConfig, archive_bytes, recover_bytes
from oligoark.framing import (
    compact_inline_frame_hint,
    decode_frame,
    decode_frame_packed,
    decode_frame_resilient,
    encode_frame,
    encode_frame_packed,
)
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
    [
        ("oligoark-152", 152),
        ("oligoark-152-compact-v3", 152),
        ("oligoark-200", 200),
        ("oligoark-248", 248),
    ],
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
        compact_framing=config.compact_framing,
        compact_index_bytes=config.compact_index_bytes,
        inline_mask_framing=config.inline_mask_framing,
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


def test_controlled_dropout_can_recover_short_final_chunk(tmp_path) -> None:
    # Nine chunks makes index 8 the first member of the final XOR group. The controlled
    # dropout profile therefore erases that short final chunk as well as index 0.
    payload = os.urandom(8 * 237 + 17)
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
    assert output.stat().st_size == len(payload)
    assert output.read_bytes() == payload


def test_fountain_dropout_recovery_uses_on_disk_rescan(tmp_path) -> None:
    payload = bytes(range(256)) * 64
    source = tmp_path / "source.bin"
    archive = tmp_path / "fountain.oab"
    output = tmp_path / "output.bin"
    source.write_bytes(payload)

    config = ArchiveConfig(
        chunk_size=237,
        rs_nsym=0,
        parity_group_size=8,
        adaptive_masks=False,
        redundancy_scheme="fountain",
        fountain_redundancy=0.75,
        min_gc_fraction=0.0,
        max_gc_fraction=1.0,
        max_homopolymer=1024,
        mask_search_limit=1,
    )
    archive_file_streaming(source, archive, config)
    report = recover_file_streaming(
        archive,
        output,
        fault=StreamingFaultProfile(dropout_rate=0.05, seed=2026),
    )

    assert report.verified_sha256 is True
    assert report.fountain_recovered_strands > 0
    assert output.read_bytes() == payload

    # A rescan must see the same noisy channel realization, never pristine fountain
    # records. With no inner RS protection this substitution regime remains unrecoverable.
    noisy = recover_file_streaming(
        archive,
        output,
        fault=StreamingFaultProfile(substitution_rate=0.001, seed=2026),
        strict=False,
    )
    assert noisy.verified_sha256 is False
    assert noisy.missing_data_strands > 0



def test_compact_frame_roundtrip_and_density_profile() -> None:
    profile = physical_strand_profile("oligoark-152-compact")
    config = profile.to_archive_config()
    payload = bytes(range(config.chunk_size))
    packed = encode_frame_packed(
        payload,
        index=1234,
        total_data=5000,
        is_parity=False,
        rs_nsym=config.rs_nsym,
        adaptive_masks=config.adaptive_masks,
        sequence_constraints=config.sequence_constraints,
        mask_search_limit=config.mask_search_limit,
        compact_framing=config.compact_framing,
        compact_index_bytes=config.compact_index_bytes,
    )
    decoded = decode_frame_packed(
        packed,
        rs_nsym=config.rs_nsym,
        mask_search_limit=config.mask_search_limit,
        compact_framing=config.compact_framing,
        compact_index_bytes=config.compact_index_bytes,
        expected_total_data=5000,
    )
    assert decoded.index == 1234
    assert decoded.payload == payload
    assert len(packed) * 4 == 152
    assert config.chunk_size >= 29


@pytest.mark.parametrize("kind", ["insertion", "deletion"])
def test_single_indel_rescue_recovers_compact_frame(kind: str) -> None:
    profile = physical_strand_profile("oligoark-248-compact")
    config = profile.to_archive_config()
    payload = bytes(range(config.chunk_size))
    packed = encode_frame_packed(
        payload,
        index=77,
        total_data=100,
        is_parity=False,
        rs_nsym=config.rs_nsym,
        adaptive_masks=config.adaptive_masks,
        sequence_constraints=config.sequence_constraints,
        mask_search_limit=config.mask_search_limit,
        compact_framing=config.compact_framing,
        compact_index_bytes=config.compact_index_bytes,
    )
    sequence = "".join("ACGT"[(byte >> shift) & 3] for byte in packed for shift in (6, 4, 2, 0))
    position = len(sequence) // 2
    if kind == "insertion":
        mutated = sequence[:position] + "A" + sequence[position:]
    else:
        mutated = sequence[:position] + sequence[position + 1 :]
    decoded = decode_frame_resilient(
        mutated,
        rs_nsym=config.rs_nsym,
        mask_search_limit=config.mask_search_limit,
        compact_framing=config.compact_framing,
        compact_index_bytes=config.compact_index_bytes,
        expected_total_data=100,
    )
    assert decoded.index == 77
    assert decoded.payload == payload


def test_inline_mask_152_profile_saves_one_byte_without_dropping_protection() -> None:
    previous = physical_strand_profile("oligoark-152-compact")
    profile = physical_strand_profile("oligoark-152-compact-v3")
    config = profile.to_archive_config()
    assert previous.chunk_size == 29
    assert profile.chunk_size == 30
    assert config.rs_nsym == 2
    assert config.inline_mask_framing is True

    payload = bytes(range(config.chunk_size))
    packed = encode_frame_packed(
        payload,
        index=123_456,
        total_data=300_000,
        is_parity=False,
        rs_nsym=config.rs_nsym,
        adaptive_masks=config.adaptive_masks,
        sequence_constraints=config.sequence_constraints,
        mask_search_limit=config.mask_search_limit,
        compact_framing=config.compact_framing,
        compact_index_bytes=config.compact_index_bytes,
        inline_mask_framing=config.inline_mask_framing,
    )
    decoded = decode_frame_packed(
        packed,
        rs_nsym=config.rs_nsym,
        mask_search_limit=config.mask_search_limit,
        compact_framing=config.compact_framing,
        compact_index_bytes=config.compact_index_bytes,
        expected_total_data=300_000,
        inline_mask_framing=config.inline_mask_framing,
    )
    assert decoded.index == 123_456
    assert decoded.payload == payload
    assert len(packed) * 4 == 152


@pytest.mark.parametrize("kind", ["insertion", "deletion"])
def test_inline_mask_152_profile_keeps_single_indel_rescue(kind: str) -> None:
    profile = physical_strand_profile("oligoark-152-compact-v3")
    config = profile.to_archive_config()
    payload = bytes(range(config.chunk_size))
    packed = encode_frame_packed(
        payload,
        index=77,
        total_data=100,
        is_parity=False,
        rs_nsym=config.rs_nsym,
        adaptive_masks=config.adaptive_masks,
        sequence_constraints=config.sequence_constraints,
        mask_search_limit=config.mask_search_limit,
        compact_framing=config.compact_framing,
        compact_index_bytes=config.compact_index_bytes,
        inline_mask_framing=config.inline_mask_framing,
    )
    sequence = "".join("ACGT"[(byte >> shift) & 3] for byte in packed for shift in (6, 4, 2, 0))
    position = len(sequence) // 2
    mutated = (
        sequence[:position] + "A" + sequence[position:]
        if kind == "insertion"
        else sequence[:position] + sequence[position + 1 :]
    )
    decoded = decode_frame_resilient(
        mutated,
        rs_nsym=config.rs_nsym,
        mask_search_limit=config.mask_search_limit,
        compact_framing=config.compact_framing,
        compact_index_bytes=config.compact_index_bytes,
        expected_total_data=100,
        inline_mask_framing=config.inline_mask_framing,
    )
    assert decoded.index == 77
    assert decoded.payload == payload


def test_inline_mask_hint_survives_late_single_indel() -> None:
    profile = physical_strand_profile("oligoark-152-compact-v3")
    config = profile.to_archive_config()
    payload = bytes(range(config.chunk_size))
    packed = encode_frame_packed(
        payload,
        index=123_456,
        total_data=300_000,
        is_parity=False,
        rs_nsym=config.rs_nsym,
        adaptive_masks=config.adaptive_masks,
        sequence_constraints=config.sequence_constraints,
        mask_search_limit=config.mask_search_limit,
        compact_framing=config.compact_framing,
        compact_index_bytes=config.compact_index_bytes,
        inline_mask_framing=config.inline_mask_framing,
    )
    sequence = "".join(
        "ACGT"[(byte >> shift) & 3]
        for byte in packed
        for shift in (6, 4, 2, 0)
    )
    position = 80
    inserted = sequence[:position] + "A" + sequence[position:]
    deleted = sequence[:position] + sequence[position + 1 :]
    for mutated in (inserted, deleted):
        hint = compact_inline_frame_hint(
            mutated,
            compact_index_bytes=config.compact_index_bytes,
        )
        assert hint is not None
        assert hint.index == 123_456
        assert hint.is_parity is False
        assert hint.is_fountain is False


def test_v3_archive_recovers_sparse_single_indels_without_bruteforce_scan() -> None:
    profile = physical_strand_profile("oligoark-152-compact-v3").with_scheme("hybrid")
    config = profile.to_archive_config()
    payload = bytes(range(256)) * 64
    archive = archive_bytes(payload, config)
    reads = list(archive.strands)
    damaged = max(1, round(len(reads) * 0.05))
    step = max(1, len(reads) // damaged)
    changed = 0
    for ordinal in range(0, len(reads), step):
        if changed >= damaged:
            break
        sequence = reads[ordinal]
        position = 40 + (ordinal % 80)
        if changed % 2:
            reads[ordinal] = sequence[:position] + "A" + sequence[position:]
        else:
            reads[ordinal] = sequence[:position] + sequence[position + 1 :]
        changed += 1
    assert changed == damaged
    assert recover_bytes(archive, reads) == payload
@pytest.mark.parametrize("alias_kind", ["same", "symlink", "hardlink"])
def test_streaming_encode_never_overwrites_its_input(tmp_path, alias_kind: str) -> None:
    source = tmp_path / "payload.bin"
    payload = os.urandom(256)
    source.write_bytes(payload)
    destination = source
    if alias_kind != "same":
        destination = tmp_path / "alias.bin"
        if alias_kind == "symlink":
            destination.symlink_to(source)
        else:
            os.link(source, destination)

    with pytest.raises(ValueError, match="distinct files"):
        archive_file_streaming(source, destination, _scale_config("none"))
    assert source.read_bytes() == payload


@pytest.mark.parametrize("alias_kind", ["same", "symlink", "hardlink"])
def test_streaming_recover_never_overwrites_its_archive(tmp_path, alias_kind: str) -> None:
    source = tmp_path / "payload.bin"
    source.write_bytes(b"immutable source archive" * 10)
    archive = tmp_path / "archive.oab"
    archive_file_streaming(source, archive, _scale_config("none"))
    original_archive = archive.read_bytes()
    destination = archive
    if alias_kind != "same":
        destination = tmp_path / "alias.oab"
        if alias_kind == "symlink":
            destination.symlink_to(archive)
        else:
            os.link(archive, destination)

    with pytest.raises(ValueError, match="distinct files"):
        recover_file_streaming(archive, destination)
    assert archive.read_bytes() == original_archive


def test_failed_strict_recovery_keeps_existing_output(tmp_path) -> None:
    source = tmp_path / "payload.bin"
    source.write_bytes(os.urandom(4096))
    archive = tmp_path / "archive.oab"
    output = tmp_path / "existing.bin"
    output.write_bytes(b"do not replace me")
    archive_file_streaming(source, archive, _scale_config("none"))

    with pytest.raises(ValueError, match="SHA-256"):
        recover_file_streaming(
            archive,
            output,
            fault=StreamingFaultProfile(dropout_rate=0.125, seed=2026),
        )
    assert output.read_bytes() == b"do not replace me"
    assert not list(tmp_path.glob(".existing.bin.*.part"))


def test_streaming_rejects_oversized_metadata_before_allocation(tmp_path) -> None:
    archive = tmp_path / "malicious.oab"
    archive.write_bytes(b"OAB2" + (1024 * 1024 + 1).to_bytes(4, "big"))
    with pytest.raises(ValueError, match="metadata exceeds"):
        recover_file_streaming(archive, tmp_path / "output.bin")
    assert not (tmp_path / "output.bin").exists()

def test_streaming_rejects_impossible_claimed_original_size(tmp_path) -> None:
    source = tmp_path / "payload.bin"
    source.write_bytes(b"short payload")
    archive = tmp_path / "archive.oab"
    archive_file_streaming(source, archive, _scale_config("none"))
    contents = archive.read_bytes()
    length = int.from_bytes(contents[4:8], "big")
    metadata = json.loads(contents[8 : 8 + length])
    metadata["original_size"] = 10**12
    modified = json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode()
    archive.write_bytes(
        contents[:4] + len(modified).to_bytes(4, "big") + modified + contents[8 + length :]
    )
    with pytest.raises(ValueError, match="invalid streaming strand counts"):
        recover_file_streaming(archive, tmp_path / "output.bin")
    assert not (tmp_path / "output.bin").exists()

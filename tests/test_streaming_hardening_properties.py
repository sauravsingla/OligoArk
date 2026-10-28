"""Fuzz recovery safety: damaged archives must not replace existing output."""

import json

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from oligoark.archive import ArchiveConfig
from oligoark.streaming import archive_file_streaming, recover_file_streaming


@given(
    field=st.sampled_from(["data_strands", "parity_strands", "fountain_strands", "strand_count"]),
    delta=st.integers(-100, 100).filter(lambda x: x != 0),
)
@settings(
    max_examples=32,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
def test_mutated_metadata_never_replaces_output(tmp_path, field, delta):
    source, archive, target = (
        tmp_path / "source.bin",
        tmp_path / "archive.oab",
        tmp_path / "existing.bin",
    )
    source.write_bytes(b"data" * 100)
    target.write_bytes(b"existing verified content")
    archive_file_streaming(source, archive, ArchiveConfig())
    data = archive.read_bytes()
    meta_size = int.from_bytes(data[4:8], "big")
    metadata = json.loads(data[8 : 8 + meta_size])
    metadata[field] += delta
    updated = json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode()
    archive.write_bytes(
        data[:4] + len(updated).to_bytes(4, "big") + updated + data[8 + meta_size :]
    )
    with pytest.raises(ValueError):
        recover_file_streaming(archive, target)
    assert target.read_bytes() == b"existing verified content"
    assert not list(tmp_path.glob(".existing.bin.*.part"))



@given(bytes_removed=st.integers(min_value=1, max_value=60))
@settings(
    max_examples=16,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
def test_truncated_archive_never_replaces_output(tmp_path, bytes_removed):
    source, archive, target = (
        tmp_path / "source.bin",
        tmp_path / "archive.oab",
        tmp_path / "existing.bin",
    )
    source.write_bytes(b"data" * 100)
    target.write_bytes(b"existing verified content")
    archive_file_streaming(source, archive, ArchiveConfig())
    data = archive.read_bytes()
    archive.write_bytes(data[:-bytes_removed])
    with pytest.raises(ValueError):
        recover_file_streaming(archive, target)
    assert target.read_bytes() == b"existing verified content"
    assert not list(tmp_path.glob(".existing.bin.*.part"))

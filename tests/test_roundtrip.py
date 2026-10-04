import os

import pytest

from oligoark.archive import ArchiveConfig, archive_bytes, recover_bytes
from oligoark.dna import bytes_to_dna, dna_to_bytes, sequence_metrics


def test_binary_dna_roundtrip() -> None:
    data = bytes(range(256))
    assert dna_to_bytes(bytes_to_dna(data)) == data


def test_archive_roundtrip() -> None:
    raw = os.urandom(1000)
    archive = archive_bytes(raw, ArchiveConfig(chunk_size=64, rs_nsym=8, parity_group_size=4))
    assert recover_bytes(archive) == raw


def test_single_dropout_per_parity_group_recovers() -> None:
    raw = bytes(range(200))
    archive = archive_bytes(raw, ArchiveConfig(chunk_size=32, rs_nsym=8, parity_group_size=3))
    reads = [strand for i, strand in enumerate(archive.strands) if i != 1]
    assert recover_bytes(archive, reads) == raw


def test_unrecoverable_dropout_raises() -> None:
    raw = bytes(range(200))
    archive = archive_bytes(raw, ArchiveConfig(chunk_size=32, rs_nsym=8, parity_group_size=3))
    reads = [strand for i, strand in enumerate(archive.strands) if i not in {0, 1}]
    with pytest.raises(ValueError, match="not recoverable"):
        recover_bytes(archive, reads)


def test_sequence_metrics() -> None:
    m = sequence_metrics("AACCGGTTTT")
    assert m.gc_fraction == pytest.approx(0.4)
    assert m.max_homopolymer == 4

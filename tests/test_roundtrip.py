import os

import pytest

from oligoark.archive import (
    ArchiveConfig,
    archive_bytes,
    recover_bytes,
    recover_from_reads,
)
from oligoark.dna import bytes_to_dna, dna_to_bytes, sequence_metrics
from oligoark.reconstruct import ReconstructionResult


class StaticReconstructor:
    def __init__(self, consensus: str) -> None:
        self.consensus = consensus

    def reconstruct(self, reads: list[str]) -> ReconstructionResult:
        return ReconstructionResult([self.consensus], [len(reads)])


def test_binary_dna_roundtrip() -> None:
    data = bytes(range(256))
    assert dna_to_bytes(bytes_to_dna(data)) == data


def test_archive_roundtrip() -> None:
    raw = os.urandom(1000)
    archive = archive_bytes(raw, ArchiveConfig(chunk_size=64, rs_nsym=8, parity_group_size=4))
    assert recover_bytes(archive) == raw


def test_single_dropout_per_parity_group_recovers() -> None:
    raw = bytes(range(200))
    cfg = ArchiveConfig(chunk_size=32, rs_nsym=8, parity_group_size=3)
    archive = archive_bytes(raw, cfg)
    reads = [strand for i, strand in enumerate(archive.strands) if i != 1]
    assert recover_bytes(archive, reads) == raw


def test_unrecoverable_dropout_raises() -> None:
    raw = bytes(range(200))
    cfg = ArchiveConfig(chunk_size=32, rs_nsym=8, parity_group_size=3)
    archive = archive_bytes(raw, cfg)
    reads = [strand for i, strand in enumerate(archive.strands) if i not in {0, 1}]
    with pytest.raises(ValueError, match="not recoverable"):
        recover_bytes(archive, reads)


def test_sequence_metrics() -> None:
    metrics = sequence_metrics("AACCGGTTTT")
    assert metrics.gc_fraction == pytest.approx(0.4)
    assert metrics.max_homopolymer == 4


def test_graph_reconstruction_is_used_when_direct_decode_fails() -> None:
    raw = b"graph reconstruction baseline"
    archive = archive_bytes(raw, ArchiveConfig(chunk_size=64, rs_nsym=0, parity_group_size=8))
    original = archive.strands[0]
    corrupted_reads = []
    for position in (12, 20, 28):
        chars = list(original)
        chars[position] = "A" if chars[position] != "A" else "C"
        corrupted_reads.append("".join(chars))

    recovered, report = recover_from_reads(archive, corrupted_reads, similarity_threshold=0.95)
    assert recovered == raw
    assert report.graph_reconstruction_used is True
    assert report.reconstruction_strategy == "GraphConsensusReconstructor"
    assert report.verified_sha256 is True


def test_custom_reconstructor_can_be_injected() -> None:
    raw = b"pluggable reconstruction"
    archive = archive_bytes(raw, ArchiveConfig(chunk_size=64, rs_nsym=0, parity_group_size=8))
    original = archive.strands[0]
    broken = original[:-4]
    recovered, report = recover_from_reads(
        archive,
        [broken],
        reconstructor=StaticReconstructor(original),
    )
    assert recovered == raw
    assert report.graph_reconstruction_used is False
    assert report.reconstruction_strategy == "StaticReconstructor"


def test_archive_validation_rejects_invalid_metadata_types() -> None:
    archive = archive_bytes(b"validation")
    archive.metadata["data_strands"] = "1"
    with pytest.raises(ValueError, match="data_strands metadata must be an integer"):
        archive.validate()


def test_archive_config_rejects_string_boolean() -> None:
    with pytest.raises(ValueError, match="adaptive_masks must be a boolean"):
        ArchiveConfig.from_mapping({"adaptive_masks": "false"})

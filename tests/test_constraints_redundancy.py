import pytest

from oligoark.archive import ArchiveConfig, archive_bytes, archive_statistics, recover_bytes
from oligoark.dna import SequenceConstraintError, SequenceConstraints
from oligoark.framing import decode_frame


def test_created_archive_satisfies_hard_constraints() -> None:
    config = ArchiveConfig(
        chunk_size=32,
        rs_nsym=8,
        parity_group_size=3,
        min_gc_fraction=0.35,
        max_gc_fraction=0.65,
        max_homopolymer=4,
        mask_search_limit=96,
    )
    archive = archive_bytes(bytes(range(128)), config)
    constraints = SequenceConstraints(0.35, 0.65, 4)
    assert all(constraints.accepts(strand) for strand in archive.strands)
    assert archive_statistics(archive).constraint_pass_rate == 1.0


def test_impossible_single_mask_constraints_raise() -> None:
    config = ArchiveConfig(
        chunk_size=32,
        rs_nsym=8,
        parity_group_size=3,
        adaptive_masks=False,
        min_gc_fraction=1.0,
        max_gc_fraction=1.0,
        max_homopolymer=1,
        mask_search_limit=1,
    )
    with pytest.raises(SequenceConstraintError):
        archive_bytes(b"constraint failure should be explicit", config)


def test_fountain_redundancy_is_in_real_archive_pipeline() -> None:
    payload = bytes(range(128))
    config = ArchiveConfig(
        chunk_size=32,
        rs_nsym=8,
        parity_group_size=3,
        redundancy_scheme="fountain",
        fountain_redundancy=3.0,
        fountain_seed=1,
        mask_search_limit=96,
    )
    archive = archive_bytes(payload, config)
    assert int(archive.metadata["fountain_strands"]) > 0
    assert int(archive.metadata["parity_strands"]) == 0

    # Remove one data strand while retaining all fountain symbols.
    reads = archive.strands[1:]
    assert recover_bytes(archive, reads) == payload
    fountain_frames = [
        decode_frame(strand, rs_nsym=config.rs_nsym)
        for strand in archive.strands[int(archive.metadata["data_strands"]) :]
    ]
    assert all(frame.is_fountain for frame in fountain_frames)


def test_hybrid_redundancy_reports_both_types() -> None:
    config = ArchiveConfig(
        chunk_size=32,
        rs_nsym=8,
        parity_group_size=3,
        redundancy_scheme="hybrid",
        fountain_redundancy=0.75,
    )
    archive = archive_bytes(b"hybrid redundancy" * 8, config)
    stats = archive_statistics(archive)
    assert stats.parity_strands > 0
    assert stats.fountain_strands > 0
    assert stats.redundancy_scheme == "hybrid"

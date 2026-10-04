from oligoark.archive import ArchiveConfig, archive_bytes, recover_bytes
from oligoark.reconstruct import (
    TraceConsensusReconstructor,
    iterative_trace_consensus,
)
from oligoark.validation import compare_reconstruction_modes


def _insert(sequence: str, position: int, base: str) -> str:
    return sequence[:position] + base + sequence[position:]


def _delete(sequence: str, position: int) -> str:
    return sequence[:position] + sequence[position + 1 :]


def _archive(payload: bytes):
    return archive_bytes(
        payload,
        ArchiveConfig(
            chunk_size=96,
            rs_nsym=0,
            parity_group_size=8,
            redundancy_scheme="none",
            mask_search_limit=128,
        ),
    )


def test_iterative_trace_consensus_recovers_reference_from_indels() -> None:
    archive = _archive(b"iterative trace consensus")
    reference = archive.strands[0]
    cluster = [
        _insert(reference, 72, "A"),
        _insert(reference, 88, "C"),
        _insert(reference, 104, "G"),
        _delete(reference, 80),
        _delete(reference, 96),
        _delete(reference, 112),
        _delete(reference, 128),
    ]
    consensus = iterative_trace_consensus(cluster, rounds=4)
    assert recover_bytes(archive, [consensus]) == b"iterative trace consensus"


def test_trace_reconstructor_emits_verified_rescue_candidates() -> None:
    archive = _archive(b"trace graph rescue")
    reference = archive.strands[0]
    reads = [
        _insert(reference, 72, "A"),
        _insert(reference, 92, "C"),
        _insert(reference, 112, "G"),
        _delete(reference, 82),
        _delete(reference, 102),
        _delete(reference, 122),
        _delete(reference, 142),
    ]
    reconstructor = TraceConsensusReconstructor(thresholds=(0.96, 0.92, 0.88), rounds=4)
    result = reconstructor.reconstruct(reads)
    assert result.node_count == len(reads)
    assert result.candidate_pairs == len(reads) * (len(reads) - 1) // 2
    assert result.consensus_reads
    assert recover_bytes(archive, reads + result.consensus_reads) == b"trace graph rescue"

    comparison = compare_reconstruction_modes(archive, reads, threshold=0.92)
    assert comparison.direct_recovered is False
    assert comparison.trace.recovered is True
    assert comparison.trace_rescued_direct_failure is True

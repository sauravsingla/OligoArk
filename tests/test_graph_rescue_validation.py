from oligoark.archive import ArchiveConfig, archive_bytes
from oligoark.validation import compare_reconstruction_modes


def _substitute(sequence: str, position: int) -> str:
    current = sequence[position]
    replacement = next(base for base in "ACGT" if base != current)
    return sequence[:position] + replacement + sequence[position + 1 :]


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


def test_substitution_reads_prove_direct_failure_then_graph_rescue() -> None:
    archive = _archive(b"substitution rescue")
    strand = archive.strands[0]
    reads = [
        _substitute(strand, position)
        for position in (72, 84, 96, 108, 120, 132, 144)
    ]
    comparison = compare_reconstruction_modes(archive, reads, threshold=0.96)
    assert comparison.direct_recovered is False
    assert comparison.alignment.recovered is True
    assert comparison.alignment.verified_sha256 is True
    assert comparison.alignment_rescued_direct_failure is True
    assert comparison.alignment.node_count == len(reads)
    assert comparison.alignment.candidate_pairs == len(reads) * (len(reads) - 1) // 2
    assert comparison.alignment.edge_count > 0
    assert comparison.alignment.component_count >= 1
    assert comparison.alignment.consensus_count >= 1

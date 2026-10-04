"""Compare direct, medoid-graph, and alignment-graph recovery."""

from oligoark import ArchiveConfig, archive_bytes, compare_reconstruction_modes


def substitute(sequence: str, position: int) -> str:
    current = sequence[position]
    replacement = next(base for base in "ACGT" if base != current)
    return sequence[:position] + replacement + sequence[position + 1 :]


payload = b"graph reconstruction diagnostic example"
archive = archive_bytes(
    payload,
    ArchiveConfig(
        chunk_size=96,
        rs_nsym=0,
        parity_group_size=8,
        redundancy_scheme="none",
        mask_search_limit=128,
    ),
)
strand = archive.strands[0]
positions = (64, 76, 88, 100, 112)
reads = [substitute(strand, position) for position in positions if position < len(strand)]

comparison = compare_reconstruction_modes(archive, reads, threshold=0.96)
print(comparison.to_dict())

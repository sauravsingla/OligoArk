"""Encode and recover an archive using integrated XOR + fountain redundancy."""

from oligoark import ArchiveConfig, archive_bytes, archive_statistics, recover_bytes

payload = b"OligoArk hybrid redundancy example" * 4
config = ArchiveConfig(
    chunk_size=48,
    rs_nsym=12,
    parity_group_size=3,
    redundancy_scheme="hybrid",
    fountain_redundancy=0.5,
    min_gc_fraction=0.35,
    max_gc_fraction=0.65,
    max_homopolymer=4,
    mask_search_limit=96,
)
archive = archive_bytes(payload, config)
assert recover_bytes(archive) == payload
stats = archive_statistics(archive)
assert stats.fountain_strands > 0
assert stats.parity_strands > 0
assert stats.constraint_pass_rate == 1.0
print(stats.to_dict())

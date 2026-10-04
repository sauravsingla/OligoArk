"""Demonstrate graph-consensus rescue on controlled substitution errors."""

from oligoark.archive import ArchiveConfig, archive_bytes, recover_from_reads

payload = b"graph reconstruction demonstration"
archive = archive_bytes(payload, ArchiveConfig(chunk_size=48, rs_nsym=0, parity_group_size=8))
original = archive.strands[0]
reads: list[str] = []
for position in (20, 40, 60):
    damaged = list(original)
    damaged[position] = "A" if damaged[position] != "A" else "C"
    reads.append("".join(damaged))

recovered, report = recover_from_reads(archive, reads, similarity_threshold=0.95)
assert recovered == payload
assert report.graph_reconstruction_used
print(report.to_dict())

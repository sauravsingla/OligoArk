from oligoark.archive import archive_bytes, recover_bytes
from oligoark.simulator import SimulationConfig, simulate_channel

payload = b"OligoArk reproducible end-to-end example"
archive = archive_bytes(payload)
reads = simulate_channel(archive.strands, SimulationConfig(seed=42))
assert recover_bytes(archive, reads) == payload
print(f"verified sha256={archive.metadata['sha256']}")

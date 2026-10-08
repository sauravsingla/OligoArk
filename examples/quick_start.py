"""Small deterministic OligoArk encode -> simulated channel -> recover example.

Run from the repository root after installing with: python -m pip install -e .
"""
import hashlib

from oligoark.archive import archive_bytes, recover_bytes
from oligoark.simulator import SimulationConfig, simulate_channel


def main() -> None:
    payload = b"Hello from OligoArk! This is a reproducible DNA storage demo."
    expected_hash = hashlib.sha256(payload).hexdigest()

    archive = archive_bytes(payload)
    reads = simulate_channel(archive.strands, SimulationConfig(seed=42))
    recovered = recover_bytes(archive, reads)
    actual_hash = hashlib.sha256(recovered).hexdigest()

    assert recovered == payload, "Recovered bytes differ from the original"
    assert actual_hash == expected_hash, "SHA-256 mismatch"

    print("PASS: original data recovered exactly")
    print(f"Bytes: {len(payload)}")
    print(f"SHA-256: {actual_hash}")
    print("Note: simulated DNA channel; no physical synthesis or sequencing.")


if __name__ == "__main__":
    main()

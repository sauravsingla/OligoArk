from oligoark.archive import ArchiveConfig, archive_bytes, recover_bytes
from oligoark.framing import decode_frame
from oligoark.simulator import SimulationConfig, simulate_channel


def test_reed_solomon_corrects_small_substitution_damage() -> None:
    archive = archive_bytes(b"hello-reed-solomon",
        ArchiveConfig(chunk_size=32, rs_nsym=12, parity_group_size=4))
    dna = list(archive.strands[0])
    dna[-5] = "A" if dna[-5] != "A" else "C"
    frame = decode_frame("".join(dna), rs_nsym=12)
    assert frame.payload == b"hello-reed-solomon"


def test_simulator_is_seeded() -> None:
    strands = ["ACGT" * 20]
    cfg = SimulationConfig(substitution_rate=0.05, insertion_rate=0.02, seed=42)
    assert simulate_channel(strands, cfg) == simulate_channel(strands, cfg)


def test_zero_noise_simulation_recovers() -> None:
    raw = b"software-only DNA channel simulation" * 5
    archive = archive_bytes(raw)
    reads = simulate_channel(archive.strands, SimulationConfig(seed=9))
    assert recover_bytes(archive, reads) == raw

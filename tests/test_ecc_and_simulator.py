import pytest

from oligoark.archive import ArchiveConfig, archive_bytes, recover_bytes
from oligoark.framing import decode_frame
from oligoark.simulator import SimulationConfig, simulate_channel


def test_reed_solomon_corrects_small_substitution_damage() -> None:
    config = ArchiveConfig(chunk_size=32, rs_nsym=12, parity_group_size=4)
    archive = archive_bytes(b"hello-reed-solomon", config)
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


def test_simulator_supports_deterministic_multi_trace_coverage() -> None:
    config = SimulationConfig(seed=17, copies_per_strand=4)
    first = simulate_channel(["ACGTACGT"], config)
    second = simulate_channel(["ACGTACGT"], config)
    assert first == second
    assert len(first) == 4
    assert first == ["ACGTACGT"] * 4

    with pytest.raises(ValueError, match="copies_per_strand"):
        simulate_channel(["ACGT"], SimulationConfig(copies_per_strand=0))

from pathlib import Path

import pytest

from oligoark import cli
from oligoark.logging_utils import configure_logging


def _run_cli(argv: list[str], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("sys.argv", ["oligoark", *argv])
    cli.main()


def test_cli_end_to_end(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = tmp_path / "payload.bin"
    archive = tmp_path / "payload.oligoark.json"
    reads = tmp_path / "reads.txt"
    recovered = tmp_path / "recovered.bin"
    recovered_reads = tmp_path / "recovered-reads.bin"
    source.write_bytes(b"cli end to end")

    _run_cli(
        [
            "archive",
            str(source),
            "--output",
            str(archive),
            "--redundancy-scheme",
            "hybrid",
            "--fountain-redundancy",
            "0.5",
        ],
        monkeypatch,
    )
    _run_cli(["inspect", str(archive)], monkeypatch)
    _run_cli(["recover", str(archive), "--output", str(recovered)], monkeypatch)
    _run_cli(["simulate", str(archive), "--output", str(reads), "--seed", "12"], monkeypatch)
    _run_cli(
        [
            "recover-reads",
            str(archive),
            str(reads),
            "--output",
            str(recovered_reads),
        ],
        monkeypatch,
    )
    assert recovered.read_bytes() == source.read_bytes()
    assert recovered_reads.read_bytes() == source.read_bytes()


def test_cli_recommend_policy_plan_and_optimizer(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    economics = tmp_path / "economics.json"
    economics.write_text(
        '{"storage_cost_index":{"ssd":0.8,"object_archive":0.2,"tape":0.3,'
        '"dna_future":0.5},"retrieval_cost_index":{"ssd":0.1,'
        '"object_archive":0.4,"tape":0.6,"dna_future":0.8}}',
        encoding="utf-8",
    )
    payload = tmp_path / "optimizer.bin"
    payload.write_bytes(b"cli optimizer payload")

    _run_cli(
        [
            "recommend",
            "--retention-years",
            "100",
            "--economics-json",
            str(economics),
        ],
        monkeypatch,
    )
    _run_cli(["policy", "--substitution", "0.01"], monkeypatch)
    _run_cli(
        [
            "plan",
            "--retention-years",
            "100",
            "--substitution",
            "0.01",
            "--dropout",
            "0.05",
            "--economics-json",
            str(economics),
        ],
        monkeypatch,
    )
    _run_cli(
        [
            "optimize-plan",
            str(payload),
            "--retention-years",
            "100",
            "--substitution",
            "0.002",
            "--max-candidates",
            "4",
            "--seeds",
            "2026",
        ],
        monkeypatch,
    )


def test_logging_configuration_validation() -> None:
    configure_logging("INFO")
    with pytest.raises(ValueError, match="Unknown logging level"):
        configure_logging("not-a-level")

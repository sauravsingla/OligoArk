"""Opt-in codec worker diagnostics retain evidence without altering trial verdicts."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

import run_matched_codec_scale as benchmark


def test_progress_log_records_worker_rss_and_exact_trial_outcomes(
    monkeypatch, tmp_path: Path
) -> None:
    progress = tmp_path / "worker-progress.jsonl"
    monkeypatch.setenv("OLIGOARK_CODEC_PROGRESS_FILE", str(progress))
    row = benchmark._condition_row(
        method="oligoark-compact-hybrid",
        payload=b"x",
        sequences=("ACGT" * 38,),
        decode=lambda _reads: b"x",
        condition_name="clean",
        trials=2,
        redundancy=0.25,
        target_max_strand_nt=152,
        data_units=1,
        encoded_nucleotides=152,
        logical_bits_per_nucleotide=0.05,
        encode_seconds=0.01,
        condition_timeout_seconds=0,
    )
    assert row["successes"] == 2
    assert all(t["sha256_verified"] is True for t in row["trial_results"])
    events = [json.loads(line) for line in progress.read_text().splitlines()]
    assert [event["event"] for event in events] == [
        "condition-start",
        "trial-start",
        "trial-exit",
        "trial-start",
        "trial-exit",
        "condition-complete",
    ]
    assert all(event["peak_worker_rss_mib"] >= 0 for event in events)
    assert all(event["pid"] > 0 and event["timestamp_unix"] > 0 for event in events)
    recorded = [event["outcome_recorded"] for event in events if event["event"] == "trial-exit"]
    assert recorded == [True, True]


def test_worker_lifecycle_preserves_machine_readable_stdout(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    progress = tmp_path / "worker-progress.jsonl"
    monkeypatch.setenv("OLIGOARK_CODEC_PROGRESS_FILE", str(progress))
    monkeypatch.setattr(
        benchmark, "_worker", lambda *args: [{"condition": "clean", "successes": 1}]
    )
    monkeypatch.setattr(
        sys, "argv",
        ["run_matched_codec_scale.py", "--worker", "--profile", "ci",
         "--method", "dna-fountain-cleanroom", "--size", "1024",
         "--conditions", "clean", "--trials", "1"],
    )
    benchmark.main()
    assert json.loads(capsys.readouterr().out) == [
        {"condition": "clean", "successes": 1}
    ]
    events = [json.loads(line)["event"] for line in progress.read_text().splitlines()]
    assert events == ["worker-start", "worker-complete"]


def test_optional_diagnostic_io_failure_does_not_reclassify_trials(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    monkeypatch.setenv(
        "OLIGOARK_CODEC_PROGRESS_FILE",
        str(tmp_path / "missing-directory" / "worker-progress.jsonl"),
    )
    benchmark._progress_event("worker-start")
    assert "codec progress logging unavailable" in capsys.readouterr().err


def test_disabled_progress_does_not_create_files(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.delenv("OLIGOARK_CODEC_PROGRESS_FILE", raising=False)
    benchmark._progress_event("worker-start")
    assert list(tmp_path.iterdir()) == []

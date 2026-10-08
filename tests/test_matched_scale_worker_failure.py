"""Regression: worker encoding failures retain every requested negative trial."""

import importlib.util
import subprocess
from pathlib import Path
from types import SimpleNamespace


def test_matched_scale_worker_failure_retains_negative_trials(monkeypatch):
    benchmark_dir = Path(__file__).resolve().parents[1] / "benchmarks"
    monkeypatch.syspath_prepend(str(benchmark_dir))
    spec = importlib.util.spec_from_file_location(
        "run_matched_codec_scale", benchmark_dir / "run_matched_codec_scale.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    def failed_worker(*args, **kwargs):
        return SimpleNamespace(returncode=1, stderr="SequenceConstraintError", stdout="")

    monkeypatch.setattr(subprocess, "run", failed_worker)
    rows = module._isolated(
        "oligoark-compact-v3-hybrid", 100 * 1024 * 1024, 0.25,
        10, ("clean", "indel-low"), 1800, 300,
    )
    assert len(rows) == 2
    for row in rows:
        assert row["encoding_failed"] is True
        assert row["successes"] == 0
        assert row["failed_trials"] == 10
        assert len(row["trial_results"]) == 10
        assert [t["trial"] for t in row["trial_results"]] == list(range(10))
        assert all(t["sha256_verified"] is False for t in row["trial_results"])
        assert all(t["outcome_unavailable"] is True for t in row["trial_results"])
        assert "SequenceConstraintError" in row["error"]


def test_matched_scale_worker_timeout_retains_unverified_trials(monkeypatch):
    benchmark_dir = Path(__file__).resolve().parents[1] / "benchmarks"
    monkeypatch.syspath_prepend(str(benchmark_dir))
    spec = importlib.util.spec_from_file_location(
        "run_matched_codec_scale_timeout", benchmark_dir / "run_matched_codec_scale.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    def timed_out_worker(*args, **kwargs):
        raise subprocess.TimeoutExpired(cmd="benchmark", timeout=10)

    monkeypatch.setattr(subprocess, "run", timed_out_worker)
    rows = module._isolated(
        "dna-fountain-cleanroom", 100 * 1024 * 1024, 0.25,
        10, ("clean", "mixed"), 1800, 300,
    )
    assert len(rows) == 2
    for row in rows:
        assert row["timed_out"] is True
        assert row["successes"] == 0
        assert row["failed_trials"] == 10
        assert len(row["trial_results"]) == 10
        assert [trial["seed"] for trial in row["trial_results"]] == [
            20_260_000 + index for index in range(10)
        ]
        assert all(not trial["sha256_verified"] for trial in row["trial_results"])
        assert all(trial["outcome_unavailable"] for trial in row["trial_results"])

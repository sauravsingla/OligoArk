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
        assert all(t["not_executed"] is True for t in row["trial_results"])
        assert "SequenceConstraintError" in row["error"]

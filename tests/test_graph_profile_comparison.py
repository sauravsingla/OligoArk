"""Validate comparison refuses incomparable or semantically changed benchmarks."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks" / "compare_graph_profiles.py"
_spec = importlib.util.spec_from_file_location("compare_graph_profiles", _SCRIPT)
assert _spec is not None and _spec.loader is not None
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
compare = _module.compare


def _profile(seconds: float = 2.0, memory: int = 100) -> dict[str, object]:
    return {
        "parameters": {"seed": 17, "strands": 12, "copies": 4, "length": 80, "repeats": 3},
        "results": {
            "graph": {
                "median_seconds": seconds,
                "timings_seconds": [seconds] * 3,
                "peak_tracemalloc_bytes": memory,
                "result_sha256": "a" * 64,
            }
        },
    }


def test_matched_comparison() -> None:
    result = compare(_profile(), _profile(1.0, 50))
    assert result["results"]["graph"]["time_ratio_baseline_over_candidate"] == 2.0
    assert result["results"]["graph"]["memory_ratio_baseline_over_candidate"] == 2.0
    assert "not establish statistical significance" in result["limitations"]


def test_rejects_parameter_mismatch() -> None:
    after = _profile()
    after["parameters"]["seed"] = 18
    with pytest.raises(ValueError, match="parameters differ"):
        compare(_profile(), after)


def test_rejects_fingerprint_mismatch() -> None:
    after = _profile()
    after["results"]["graph"]["result_sha256"] = "b" * 64
    with pytest.raises(ValueError, match="fingerprint changed"):
        compare(_profile(), after)


def test_rejects_invalid_timing() -> None:
    after = _profile()
    after["results"]["graph"]["median_seconds"] = 0.0
    with pytest.raises(ValueError, match="invalid timings"):
        compare(_profile(), after)


def test_rejects_missing_case() -> None:
    after = _profile()
    after["results"].clear()
    with pytest.raises(ValueError, match="benchmark results are missing"):
        compare(_profile(), after)

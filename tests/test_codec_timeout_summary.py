"""Regression checks for timeout summary, including failed and incomplete trials."""
from __future__ import annotations

import importlib.util
from pathlib import Path


def _module():
    path = Path(__file__).resolve().parents[1] / "benchmarks" / "summarize_codec_timeouts.py"
    spec = importlib.util.spec_from_file_location("codec_timeout_summary", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_preserves_timeouts_failures_and_hotspot_counts() -> None:
    rows = [{"method": "compact", "condition": "indel-low", "trial_results": [
        {"sha256_verified": False, "timed_out": True, "timeout_stage": "decode",
         "timeout_stack": [{"file": "x.py", "line": 5, "function": "slow"}]},
        {"sha256_verified": False, "timed_out": True, "timeout_stage": "decode",
         "timeout_stack": [{"file": "x.py", "line": 5, "function": "slow"}]},
        {"sha256_verified": False, "error": "checksum mismatch"},
        {"sha256_verified": True},
    ]}]
    result = _module().summarize({"rows": rows})["groups"]["compact / indel-low"]
    assert result["trials"] == 4
    assert result["successes"] == 1
    assert result["timeouts"] == 2
    assert result["timeout_stages"] == {"decode": 2}
    assert result["interrupted_locations"] == {"x.py:5:slow": 2}


def test_empty_stack_and_empty_results() -> None:
    module = _module()
    assert module.summarize([])["groups"] == {}
    group = module.summarize([{"trial_results": [
        {"timed_out": True, "timeout_stack": []},
    ]}])["groups"]["unknown / unknown"]
    assert group["timeouts"] == 1
    assert group["interrupted_locations"] == {}

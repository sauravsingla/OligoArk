from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

from compare_rs_lookup_e2e import summarize, validate_rows


def _rows(*, timeout: bool = False) -> list[dict[str, object]]:
    payload_sha = hashlib.sha256(b"payload").hexdigest()
    return [
        {
            "method": "oligoark-compact-v3-hybrid",
            "size_bytes": 10 * 1024 * 1024,
            "condition": condition,
            "trials": 2,
            "target_max_strand_nt": 152,
            "redundancy_budget": 0.25,
            "successes": 1 if timeout else 2,
            "failed_trials": 1 if timeout else 0,
            "completed_trials": 1 if timeout else 2,
            "timeout_trials": 1 if timeout else 0,
            "payload_sha256": payload_sha,
            "strand_count": 42,
            "encoded_nucleotides": 6384,
            "rs_nsym": 12,
            "encode_seconds": 1.5,
            "peak_rss_mib": 100.0,
            "trial_results": [
                {
                    "trial": 0, "seed": 20260000, "success": True,
                    "sha256_verified": True, "decode_seconds": 0.5,
                },
                (
                    {
                        "trial": 1, "seed": 20260001, "success": False,
                        "sha256_verified": False, "timed_out": True, "timeout_seconds": 75,
                    }
                    if timeout else
                    {
                        "trial": 1, "seed": 20260001, "success": True,
                        "sha256_verified": True, "decode_seconds": 0.6,
                    }
                ),
            ],
        }
        for condition in ("indel-low", "mixed")
    ]


def test_validator_preserves_real_timeouts() -> None:
    report = validate_rows(_rows(timeout=True), method="oligoark-compact-v3-hybrid", trials=2)
    assert report["indel-low"]["timeouts"] == 1
    assert report["indel-low"]["failures"] == 1


def test_validator_rejects_unavailable_or_missing_trials() -> None:
    rows = _rows()
    rows[0]["trial_results"][0]["outcome_unavailable"] = True  # type: ignore[index]
    with pytest.raises(ValueError, match="unavailable"):
        validate_rows(rows, method="oligoark-compact-v3-hybrid", trials=2)
    rows = _rows()
    rows[0]["trial_results"].pop()  # type: ignore[union-attr]
    with pytest.raises(ValueError, match="missing trial records"):
        validate_rows(rows, method="oligoark-compact-v3-hybrid", trials=2)


def test_paired_repeats_and_complete_denominators(tmp_path: Path) -> None:
    for repeat in range(2):
        root = tmp_path / f"repeat-{repeat}"
        root.mkdir()
        for version in ("baseline", "optimized"):
            (root / f"{version}.json").write_text(json.dumps(_rows(timeout=version == "baseline")))
            (root / f"{version}.exit-code").write_text("0\n")
    report = summarize(
        tmp_path, method="oligoark-compact-v3-hybrid", trials=2, repeats=2,
        baseline_sha="a" * 40, optimized_sha="b" * 40,
    )
    assert report["order_by_repeat"] == [
        ["baseline", "optimized"], ["optimized", "baseline"],
    ]
    assert report["results"]["baseline"]["mixed"]["timeout_trials"] == 2
    assert report["results"]["optimized"]["mixed"]["exact_recoveries"] == 4
    (tmp_path / "repeat-1" / "optimized.exit-code").write_text("124\n")
    with pytest.raises(ValueError, match="worker did not complete"):
        summarize(
            tmp_path, method="oligoark-compact-v3-hybrid", trials=2, repeats=2,
            baseline_sha="a" * 40, optimized_sha="b" * 40,
        )

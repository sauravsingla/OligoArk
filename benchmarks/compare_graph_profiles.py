"""Compare matched graph/consensus profiles without claiming unverified speedups.

Only compares runs with identical parameters and output fingerprints. Timing
ratios are descriptive, not evidence of significance or hardware equivalence.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def compare(baseline: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    """Return per-case timing/memory ratios after strict semantic validation."""
    if baseline.get("parameters") != candidate.get("parameters"):
        raise ValueError("benchmark parameters differ; cannot compare")
    params = baseline.get("parameters")
    if not isinstance(params, dict) or not params:
        raise ValueError("benchmark parameters are missing")
    before = baseline.get("results")
    after = candidate.get("results")
    if (
        not isinstance(before, dict)
        or not before
        or not isinstance(after, dict)
        or not after
    ):
        raise ValueError("benchmark results are missing")
    if before.keys() != after.keys():
        raise ValueError("benchmark cases differ; cannot compare")
    results: dict[str, dict[str, float | int | str]] = {}
    for case in sorted(before):
        old, new = before[case], after[case]
        if not isinstance(old, dict) or not isinstance(new, dict):
            raise ValueError(f"invalid case: {case}")
        old_hash = old.get("result_sha256")
        if not isinstance(old_hash, str) or len(old_hash) != 64:
            raise ValueError(f"invalid baseline fingerprint for {case}")
        if old_hash != new.get("result_sha256"):
            raise ValueError(f"result fingerprint changed for {case}; verify semantics first")
        times = []
        for row in (old, new):
            value = row.get("median_seconds")
            samples = row.get("timings_seconds")
            if (
                not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value <= 0
                or not isinstance(samples, list)
                or not samples
                or any(
                    not isinstance(t, (int, float)) or not math.isfinite(t) or t <= 0
                    for t in samples
                )
            ):
                raise ValueError(f"invalid timings for {case}")
            times.append(float(value))
        old_memory = old.get("peak_tracemalloc_bytes")
        new_memory = new.get("peak_tracemalloc_bytes")
        if (
            not isinstance(old_memory, int)
            or not isinstance(new_memory, int)
            or old_memory <= 0
            or new_memory <= 0
        ):
            raise ValueError(f"invalid traced memory for {case}")
        results[case] = {
            "baseline_median_seconds": times[0],
            "candidate_median_seconds": times[1],
            "time_ratio_baseline_over_candidate": times[0] / times[1],
            "baseline_peak_tracemalloc_bytes": old_memory,
            "candidate_peak_tracemalloc_bytes": new_memory,
            "memory_ratio_baseline_over_candidate": old_memory / new_memory,
            "result_sha256": old_hash,
        }
    return {
        "parameters": params,
        "results": results,
        "limitations": (
            "Descriptive matched-parameter ratios only. Check identical benchmark code, "
            "Python version, optional acceleration, CPU allocation and runner environment "
            "before interpreting timing; traced Python allocations are not process RSS. "
            "These ratios do not establish statistical significance or exact archive recovery."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path, default=Path("performance-results/comparison.json"))
    args = parser.parse_args()
    baseline = json.loads(args.baseline.read_text(encoding="utf-8"))
    candidate = json.loads(args.candidate.read_text(encoding="utf-8"))
    report = compare(baseline, candidate)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

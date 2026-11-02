"""Validate paired unprofiled codec outcomes across pinned RS-lookup revisions.

This report is research evidence, not a speedup or recovery-improvement gate.
Missing/partial workers are *not* counted as verified recovery failures.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics
from pathlib import Path
from typing import Any

CONDITIONS = ("indel-low", "mixed")
VERSIONS = ("baseline", "optimized")
SIZE = 10 * 1024 * 1024
DEADLINE = 75


def validate_rows(rows: Any, *, method: str, trials: int) -> dict[str, Any]:
    """Reject incomplete evidence while retaining genuine failed recoveries."""
    if not isinstance(rows, list) or len(rows) != len(CONDITIONS):
        raise ValueError("expected exactly two condition rows")
    if not all(isinstance(row, dict) for row in rows):
        raise ValueError("invalid condition row")
    if {row.get("condition") for row in rows} != set(CONDITIONS):
        raise ValueError("missing or duplicated indel-low/mixed condition")
    validated: dict[str, Any] = {}
    for row in rows:
        condition = row["condition"]
        if row.get("method") != method or row.get("size_bytes") != SIZE:
            raise ValueError(f"{condition}: method or size mismatch")
        if row.get("trials") != trials or row.get("timeout_seconds", DEADLINE) != DEADLINE:
            raise ValueError(f"{condition}: trials or per-trial deadline mismatch")
        if row.get("target_max_strand_nt") != 152 or row.get("redundancy_budget") != 0.25:
            raise ValueError(f"{condition}: strand ceiling or redundancy mismatch")
        trial_rows = row.get("trial_results")
        if not isinstance(trial_rows, list) or len(trial_rows) != trials:
            raise ValueError(f"{condition}: missing trial records")
        if sorted(item.get("trial") for item in trial_rows) != list(range(trials)):
            raise ValueError(f"{condition}: missing or duplicate trial indexes")
        successes = 0
        timeouts = 0
        successful_seconds: list[float] = []
        for item in trial_rows:
            if item.get("seed") != 20_260_000 + item["trial"]:
                raise ValueError(f"{condition}: seed mismatch")
            if item.get("outcome_unavailable") or item.get("not_executed"):
                raise ValueError(f"{condition}: worker outcome unavailable")
            success = item.get("success")
            verified = item.get("sha256_verified")
            if (
                not isinstance(success, bool)
                or not isinstance(verified, bool)
                or success != verified
            ):
                raise ValueError(f"{condition}: missing/inconsistent exact recovery verdict")
            timed_out = item.get("timed_out", False)
            if not isinstance(timed_out, bool):
                raise ValueError(f"{condition}: invalid timeout flag")
            if timed_out:
                timeouts += 1
                if success or item.get("timeout_seconds") != DEADLINE:
                    raise ValueError(f"{condition}: timeout incorrectly recorded")
            if success:
                successes += 1
                seconds = item.get("decode_seconds")
                if not isinstance(seconds, (int, float)) or seconds < 0:
                    raise ValueError(f"{condition}: missing successful decode duration")
                successful_seconds.append(float(seconds))
        if (
            row.get("successes") != successes
            or row.get("failed_trials") != trials - successes
            or row.get("timeout_trials") != timeouts
            or row.get("completed_trials") != trials - timeouts
        ):
            raise ValueError(f"{condition}: aggregate accounting mismatch")
        payload_hash = row.get("payload_sha256")
        if not isinstance(payload_hash, str) or len(payload_hash) != 64:
            raise ValueError(f"{condition}: missing payload SHA-256")
        try:
            int(payload_hash, 16)
        except ValueError as exc:
            raise ValueError(f"{condition}: invalid payload SHA-256") from exc
        validated[condition] = {
            "successes": successes,
            "failures": trials - successes,
            "timeouts": timeouts,
            "successful_decode_seconds": successful_seconds,
            "encode_seconds": row.get("encode_seconds"),
            "peak_rss_mib": row.get("peak_rss_mib"),
            "payload_sha256": payload_hash,
            "strand_count": row.get("strand_count"),
            "encoded_nucleotides": row.get("encoded_nucleotides"),
            "rs_nsym": row.get("rs_nsym"),
        }
    if validated[CONDITIONS[0]]["payload_sha256"] != validated[CONDITIONS[1]]["payload_sha256"]:
        raise ValueError("payload differs between conditions")
    return validated


def summarize(
    root: Path, *, method: str, trials: int, repeats: int,
    baseline_sha: str, optimized_sha: str,
) -> dict[str, Any]:
    if repeats < 2 or trials < 1:
        raise ValueError("at least two alternating repeats and one trial required")
    samples: dict[str, dict[str, list[dict[str, Any]]]] = {
        version: {condition: [] for condition in CONDITIONS} for version in VERSIONS
    }
    payload_hashes: set[str] = set()
    strand_fingerprints: set[tuple[Any, Any, Any]] = set()
    order: list[list[str]] = []
    for repeat in range(repeats):
        versions = list(VERSIONS if repeat % 2 == 0 else reversed(VERSIONS))
        order.append(versions)
        for version in versions:
            path = root / f"repeat-{repeat}" / f"{version}.json"
            status_path = root / f"repeat-{repeat}" / f"{version}.exit-code"
            if not path.is_file() or not status_path.is_file():
                raise ValueError(f"missing worker evidence: {path}")
            if status_path.read_text().strip() != "0":
                raise ValueError(f"worker did not complete: {status_path}")
            try:
                rows = json.loads(path.read_text())
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(f"incomplete worker JSON: {path}") from exc
            validated = validate_rows(rows, method=method, trials=trials)
            for condition in CONDITIONS:
                sample = validated[condition]
                samples[version][condition].append(sample)
                payload_hashes.add(sample["payload_sha256"])
                strand_fingerprints.add((
                    sample["strand_count"], sample["encoded_nucleotides"], sample["rs_nsym"]
                ))
    if len(payload_hashes) != 1 or len(strand_fingerprints) != 1:
        raise ValueError("paired revisions differ in payload or encoded strand layout")
    summary: dict[str, Any] = {}
    for version in VERSIONS:
        summary[version] = {}
        for condition in CONDITIONS:
            records = samples[version][condition]
            completed_seconds = [
                t for row in records for t in row["successful_decode_seconds"]
            ]
            encode_seconds = [float(row["encode_seconds"]) for row in records]
            rss = [float(row["peak_rss_mib"]) for row in records]
            summary[version][condition] = {
                "recorded_trial_observations": trials * repeats,
                "independent_trial_seeds": trials,
                "exact_recoveries": sum(int(row["successes"]) for row in records),
                "failed_trials": sum(int(row["failures"]) for row in records),
                "timeout_trials": sum(int(row["timeouts"]) for row in records),
                "successful_decode_seconds": completed_seconds,
                "median_successful_decode_seconds": (
                    statistics.median(completed_seconds) if completed_seconds else None
                ),
                "encode_seconds_per_repeat": encode_seconds,
                "peak_worker_rss_mib_per_repeat": rss,
            }
    return {
        "method": method,
        "baseline_sha": baseline_sha,
        "optimized_sha": optimized_sha,
        "payload_sha256": next(iter(payload_hashes)),
        "size_bytes": SIZE,
        "trial_seeds": [20_260_000 + trial for trial in range(trials)],
        "per_trial_deadline_seconds": DEADLINE,
        "conditions": list(CONDITIONS),
        "order_by_repeat": order,
        "unprofiled": True,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "results": summary,
        "interpretation": (
            "Research-only paired unprofiled measurements. Failed trials and timeouts remain "
            "in denominators. Repeated seeds are paired timing observations, not "
            "independent recovery trials. Successful-only decode medians exclude "
            "unsuccessful trials and cannot establish an overall speedup. No statistical "
            "significance or end-to-end performance improvement is claimed from this "
            "small sample. Peak RSS is per worker, not a matched allocation microbenchmark."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--method", required=True)
    parser.add_argument("--trials", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--baseline-sha", required=True)
    parser.add_argument("--optimized-sha", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = summarize(
        args.root,
        method=args.method,
        trials=args.trials,
        repeats=args.repeats,
        baseline_sha=args.baseline_sha,
        optimized_sha=args.optimized_sha,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

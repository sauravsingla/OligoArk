"""Summarize timeout call-stack hotspots in matched-codec research JSON.

Usage: python benchmarks/summarize_codec_timeouts.py matched-codec-results/results.json
The report is diagnostic only: a deadline stack is one sample per timed-out
trial, not a statistical profiler or proof of a performance improvement.
"""
from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path
from typing import Any


def summarize(payload: Any) -> dict[str, object]:
    """Count all failures and deadline-stack locations without dropping trials."""
    if isinstance(payload, dict):
        rows = payload.get("rows", payload.get("results", []))
    else:
        rows = payload
    if not isinstance(rows, list):
        raise ValueError("expected a list of result rows")
    groups: dict[str, dict[str, object]] = {}
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("invalid result row")
        key = f"{row.get('method', 'unknown')} / {row.get('condition', 'unknown')}"
        group = groups.setdefault(key, {
            "trials": 0, "successes": 0, "timeouts": 0,
            "timeout_stages": collections.Counter(),
            "interrupted_locations": collections.Counter(),
        })
        trials = row.get("trial_results", [])
        if not isinstance(trials, list):
            raise ValueError("invalid trial_results")
        for trial in trials:
            if not isinstance(trial, dict):
                raise ValueError("invalid trial")
            group["trials"] += 1
            group["successes"] += int(trial.get("sha256_verified", False))
            if not trial.get("timed_out", False):
                continue
            group["timeouts"] += 1
            group["timeout_stages"][str(trial.get("timeout_stage", "unknown"))] += 1
            stack = trial.get("timeout_stack", [])
            if stack:
                frame = stack[-1]
                location = f"{frame.get('file', '?')}:{frame.get('line', '?')}:{frame.get('function', '?')}"
                group["interrupted_locations"][location] += 1
    return {
        "note": "Single interrupted frames are diagnostic, not representative CPU profiles.",
        "groups": {
            key: {
                "trials": group["trials"],
                "successes": group["successes"],
                "timeouts": group["timeouts"],
                "timeout_stages": dict(group["timeout_stages"].most_common()),
                "interrupted_locations": dict(group["interrupted_locations"].most_common()),
            }
            for key, group in sorted(groups.items())
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = summarize(json.loads(args.results.read_text()))
    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()

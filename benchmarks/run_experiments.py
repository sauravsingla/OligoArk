"""Run OligoArk multi-seed ablation experiments and save reproducible artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import platform
from dataclasses import asdict
from pathlib import Path

from oligoark import __version__
from oligoark.experiments import (
    aggregate_experiments,
    publication_profile,
    run_experiments,
    smoke_profile,
)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("cannot write empty experiment CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_plots(summary: list[dict[str, object]], output: Path) -> None:
    import matplotlib.pyplot as plt

    labels = [
        f"{row['scenario']}\n{row['strategy']}\n{row['payload_size']}B"
        for row in summary
    ]
    recovery = [float(row["recovery_rate"]) for row in summary]
    overhead = [float(row["mean_overhead_ratio"]) for row in summary]

    plt.figure(figsize=(max(12, len(labels) * 0.45), 5))
    plt.bar(labels, recovery)
    plt.ylabel("SHA-256 verified recovery rate")
    plt.ylim(0, 1.05)
    plt.xticks(rotation=65, ha="right")
    plt.tight_layout()
    plt.savefig(output / "experiment_recovery.png", dpi=160)
    plt.close()

    plt.figure(figsize=(max(12, len(labels) * 0.45), 5))
    plt.bar(labels, overhead)
    plt.ylabel("Encoded nt / ideal 2-bit nt")
    plt.xticks(rotation=65, ha="right")
    plt.tight_layout()
    plt.savefig(output / "experiment_overhead.png", dpi=160)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("smoke", "publication"), default="smoke")
    parser.add_argument("--output-dir", default="experiment-results")
    args = parser.parse_args()

    profile = smoke_profile() if args.profile == "smoke" else publication_profile()
    records = run_experiments(profile)
    summaries = aggregate_experiments(records)

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    raw_rows = [record.to_dict() for record in records]
    summary_rows = [summary.to_dict() for summary in summaries]
    (output / "raw.json").write_text(json.dumps(raw_rows, indent=2), encoding="utf-8")
    (output / "summary.json").write_text(json.dumps(summary_rows, indent=2), encoding="utf-8")
    _write_csv(output / "raw.csv", raw_rows)
    _write_csv(output / "summary.csv", summary_rows)

    metadata = {
        "oligoark_version": __version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "profile": args.profile,
        "seeds": list(profile.seeds),
        "payload_sizes": list(profile.payload_sizes),
        "strategies": list(profile.strategies),
        "scenarios": [asdict(scenario) for scenario in profile.scenarios],
        "claim_scope": "software simulation only; no wet-lab performance is implied",
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    _write_plots(summary_rows, output)
    print(json.dumps({"metadata": metadata, "summary": summary_rows}, indent=2))


if __name__ == "__main__":
    main()

"""Run held-out OligoArk ablation experiments and save reproducible research artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
from dataclasses import asdict, replace
from pathlib import Path

from oligoark import __version__
from oligoark.experiments import (
    CalibrationRecord,
    aggregate_experiments,
    paired_strategy_effects,
    publication_profile,
    run_experiment_bundle,
    smoke_profile,
)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("cannot write empty experiment CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _calibration_candidate_rows(
    calibrations: list[CalibrationRecord],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for calibration in calibrations:
        optimization = calibration.optimization
        selected_key = (
            optimization.best_config,
            optimization.reconstruction_mode,
            optimization.best_score,
        )
        for candidate_index, evaluation in enumerate(optimization.evaluations):
            config = asdict(evaluation.config)
            objective = asdict(evaluation.objective)
            row: dict[str, object] = {
                "scenario": calibration.scenario,
                "payload_size": calibration.payload_size,
                "calibration_payload_size": calibration.calibration_payload_size,
                "calibration_seeds": ",".join(map(str, optimization.calibration_seeds)),
                "search_method": optimization.search_method,
                "search_seed": optimization.search_seed,
                "total_possible_candidates": optimization.total_possible_candidates,
                "evaluated_candidates": optimization.evaluated_candidates,
                "rejected_candidates": optimization.rejected_candidates,
                "candidate_index": candidate_index,
                "selected": (
                    evaluation.config,
                    evaluation.reconstruction_mode,
                    evaluation.score,
                )
                == selected_key,
                "reconstruction_mode": evaluation.reconstruction_mode,
                "trials": evaluation.trials,
                "verified_successes": evaluation.verified_successes,
                "recovery_rate": evaluation.recovery_rate,
                "fold_recovery_rates": ",".join(
                    f"{value:.6f}" for value in evaluation.fold_recovery_rates
                ),
                "recovery_instability": evaluation.recovery_instability,
                "encoded_nucleotides": evaluation.encoded_nucleotides,
                "overhead_ratio": evaluation.overhead_ratio,
                "redundancy_ratio": evaluation.redundancy_ratio,
                "mean_runtime_seconds": evaluation.mean_runtime_seconds,
                "graph_recovery_count": evaluation.graph_recovery_count,
                "score": evaluation.score,
                "rejected_reason": evaluation.rejected_reason or "",
            }
            row.update({f"config_{key}": value for key, value in config.items()})
            row.update({f"objective_{key}": value for key, value in objective.items()})
            rows.append(row)
    return rows


def _error_rate(row: dict[str, object]) -> float:
    return sum(
        float(row[name])
        for name in (
            "substitution_rate",
            "insertion_rate",
            "deletion_rate",
            "dropout_rate",
        )
    )


def write_plots(
    raw: list[dict[str, object]],
    summary: list[dict[str, object]],
    effects: list[dict[str, object]],
    output: Path,
) -> None:
    import matplotlib.pyplot as plt

    labels = [
        f"{row['scenario']}\n{row['strategy']}\n{row['payload_size']}B"
        for row in summary
    ]
    recovery = [float(row["recovery_rate"]) for row in summary]
    overhead = [float(row["mean_overhead_ratio"]) for row in summary]
    runtime = [float(row["mean_runtime_seconds"]) for row in summary]
    rescue = [float(row["graph_rescue_rate"]) for row in summary]

    def save_bar(
        values: list[float],
        ylabel: str,
        filename: str,
        ylim: tuple[float, float] | None = None,
    ) -> None:
        plt.figure(figsize=(max(12, len(labels) * 0.45), 5))
        plt.bar(labels, values)
        plt.ylabel(ylabel)
        if ylim is not None:
            plt.ylim(*ylim)
        plt.xticks(rotation=65, ha="right")
        plt.tight_layout()
        plt.savefig(output / filename, dpi=160)
        plt.close()

    save_bar(recovery, "SHA-256 verified recovery rate", "experiment_recovery.png", (0, 1.05))
    save_bar(overhead, "Encoded nt / ideal 2-bit nt", "experiment_overhead.png")
    save_bar(runtime, "Mean runtime (seconds)", "runtime_by_strategy.png")
    save_bar(rescue, "Direct-failure rescue rate", "graph_rescue_rate.png", (0, 1.05))

    plt.figure(figsize=(8, 5))
    for strategy in sorted({str(row["strategy"]) for row in raw}):
        subset = [row for row in raw if row["strategy"] == strategy]
        plt.scatter(
            [_error_rate(row) for row in subset],
            [1.0 if bool(row["recovered"]) else 0.0 for row in subset],
            label=strategy,
            alpha=0.6,
        )
    plt.xlabel("Configured aggregate software error/dropout rate")
    plt.ylabel("SHA-256 verified recovery")
    plt.legend(fontsize="small")
    plt.tight_layout()
    plt.savefig(output / "recovery_vs_error_rate.png", dpi=160)
    plt.close()

    plt.figure(figsize=(8, 5))
    for strategy in sorted({str(row["strategy"]) for row in summary}):
        subset = [row for row in summary if row["strategy"] == strategy]
        plt.scatter(
            [float(row["mean_overhead_ratio"]) for row in subset],
            [float(row["recovery_rate"]) for row in subset],
            label=strategy,
        )
    plt.xlabel("Mean encoded nucleotide overhead ratio")
    plt.ylabel("Recovery rate")
    plt.legend(fontsize="small")
    plt.tight_layout()
    plt.savefig(output / "overhead_vs_recovery.png", dpi=160)
    plt.close()

    plt.figure(figsize=(8, 5))
    for strategy in sorted({str(row["strategy"]) for row in summary}):
        subset = [row for row in summary if row["strategy"] == strategy]
        plt.scatter(
            [float(row["mean_runtime_seconds"]) for row in subset],
            [float(row["recovery_rate"]) for row in subset],
            label=strategy,
        )
    plt.xlabel("Mean runtime (seconds)")
    plt.ylabel("Recovery rate")
    plt.legend(fontsize="small")
    plt.tight_layout()
    plt.savefig(output / "runtime_vs_recovery.png", dpi=160)
    plt.close()

    effect_labels = [
        f"{row['scenario']}\n{row['strategy']}\n{row['payload_size']}B"
        for row in effects
    ]
    effect_values = [float(row["recovery_rate_difference"]) for row in effects]
    plt.figure(figsize=(max(12, len(effect_labels) * 0.45), 5))
    plt.bar(effect_labels, effect_values)
    plt.axhline(0.0)
    plt.ylabel("Paired recovery-rate difference vs fixed")
    plt.xticks(rotation=65, ha="right")
    plt.tight_layout()
    plt.savefig(output / "strategy_ablation.png", dpi=160)
    plt.close()

    combined = [row for row in raw if row["strategy"] == "combined"]
    generalization: dict[tuple[str, int], list[dict[str, object]]] = {}
    for row in combined:
        generalization.setdefault((str(row["scenario"]), int(row["payload_size"])), []).append(row)
    gen_labels: list[str] = []
    gen_calibration: list[float] = []
    gen_heldout: list[float] = []
    for (scenario, payload_size), rows in sorted(generalization.items()):
        scores = [
            float(row["selection_calibration_recovery_rate"])
            for row in rows
            if row["selection_calibration_recovery_rate"] is not None
        ]
        gen_labels.append(f"{scenario}\n{payload_size}B")
        gen_calibration.append(sum(scores) / len(scores) if scores else 0.0)
        gen_heldout.append(sum(float(bool(row["recovered"])) for row in rows) / len(rows))
    positions = list(range(len(gen_labels)))
    width = 0.4
    plt.figure(figsize=(max(10, len(gen_labels) * 0.55), 5))
    plt.bar(
        [position - width / 2 for position in positions],
        gen_calibration,
        width,
        label="calibration recovery",
    )
    plt.bar(
        [position + width / 2 for position in positions],
        gen_heldout,
        width,
        label="held-out recovery",
    )
    plt.ylabel("SHA-256 verified recovery rate")
    plt.xticks(positions, gen_labels, rotation=65, ha="right")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output / "optimizer_generalization.png", dpi=160)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("smoke", "publication"), default="smoke")
    parser.add_argument("--output-dir", default="experiment-results")
    parser.add_argument("--payload-size", type=int)
    parser.add_argument("--scenario")
    parser.add_argument("--strategies")
    parser.add_argument("--seed-shard-index", type=int)
    parser.add_argument("--seed-shard-count", type=int)
    args = parser.parse_args()

    profile = smoke_profile() if args.profile == "smoke" else publication_profile()
    if args.payload_size is not None:
        if args.payload_size not in profile.payload_sizes:
            raise ValueError("requested payload size is not part of the selected profile")
        profile = replace(profile, payload_sizes=(args.payload_size,))
    if args.scenario is not None:
        selected = tuple(
            scenario for scenario in profile.scenarios if scenario.name == args.scenario
        )
        if not selected:
            raise ValueError("requested scenario is not part of the selected profile")
        profile = replace(profile, scenarios=selected)
    if args.strategies is not None:
        requested_strategies = tuple(
            item.strip() for item in args.strategies.split(",") if item.strip()
        )
        if not requested_strategies:
            raise ValueError("at least one strategy must be requested")
        unknown = sorted(set(requested_strategies) - set(profile.strategies))
        if unknown:
            raise ValueError(f"unknown requested strategy or strategies: {unknown}")
        profile = replace(profile, strategies=requested_strategies)
    if (args.seed_shard_index is None) != (args.seed_shard_count is None):
        raise ValueError("seed shard index and count must be supplied together")
    if args.seed_shard_count is not None:
        if args.seed_shard_count < 1:
            raise ValueError("seed shard count must be positive")
        if args.seed_shard_index is None or not 0 <= args.seed_shard_index < args.seed_shard_count:
            raise ValueError("seed shard index must be within the shard count")
        selected_seeds = tuple(
            seed
            for index, seed in enumerate(profile.evaluation_seeds)
            if index % args.seed_shard_count == args.seed_shard_index
        )
        if not selected_seeds:
            raise ValueError("seed shard selected no evaluation seeds")
        profile = replace(profile, seeds=selected_seeds)
    bundle = run_experiment_bundle(profile)
    records = list(bundle.records)
    calibrations = list(bundle.calibrations)
    summaries = aggregate_experiments(records)
    effects = paired_strategy_effects(records)

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    raw_rows = [record.to_dict() for record in records]
    summary_rows = [summary.to_dict() for summary in summaries]
    effect_rows = [effect.to_dict() for effect in effects]
    calibration_rows = [calibration.to_dict() for calibration in calibrations]
    calibration_candidate_rows = _calibration_candidate_rows(calibrations)
    (output / "raw.json").write_text(json.dumps(raw_rows, indent=2), encoding="utf-8")
    (output / "summary.json").write_text(json.dumps(summary_rows, indent=2), encoding="utf-8")
    (output / "paired-effects.json").write_text(
        json.dumps(effect_rows, indent=2), encoding="utf-8"
    )
    if calibration_rows:
        (output / "calibration.json").write_text(
            json.dumps(calibration_rows, indent=2), encoding="utf-8"
        )
    _write_csv(output / "raw.csv", raw_rows)
    _write_csv(output / "summary.csv", summary_rows)
    if effect_rows:
        _write_csv(output / "paired-effects.csv", effect_rows)
    if calibration_candidate_rows:
        _write_csv(output / "calibration-candidates.csv", calibration_candidate_rows)

    metadata = {
        "oligoark_version": __version__,
        "git_commit": os.environ.get("GITHUB_SHA", "unknown-local"),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "profile": args.profile,
        "calibration_seeds": list(profile.calibration_seeds),
        "evaluation_seeds": list(profile.evaluation_seeds),
        "seed_shard_index": args.seed_shard_index,
        "seed_shard_count": args.seed_shard_count,
        "seed_sets_disjoint": not bool(
            set(profile.calibration_seeds) & set(profile.evaluation_seeds)
        ),
        "payload_sizes": list(profile.payload_sizes),
        "strategies": list(profile.strategies),
        "optimizer_search_method": profile.optimizer_search_method,
        "optimizer_search_seed": profile.optimizer_search_seed,
        "optimizer_max_candidates": profile.optimizer_max_candidates,
        "calibration_payload_limit_bytes": max(
            (item.calibration_payload_size for item in calibrations),
            default=0,
        ),
        "calibration_payload_variants": 3,
        "calibration_record_count": len(calibrations),
        "scenarios": [asdict(scenario) for scenario in profile.scenarios],
        "claim_scope": "software simulation only; no wet-lab performance is implied",
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    write_plots(raw_rows, summary_rows, effect_rows, output)
    print(
        json.dumps(
            {
                "metadata": metadata,
                "summary": summary_rows,
                "paired_effects": effect_rows,
                "calibrations": calibration_rows,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

"""Aggregate deterministic publication shards into one analysis artifact."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import cast

from run_experiments import write_plots

from oligoark.experiments import (
    ExperimentRecord,
    aggregate_experiments,
    paired_strategy_effects,
)


def _record(raw: dict[str, object]) -> ExperimentRecord:
    return ExperimentRecord(
        strategy=str(raw["strategy"]),
        scenario=str(raw["scenario"]),
        seed=int(cast(int, raw["seed"])),
        payload_size=int(cast(int, raw["payload_size"])),
        payload_sha256=str(raw["payload_sha256"]),
        recovered=bool(raw["recovered"]),
        encoded_nucleotides=int(cast(int, raw["encoded_nucleotides"])),
        overhead_ratio=float(cast(float, raw["overhead_ratio"])),
        strand_count=int(cast(int, raw["strand_count"])),
        read_count=int(cast(int, raw["read_count"])),
        runtime_seconds=float(cast(float, raw["runtime_seconds"])),
        graph_reconstruction_used=bool(raw["graph_reconstruction_used"]),
        redundancy_scheme=str(raw["redundancy_scheme"]),
        rs_nsym=int(cast(int, raw["rs_nsym"])),
        chunk_size=int(cast(int, raw["chunk_size"])),
        parity_group_size=int(cast(int, raw["parity_group_size"])),
        adaptive_masks=bool(raw["adaptive_masks"]),
        substitution_rate=float(cast(float, raw["substitution_rate"])),
        insertion_rate=float(cast(float, raw["insertion_rate"])),
        deletion_rate=float(cast(float, raw["deletion_rate"])),
        dropout_rate=float(cast(float, raw["dropout_rate"])),
        duplicate_rate=float(cast(float, raw["duplicate_rate"])),
        calibration_seeds=tuple(
            int(value) for value in cast(list[int], raw.get("calibration_seeds", []))
        ),
        selection_score=(
            float(cast(float, raw["selection_score"]))
            if raw.get("selection_score") is not None
            else None
        ),
        selection_search_method=(
            str(raw["selection_search_method"])
            if raw.get("selection_search_method") is not None
            else None
        ),
        reconstruction_mode=str(raw.get("reconstruction_mode", "direct")),
        copies_per_strand=int(cast(int, raw.get("copies_per_strand", 1))),
    )


def _calibration_signature(row: dict[str, object]) -> str:
    optimization = cast(dict[str, object], row["optimization"])
    evaluations = cast(list[dict[str, object]], optimization["evaluations"])
    stable = {
        "scenario": row["scenario"],
        "payload_size": row["payload_size"],
        "calibration_payload_size": row["calibration_payload_size"],
        "best_config": optimization["best_config"],
        "reconstruction_mode": optimization["reconstruction_mode"],
        "best_score": optimization["best_score"],
        "search_method": optimization["search_method"],
        "search_seed": optimization["search_seed"],
        "calibration_seeds": optimization["calibration_seeds"],
        "evaluations": [
            {
                "config": item["config"],
                "reconstruction_mode": item["reconstruction_mode"],
                "trials": item["trials"],
                "verified_successes": item["verified_successes"],
                "recovery_rate": item["recovery_rate"],
                "fold_recovery_rates": item["fold_recovery_rates"],
                "recovery_instability": item["recovery_instability"],
                "encoded_nucleotides": item["encoded_nucleotides"],
                "overhead_ratio": item["overhead_ratio"],
                "redundancy_ratio": item["redundancy_ratio"],
                "graph_recovery_count": item["graph_recovery_count"],
                "score": item["score"],
                "rejected_reason": item["rejected_reason"],
            }
            for item in evaluations
        ],
    }
    return json.dumps(stable, sort_keys=True)


def _candidate_signature(row: dict[str, object]) -> tuple[tuple[str, str], ...]:
    ignored = {
        "mean_runtime_seconds",
        "objective_runtime_penalty",
        "objective_retrieval_penalty",
    }
    return tuple(
        sorted(
            (key, str(value))
            for key, value in row.items()
            if key not in ignored
        )
    )


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    root = Path("publication-shards")
    raw_paths = sorted(root.glob("**/raw.json"))
    metadata_paths = sorted(root.glob("**/metadata.json"))
    calibration_paths = sorted(root.glob("**/calibration.json"))
    candidate_paths = sorted(root.glob("**/calibration-candidates.csv"))
    if not raw_paths:
        raise ValueError("no publication shard raw.json files were found")

    raw_rows: list[dict[str, object]] = []
    for path in raw_paths:
        value = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(value, list):
            raise ValueError(f"{path} must contain a list")
        raw_rows.extend(cast(list[dict[str, object]], value))

    calibration_by_key: dict[tuple[str, int], dict[str, object]] = {}
    for path in calibration_paths:
        value = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(value, list):
            raise ValueError(f"{path} must contain a list")
        for row in cast(list[dict[str, object]], value):
            key = (str(row["scenario"]), int(cast(int, row["payload_size"])))
            existing = calibration_by_key.get(key)
            if (
                existing is not None
                and _calibration_signature(existing) != _calibration_signature(row)
            ):
                raise ValueError(f"calibration shards disagree for {key}")
            if existing is None:
                calibration_by_key[key] = row
    calibration_rows = [
        calibration_by_key[key] for key in sorted(calibration_by_key)
    ]

    candidate_by_key: dict[tuple[str, int, str], dict[str, object]] = {}
    for path in candidate_paths:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                key = (
                    str(row["scenario"]),
                    int(row["payload_size"]),
                    str(row["candidate_index"]),
                )
                current = dict(row)
                existing = candidate_by_key.get(key)
                if (
                    existing is not None
                    and _candidate_signature(existing) != _candidate_signature(current)
                ):
                    raise ValueError(f"calibration candidate shards disagree for {key}")
                if existing is None:
                    candidate_by_key[key] = current
    candidate_rows = [
        candidate_by_key[key] for key in sorted(candidate_by_key)
    ]

    records = [_record(row) for row in raw_rows]
    summaries = aggregate_experiments(records)
    effects = paired_strategy_effects(records)
    summary_rows = [item.to_dict() for item in summaries]
    effect_rows = [item.to_dict() for item in effects]

    metadata_values = [
        json.loads(path.read_text(encoding="utf-8")) for path in metadata_paths
    ]
    if not metadata_values:
        raise ValueError("publication shard metadata is missing")
    first = metadata_values[0]
    payload_sizes = sorted(
        {
            int(size)
            for metadata in metadata_values
            for size in metadata.get("payload_sizes", [])
        }
    )
    evaluation_seeds = sorted(
        {
            int(seed)
            for metadata in metadata_values
            for seed in metadata.get("evaluation_seeds", [])
        }
    )
    calibration_seed_sets = {
        tuple(int(seed) for seed in metadata.get("calibration_seeds", []))
        for metadata in metadata_values
    }
    if len(calibration_seed_sets) != 1:
        raise ValueError("publication shards disagree on calibration seeds")
    strategy_order = (
        "fixed",
        "adaptive",
        "adaptive_fountain",
        "adaptive_medoid",
        "adaptive_graph",
        "adaptive_trace",
        "combined",
    )
    observed_strategies = {
        str(strategy)
        for metadata in metadata_values
        for strategy in metadata.get("strategies", [])
    }
    strategies = [
        strategy for strategy in strategy_order if strategy in observed_strategies
    ]
    scenario_by_name: dict[str, object] = {}
    for metadata in metadata_values:
        for scenario in metadata.get("scenarios", []):
            if isinstance(scenario, dict) and isinstance(scenario.get("name"), str):
                scenario_by_name[str(scenario["name"])] = scenario
    calibration_payload_limit = max(
        int(metadata.get("calibration_payload_limit_bytes", 0))
        for metadata in metadata_values
    )
    calibration_payload_variants = max(
        int(metadata.get("calibration_payload_variants", 0))
        for metadata in metadata_values
    )
    metadata = {
        **first,
        "payload_sizes": payload_sizes,
        "strategies": strategies,
        "calibration_payload_limit_bytes": calibration_payload_limit,
        "calibration_payload_variants": calibration_payload_variants,
        "evaluation_seeds": evaluation_seeds,
        "calibration_seeds": list(next(iter(calibration_seed_sets))),
        "seed_shard_index": None,
        "seed_shard_count": max(
            (
                int(metadata["seed_shard_count"])
                for metadata in metadata_values
                if metadata.get("seed_shard_count") is not None
            ),
            default=1,
        ),
        "scenarios": [
            scenario_by_name[name] for name in sorted(scenario_by_name)
        ],
        "shard_count": len(raw_paths),
        "raw_trial_count": len(records),
        "calibration_record_count": len(calibration_rows),
        "calibration_candidate_count": len(candidate_rows),
        "aggregation": (
            "merged from payload/scenario/evaluation-seed shards without dropping failures"
        ),
    }

    output = Path("publication-results")
    output.mkdir(parents=True, exist_ok=True)
    (output / "raw.json").write_text(json.dumps(raw_rows, indent=2), encoding="utf-8")
    (output / "summary.json").write_text(
        json.dumps(summary_rows, indent=2), encoding="utf-8"
    )
    (output / "paired-effects.json").write_text(
        json.dumps(effect_rows, indent=2), encoding="utf-8"
    )
    (output / "metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    (output / "calibration.json").write_text(
        json.dumps(calibration_rows, indent=2), encoding="utf-8"
    )
    _write_csv(output / "raw.csv", raw_rows)
    _write_csv(output / "summary.csv", summary_rows)
    _write_csv(output / "paired-effects.csv", effect_rows)
    if candidate_rows:
        _write_csv(output / "calibration-candidates.csv", candidate_rows)
    write_plots(raw_rows, summary_rows, effect_rows, output)
    print(
        json.dumps(
            {
                "shards": len(raw_paths),
                "raw_trials": len(records),
                "summaries": len(summary_rows),
                "effects": len(effect_rows),
                "calibrations": len(calibration_rows),
                "calibration_candidates": len(candidate_rows),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

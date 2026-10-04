"""Train transparent policy models from experiment output and evaluate held-out seeds."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import cast

from oligoark.experiments import ExperimentRecord
from oligoark.learning_eval import evaluate_learning_from_records


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-dir", default="experiment-results")
    parser.add_argument("--output-dir", default="learning-results")
    args = parser.parse_args()

    experiment_dir = Path(args.experiment_dir)
    raw = json.loads((experiment_dir / "raw.json").read_text(encoding="utf-8"))
    metadata = json.loads((experiment_dir / "metadata.json").read_text(encoding="utf-8"))
    if not isinstance(raw, list) or not all(isinstance(row, dict) for row in raw):
        raise ValueError("experiment raw.json must contain a list of objects")
    evaluation_seeds = tuple(int(seed) for seed in metadata["evaluation_seeds"])
    if len(evaluation_seeds) < 2:
        raise ValueError("learning evaluation requires at least two held-out experiment seeds")
    if len(evaluation_seeds) >= 6:
        training_count = max(2, (len(evaluation_seeds) * 2) // 5)
        validation_count = max(1, len(evaluation_seeds) // 5)
        training_seeds = evaluation_seeds[:training_count]
        validation_seeds = evaluation_seeds[
            training_count : training_count + validation_count
        ]
        test_seeds = evaluation_seeds[training_count + validation_count :]
    else:
        split = max(1, len(evaluation_seeds) // 2)
        training_seeds = evaluation_seeds[:split]
        validation_seeds = ()
        test_seeds = evaluation_seeds[split:]
    if not test_seeds:
        raise ValueError("learning evaluation needs non-empty test seeds")

    scenario_items = metadata.get("scenarios")
    if not isinstance(scenario_items, list):
        raise ValueError("experiment metadata must include scenario definitions")
    scenario_names = tuple(
        str(item["name"])
        for item in scenario_items
        if isinstance(item, dict) and "name" in item
    )
    if len(scenario_names) < 2:
        raise ValueError("learning evaluation needs at least two channel scenarios")
    scenario_split = max(1, len(scenario_names) // 2)
    training_scenarios = scenario_names[:scenario_split]
    test_scenarios = scenario_names[scenario_split:]
    if not test_scenarios:
        raise ValueError("learning evaluation needs non-empty held-out scenarios")

    records = [_record(cast(dict[str, object], row)) for row in raw]
    result = evaluate_learning_from_records(
        records,
        training_seeds=training_seeds,
        validation_seeds=validation_seeds,
        test_seeds=test_seeds,
        training_scenarios=training_scenarios,
        test_scenarios=test_scenarios,
    )
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    (output / "evaluation.json").write_text(
        json.dumps(result.to_dict(), indent=2), encoding="utf-8"
    )
    summary_rows = [summary.to_dict() for summary in result.summaries]
    with (output / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    (output / "linear-model.json").write_text(
        json.dumps(result.linear_model_state, indent=2), encoding="utf-8"
    )
    (output / "kernel-model.json").write_text(
        json.dumps(result.kernel_model_state, indent=2), encoding="utf-8"
    )
    print(json.dumps(result.to_dict(), indent=2))


if __name__ == "__main__":
    main()

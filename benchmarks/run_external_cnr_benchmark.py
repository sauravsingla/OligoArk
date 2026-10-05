"""Fast physical-read benchmark on Microsoft's Clustered Nanopore Reads dataset."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from oligoark.reconstruct import (
    GraphConsensusReconstructor,
    edit_distance,
    global_align,
    iterative_trace_consensus,
    medoid_consensus,
    multistart_trace_consensus,
)

DATASET_NAME = "Microsoft Clustered Nanopore Reads (CNR)"
DATASET_REPOSITORY = "microsoft/clustered-nanopore-reads-dataset"
DATASET_COMMIT = "6938f44796185902a08381943c2895782886c5c3"
CENTERS_GIT_BLOB = "2d72f78cd86f2172578aaae297113a4baa14975c"
CLUSTERS_GIT_BLOB = "95a71c77cb82f6d501d6aa99c7415591565c5ae3"
BBS_REPOSITORY = "GZHoffie/bbs"
BBS_COMMIT = "3e4ab46871929819e4f3e34a831c57cac88bb456"
TARGET_LENGTH = 110
DEFAULT_SUBSET_SIZE = 96
DEFAULT_COVERAGES = (1, 5, 10)
DEFAULT_SEED = 20261005
DEFAULT_CALIBRATION_SEED = 20261006
DEFAULT_CALIBRATION_SIZE = 48
CALIBRATION_COVERAGES = (5, 10)
MULTISTART_CANDIDATES: tuple[dict[str, object], ...] = (
    {
        "name": "a1-r2-bi",
        "anchors": 1,
        "rounds": 2,
        "bidirectional": True,
        "length_penalty": 1.0,
    },
    {
        "name": "a2-r1-bi",
        "anchors": 2,
        "rounds": 1,
        "bidirectional": True,
        "length_penalty": 1.0,
    },
    {
        "name": "a2-r2-bi",
        "anchors": 2,
        "rounds": 2,
        "bidirectional": True,
        "length_penalty": 1.0,
    },
    {
        "name": "a3-r1-bi",
        "anchors": 3,
        "rounds": 1,
        "bidirectional": True,
        "length_penalty": 1.0,
    },
    {
        "name": "a2-r2-forward",
        "anchors": 2,
        "rounds": 2,
        "bidirectional": False,
        "length_penalty": 1.0,
    },
)
OLIGOARK_METHODS = (
    "direct",
    "medoid",
    "graph_alignment",
    "iterative_trace",
    "multistart_trace",
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _validate_dna(sequence: str, *, expected_length: int | None = None) -> str:
    normalized = sequence.strip().upper()
    if not normalized or any(base not in "ACGT" for base in normalized):
        raise ValueError("expected a non-empty DNA sequence over A/C/G/T")
    if expected_length is not None and len(normalized) != expected_length:
        raise ValueError(
            f"expected DNA sequence of length {expected_length}, got {len(normalized)}"
        )
    return normalized


def load_centers(path: Path) -> list[str]:
    centers = [
        _validate_dna(line, expected_length=TARGET_LENGTH)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(centers) != 10_000:
        raise ValueError(f"expected 10000 CNR centers, got {len(centers)}")
    return centers


def load_clusters(path: Path, *, expected_count: int) -> list[list[str]]:
    clusters: list[list[str]] = []
    current: list[str] = []
    started = False
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("==="):
                if started:
                    clusters.append(current)
                    current = []
                else:
                    started = True
                continue
            if line:
                current.append(_validate_dna(line))
                started = True
    if current or len(clusters) < expected_count:
        clusters.append(current)
    if len(clusters) != expected_count:
        raise ValueError(
            f"expected {expected_count} CNR clusters, parsed {len(clusters)}"
        )
    return clusters


def _rank_key(seed: int, *values: object) -> bytes:
    joined = ":".join([str(seed), *(str(value) for value in values)])
    return hashlib.sha256(joined.encode("utf-8")).digest()


def select_subset(
    centers: list[str],
    clusters: list[list[str]],
    *,
    subset_size: int,
    max_coverage: int,
    seed: int,
    excluded_indices: set[int] | None = None,
) -> list[dict[str, Any]]:
    excluded = excluded_indices or set()
    eligible = [
        index
        for index, cluster in enumerate(clusters)
        if len(cluster) >= max_coverage and index not in excluded
    ]
    if len(eligible) < subset_size:
        raise ValueError(
            f"only {len(eligible)} clusters have at least {max_coverage} reads; "
            f"cannot select {subset_size}"
        )
    selected = sorted(eligible, key=lambda index: _rank_key(seed, index))[:subset_size]
    records: list[dict[str, Any]] = []
    for index in selected:
        ranked_reads = sorted(
            enumerate(clusters[index]),
            key=lambda item: _rank_key(seed, index, item[0], item[1]),
        )
        records.append(
            {
                "cluster_index": index + 1,
                "reference": centers[index],
                "available_reads": len(clusters[index]),
                "ranked_reads": [read for _, read in ranked_reads],
            }
        )
    return records


def _records_for_coverage(
    records: list[dict[str, Any]], coverage: int
) -> list[dict[str, Any]]:
    return [
        {
            "cluster_index": record["cluster_index"],
            "reference": record["reference"],
            "reads": record["ranked_reads"][:coverage],
        }
        for record in records
    ]


def write_microsoft_clusters(path: Path, records: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write("========================================\n")
            for read in record["reads"]:
                handle.write(f"{read}\n")


def _choose_graph_candidate(reads: list[str]) -> str:
    result = GraphConsensusReconstructor().reconstruct(reads)
    if not result.consensus_reads:
        return ""
    best_index = max(
        range(len(result.consensus_reads)),
        key=lambda index: (result.cluster_sizes[index], -index),
    )
    return result.consensus_reads[best_index]


def _reconstruct(
    method: str,
    reads: list[str],
    consensus_config: dict[str, object] | None = None,
) -> str:
    if method == "direct":
        return reads[0]
    if method == "medoid":
        return medoid_consensus(reads)
    if method == "graph_alignment":
        return _choose_graph_candidate(reads)
    if method == "iterative_trace":
        return iterative_trace_consensus(reads, rounds=3)
    if method == "multistart_trace":
        if consensus_config is None:
            raise ValueError("multistart_trace requires a frozen consensus_config")
        return multistart_trace_consensus(
            reads,
            target_length=TARGET_LENGTH,
            anchors=int(consensus_config["anchors"]),
            rounds=int(consensus_config["rounds"]),
            bidirectional=bool(consensus_config["bidirectional"]),
            length_penalty=float(consensus_config["length_penalty"]),
        )
    raise ValueError(f"unknown worker method: {method}")


def _case_record(
    method: str,
    coverage: int,
    cluster_index: int,
    reference: str,
    reconstruction: str,
) -> dict[str, Any]:
    distance = edit_distance(reference, reconstruction)
    return {
        "method": method,
        "method_family": "oligoark",
        "coverage": coverage,
        "cluster_index": cluster_index,
        "exact": reconstruction == reference,
        "edit_distance": distance,
        "normalized_edit_distance": round(distance / len(reference), 8),
        "reference": reference,
        "reconstruction": reconstruction,
        "reconstructed_length": len(reconstruction),
    }


def run_worker(method: str, subset_json: Path, output: Path) -> None:
    payload = json.loads(subset_json.read_text(encoding="utf-8"))
    coverage = int(payload["coverage"])
    consensus_config = payload.get("consensus_config")
    rows: list[dict[str, Any]] = []
    for record in payload["records"]:
        reads = [str(read) for read in record["reads"]]
        reconstruction = _reconstruct(method, reads, consensus_config)
        rows.append(
            _case_record(
                method,
                coverage,
                int(record["cluster_index"]),
                str(record["reference"]),
                reconstruction,
            )
        )
    output.write_text(json.dumps(rows, indent=2), encoding="utf-8")



def _evaluate_calibration_method(
    records: list[dict[str, Any]],
    *,
    coverages: tuple[int, ...],
    method: str,
    consensus_config: dict[str, object] | None = None,
) -> dict[str, Any]:
    """Evaluate a candidate only on the disjoint calibration references."""
    started = time.perf_counter()
    successes = 0
    distances: list[int] = []
    per_coverage: list[dict[str, Any]] = []
    for coverage in coverages:
        coverage_successes = 0
        coverage_distances: list[int] = []
        for record in _records_for_coverage(records, coverage):
            reads = [str(read) for read in record["reads"]]
            reference = str(record["reference"])
            reconstructed = _reconstruct(method, reads, consensus_config)
            distance = edit_distance(reference, reconstructed)
            coverage_successes += int(reconstructed == reference)
            coverage_distances.append(distance)
        successes += coverage_successes
        distances.extend(coverage_distances)
        per_coverage.append(
            {
                "coverage": coverage,
                "successes": coverage_successes,
                "trials": len(records),
                "mean_edit_distance": round(
                    statistics.fmean(coverage_distances),
                    8,
                ),
            }
        )
    elapsed = time.perf_counter() - started
    return {
        "method": method,
        "config": consensus_config,
        "successes": successes,
        "trials": len(records) * len(coverages),
        "exact_recovery_rate": round(
            successes / max(1, len(records) * len(coverages)),
            8,
        ),
        "mean_edit_distance": round(statistics.fmean(distances), 8),
        "runtime_seconds": round(elapsed, 6),
        "per_coverage": per_coverage,
    }


def _calibrate_multistart(
    records: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    """Choose the lightweight configuration before touching held-out references."""
    baseline = _evaluate_calibration_method(
        records,
        coverages=CALIBRATION_COVERAGES,
        method="iterative_trace",
    )
    candidates = [
        _evaluate_calibration_method(
            records,
            coverages=CALIBRATION_COVERAGES,
            method="multistart_trace",
            consensus_config=dict(config),
        )
        for config in MULTISTART_CANDIDATES
    ]
    runtime_budget = max(
        float(baseline["runtime_seconds"]) * 4.0,
        float(baseline["runtime_seconds"]) + 2.0,
    )
    practical = [
        candidate
        for candidate in candidates
        if float(candidate["runtime_seconds"]) <= runtime_budget
    ]
    if not practical:
        raise RuntimeError("all multistart calibration candidates exceeded the runtime budget")
    winner = min(
        practical,
        key=lambda candidate: (
            -int(candidate["successes"]),
            float(candidate["mean_edit_distance"]),
            float(candidate["runtime_seconds"]),
            str(candidate["config"]),
        ),
    )
    improved = (
        int(winner["successes"]) > int(baseline["successes"])
        or (
            int(winner["successes"]) == int(baseline["successes"])
            and float(winner["mean_edit_distance"])
            < float(baseline["mean_edit_distance"])
        )
    )
    if not improved:
        raise RuntimeError(
            "no lightweight multistart candidate improved calibration accuracy/edit distance"
        )
    return winner, baseline, candidates


def _timed_subprocess(command: list[str], timing_path: Path) -> tuple[float, float | None]:
    gnu_time = Path("/usr/bin/time")
    if gnu_time.exists():
        wrapped = [
            str(gnu_time),
            "-f",
            "%e,%M",
            "-o",
            str(timing_path),
            *command,
        ]
        subprocess.run(wrapped, check=True)
        elapsed_text, rss_text = timing_path.read_text(encoding="utf-8").strip().split(",")
        return float(elapsed_text), round(float(rss_text) / 1024.0, 4)
    started = time.perf_counter()
    subprocess.run(command, check=True)
    return round(time.perf_counter() - started, 6), None


def _parse_bbs_output(
    path: Path, records: list[dict[str, Any]], coverage: int
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        parsed = list(csv.DictReader(handle))
    if len(parsed) != len(records):
        raise ValueError(
            f"BBS returned {len(parsed)} reconstructed clusters for {len(records)} inputs"
        )
    for output_row, record in zip(parsed, records, strict=True):
        reconstruction = str(output_row["reconstruction_result"]).strip().upper()
        reference = str(record["reference"])
        distance = edit_distance(reference, reconstruction)
        rows.append(
            {
                "method": "external_bbs",
                "method_family": "external",
                "coverage": coverage,
                "cluster_index": int(record["cluster_index"]),
                "exact": reconstruction == reference,
                "edit_distance": distance,
                "normalized_edit_distance": round(distance / len(reference), 8),
                "reference": reference,
                "reconstruction": reconstruction,
                "reconstructed_length": len(reconstruction),
                "bbs_k": int(output_row["k"]),
                "bbs_path_weight": float(output_row["path_weight"]),
                "bbs_confidence": float(output_row["confidence"]),
            }
        )
    return rows


def wilson_interval(successes: int, trials: int) -> tuple[float, float]:
    if trials < 1:
        return (0.0, 0.0)
    z = 1.959963984540054
    p = successes / trials
    denominator = 1.0 + z * z / trials
    center = (p + z * z / (2 * trials)) / denominator
    spread = (
        z
        * math.sqrt(p * (1.0 - p) / trials + z * z / (4 * trials * trials))
        / denominator
    )
    return (max(0.0, center - spread), min(1.0, center + spread))


def mcnemar_exact(left: list[dict[str, Any]], right: list[dict[str, Any]]) -> dict[str, Any]:
    left_by_cluster = {int(row["cluster_index"]): bool(row["exact"]) for row in left}
    right_by_cluster = {int(row["cluster_index"]): bool(row["exact"]) for row in right}
    if left_by_cluster.keys() != right_by_cluster.keys():
        raise ValueError("paired methods must contain identical cluster IDs")
    left_only = sum(
        left_by_cluster[index] and not right_by_cluster[index] for index in left_by_cluster
    )
    right_only = sum(
        right_by_cluster[index] and not left_by_cluster[index] for index in left_by_cluster
    )
    discordant = left_only + right_only
    if discordant == 0:
        p_value = 1.0
    else:
        tail = sum(
            math.comb(discordant, value) for value in range(min(left_only, right_only) + 1)
        ) / (2**discordant)
        p_value = min(1.0, 2.0 * tail)
    return {
        "left_only_successes": left_only,
        "right_only_successes": right_only,
        "discordant_pairs": discordant,
        "two_sided_exact_p": p_value,
    }


def summarize_rows(
    rows: list[dict[str, Any]], elapsed_seconds: float, peak_rss_mb: float | None
) -> dict[str, Any]:
    successes = sum(bool(row["exact"]) for row in rows)
    trials = len(rows)
    low, high = wilson_interval(successes, trials)
    edit_distances = [int(row["edit_distance"]) for row in rows]
    normalized = [float(row["normalized_edit_distance"]) for row in rows]
    reference_concat = "".join(str(row["reference"]) for row in rows)
    reconstructed_concat = "".join(str(row["reconstruction"]) for row in rows)
    return {
        "method": rows[0]["method"],
        "method_family": rows[0]["method_family"],
        "coverage": int(rows[0]["coverage"]),
        "successes": successes,
        "trials": trials,
        "exact_recovery_rate": round(successes / trials, 8),
        "recovery_ci95_low": round(low, 8),
        "recovery_ci95_high": round(high, 8),
        "mean_edit_distance": round(statistics.fmean(edit_distances), 8),
        "median_edit_distance": round(float(statistics.median(edit_distances)), 8),
        "mean_normalized_edit_distance": round(statistics.fmean(normalized), 8),
        "elapsed_seconds": round(elapsed_seconds, 6),
        "seconds_per_cluster": round(elapsed_seconds / trials, 8),
        "peak_rss_mb": peak_rss_mb,
        "reference_subset_sha256": _sha256_text(reference_concat),
        "reconstructed_subset_sha256": _sha256_text(reconstructed_concat),
        "subset_sha256_exact_match": reference_concat == reconstructed_concat,
    }


def _error_profile(records: list[dict[str, Any]], max_coverage: int) -> dict[str, Any]:
    substitutions = 0
    insertions = 0
    deletions = 0
    bases = 0
    reads = 0
    for record in records:
        reference = str(record["reference"])
        for read in record["ranked_reads"][:max_coverage]:
            aligned_reference, aligned_read = global_align(reference, str(read))
            for reference_base, read_base in zip(aligned_reference, aligned_read, strict=True):
                if reference_base == "-":
                    insertions += 1
                elif read_base == "-":
                    deletions += 1
                    bases += 1
                else:
                    bases += 1
                    substitutions += int(reference_base != read_base)
            reads += 1
    denominator = max(1, bases)
    return {
        "reads_profiled": reads,
        "reference_bases_profiled": bases,
        "substitutions": substitutions,
        "insertions": insertions,
        "deletions": deletions,
        "substitution_rate_per_reference_base": round(substitutions / denominator, 8),
        "insertion_rate_per_reference_base": round(insertions / denominator, 8),
        "deletion_rate_per_reference_base": round(deletions / denominator, 8),
    }


def _git_commit() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _relation(
    oligoark_summary: dict[str, Any],
    bbs_summary: dict[str, Any],
    paired: dict[str, Any],
) -> str:
    left = float(oligoark_summary["exact_recovery_rate"])
    right = float(bbs_summary["exact_recovery_rate"])
    p_value = float(paired["two_sided_exact_p"])
    if p_value < 0.05 and left > right:
        return "statistically_outperforms_bbs_on_exact_recovery"
    if p_value < 0.05 and left < right:
        return "statistically_underperforms_bbs_on_exact_recovery"
    return "no_statistically_significant_exact_recovery_difference"


def run_benchmark(args: argparse.Namespace) -> None:
    if not args.centers or not args.clusters or not args.bbs_bin:
        raise ValueError("--centers, --clusters and --bbs-bin are required")
    centers_path = Path(args.centers)
    clusters_path = Path(args.clusters)
    bbs_bin = Path(args.bbs_bin)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    coverages = tuple(sorted(set(args.coverages)))
    if not coverages or coverages[0] < 1:
        raise ValueError("coverages must contain positive integers")
    max_coverage = max(coverages)

    centers = load_centers(centers_path)
    clusters = load_clusters(clusters_path, expected_count=len(centers))
    selected = select_subset(
        centers,
        clusters,
        subset_size=args.subset_size,
        max_coverage=max_coverage,
        seed=args.seed,
    )
    held_out_zero_based = {
        int(record["cluster_index"]) - 1
        for record in selected
    }
    calibration_records = select_subset(
        centers,
        clusters,
        subset_size=args.calibration_size,
        max_coverage=max_coverage,
        seed=args.calibration_seed,
        excluded_indices=held_out_zero_based,
    )
    calibration_ids = {
        int(record["cluster_index"])
        for record in calibration_records
    }
    held_out_ids = {
        int(record["cluster_index"])
        for record in selected
    }
    if calibration_ids & held_out_ids:
        raise RuntimeError("calibration and held-out cluster IDs must be disjoint")

    calibration_winner, calibration_baseline, calibration_candidates = (
        _calibrate_multistart(calibration_records)
    )
    frozen_multistart_config = dict(calibration_winner["config"])

    selected_metadata = [
        {
            "cluster_index": int(record["cluster_index"]),
            "available_reads": int(record["available_reads"]),
            "reference_sha256": _sha256_text(str(record["reference"])),
        }
        for record in selected
    ]
    calibration_metadata = [
        {
            "cluster_index": int(record["cluster_index"]),
            "available_reads": int(record["available_reads"]),
            "reference_sha256": _sha256_text(str(record["reference"])),
        }
        for record in calibration_records
    ]
    (output_dir / "selected-clusters.json").write_text(
        json.dumps(selected_metadata, indent=2), encoding="utf-8"
    )
    (output_dir / "calibration-clusters.json").write_text(
        json.dumps(calibration_metadata, indent=2), encoding="utf-8"
    )
    (output_dir / "calibration.json").write_text(
        json.dumps(
            {
                "baseline": calibration_baseline,
                "candidates": calibration_candidates,
                "winner": calibration_winner,
                "held_out_overlap_count": 0,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    all_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    bbs_repeat_summaries: list[dict[str, Any]] = []
    if args.bbs_repeats < 1:
        raise ValueError("--bbs-repeats must be positive")
    timing_dir = output_dir / "timing"
    timing_dir.mkdir(exist_ok=True)

    for coverage in coverages:
        coverage_records = _records_for_coverage(selected, coverage)
        subset_json = output_dir / f"subset-coverage-{coverage}.json"
        subset_json.write_text(
            json.dumps(
                {
                    "coverage": coverage,
                    "records": coverage_records,
                    "consensus_config": frozen_multistart_config,
                }
            ),
            encoding="utf-8",
        )
        clusters_subset = output_dir / f"clusters-coverage-{coverage}.txt"
        write_microsoft_clusters(clusters_subset, coverage_records)

        for method in OLIGOARK_METHODS:
            worker_output = output_dir / f"raw-{method}-coverage-{coverage}.json"
            elapsed, peak_rss = _timed_subprocess(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--worker-method",
                    method,
                    "--subset-json",
                    str(subset_json),
                    "--worker-output",
                    str(worker_output),
                ],
                timing_dir / f"{method}-coverage-{coverage}.txt",
            )
            rows = json.loads(worker_output.read_text(encoding="utf-8"))
            all_rows.extend(rows)
            summaries.append(summarize_rows(rows, elapsed, peak_rss))

        for repeat in range(1, args.bbs_repeats + 1):
            bbs_output = output_dir / f"bbs-coverage-{coverage}-repeat-{repeat}.csv"
            elapsed, peak_rss = _timed_subprocess(
                [
                    str(bbs_bin),
                    str(clusters_subset),
                    "-l",
                    str(TARGET_LENGTH),
                    "-t",
                    "1",
                    "-o",
                    str(bbs_output),
                ],
                timing_dir / f"external-bbs-coverage-{coverage}-repeat-{repeat}.txt",
            )
            bbs_rows = _parse_bbs_output(bbs_output, coverage_records, coverage)
            bbs_summary = summarize_rows(bbs_rows, elapsed, peak_rss)
            bbs_summary["repeat"] = repeat
            bbs_repeat_summaries.append(bbs_summary)
            if repeat == 1:
                all_rows.extend(bbs_rows)
                summaries.append(bbs_summary)

    comparisons: list[dict[str, Any]] = []
    for coverage in coverages:
        bbs_rows = [
            row
            for row in all_rows
            if row["method"] == "external_bbs" and row["coverage"] == coverage
        ]
        bbs_summary = next(
            row
            for row in summaries
            if row["method"] == "external_bbs" and row["coverage"] == coverage
        )
        for method in OLIGOARK_METHODS:
            method_rows = [
                row
                for row in all_rows
                if row["method"] == method and row["coverage"] == coverage
            ]
            method_summary = next(
                row
                for row in summaries
                if row["method"] == method and row["coverage"] == coverage
            )
            paired = mcnemar_exact(method_rows, bbs_rows)
            comparisons.append(
                {
                    "coverage": coverage,
                    "oligoark_method": method,
                    "external_method": "external_bbs",
                    "exact_recovery_rate_difference": round(
                        float(method_summary["exact_recovery_rate"])
                        - float(bbs_summary["exact_recovery_rate"]),
                        8,
                    ),
                    **paired,
                    "conclusion": _relation(method_summary, bbs_summary, paired),
                }
            )


    prior_comparisons: list[dict[str, Any]] = []
    for coverage in coverages:
        multistart_rows = [
            row
            for row in all_rows
            if row["method"] == "multistart_trace" and row["coverage"] == coverage
        ]
        multistart_summary = next(
            row
            for row in summaries
            if row["method"] == "multistart_trace" and row["coverage"] == coverage
        )
        for baseline_method in ("iterative_trace", "graph_alignment"):
            baseline_rows = [
                row
                for row in all_rows
                if row["method"] == baseline_method and row["coverage"] == coverage
            ]
            baseline_summary = next(
                row
                for row in summaries
                if row["method"] == baseline_method and row["coverage"] == coverage
            )
            paired = mcnemar_exact(multistart_rows, baseline_rows)
            lower_edit = 0
            equal_edit = 0
            higher_edit = 0
            baseline_by_cluster = {
                int(row["cluster_index"]): int(row["edit_distance"])
                for row in baseline_rows
            }
            for row in multistart_rows:
                cluster_index = int(row["cluster_index"])
                candidate_distance = int(row["edit_distance"])
                baseline_distance = baseline_by_cluster[cluster_index]
                if candidate_distance < baseline_distance:
                    lower_edit += 1
                elif candidate_distance == baseline_distance:
                    equal_edit += 1
                else:
                    higher_edit += 1
            prior_comparisons.append(
                {
                    "coverage": coverage,
                    "method": "multistart_trace",
                    "baseline": baseline_method,
                    "exact_recovery_rate_difference": round(
                        float(multistart_summary["exact_recovery_rate"])
                        - float(baseline_summary["exact_recovery_rate"]),
                        8,
                    ),
                    "mean_edit_distance_difference": round(
                        float(multistart_summary["mean_edit_distance"])
                        - float(baseline_summary["mean_edit_distance"]),
                        8,
                    ),
                    "lower_edit_distance_pairs": lower_edit,
                    "equal_edit_distance_pairs": equal_edit,
                    "higher_edit_distance_pairs": higher_edit,
                    **paired,
                }
            )

    empty_clusters = sum(not cluster for cluster in clusters)
    dataset_metadata = {
        "name": DATASET_NAME,
        "repository": DATASET_REPOSITORY,
        "commit": DATASET_COMMIT,
        "centers_git_blob": CENTERS_GIT_BLOB,
        "clusters_git_blob": CLUSTERS_GIT_BLOB,
        "centers_sha256": _sha256_file(centers_path),
        "clusters_sha256": _sha256_file(clusters_path),
        "reference_count": len(centers),
        "read_count": sum(len(cluster) for cluster in clusters),
        "empty_cluster_count": empty_clusters,
        "empty_cluster_rate": round(empty_clusters / len(clusters), 8),
        "sequencing_platform": "Oxford Nanopore Technologies MinION",
        "reference_length": TARGET_LENGTH,
        "reference_mapping": "cluster order maps explicitly to Centers.txt order",
    }
    metadata = {
        "oligoark_commit": _git_commit(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "dataset": dataset_metadata,
        "external_baseline": {
            "name": "Bidirectional Beam Search (BBS)",
            "repository": BBS_REPOSITORY,
            "commit": BBS_COMMIT,
            "configuration": "official defaults, one CPU thread, target length 110",
            "repeats_per_coverage": args.bbs_repeats,
            "reproducibility_note": (
                "The official implementation uses randomized Rust HashMap iteration and "
                "does not define a deterministic tie-break for equal-score candidates. "
                "Repeated official-code executions are therefore recorded without patching "
                "the external algorithm."
            ),
        },
        "subset": {
            "selection": (
                "deterministic SHA-256 rank among clusters with >= max tested coverage; "
                "selection uses cluster index only, never reconstruction quality"
            ),
            "seed": args.seed,
            "subset_size": args.subset_size,
            "coverages": list(coverages),
            "max_coverage": max_coverage,
            "eligible_cluster_count": sum(
                len(cluster) >= max_coverage for cluster in clusters
            ),
        },
        "calibration": {
            "selection": (
                "separate deterministic SHA-256 ranked clusters with all held-out "
                "cluster IDs excluded before parameter selection"
            ),
            "seed": args.calibration_seed,
            "subset_size": args.calibration_size,
            "coverages": list(CALIBRATION_COVERAGES),
            "held_out_overlap_count": 0,
            "baseline": calibration_baseline,
            "candidates": calibration_candidates,
            "winner": calibration_winner,
            "runtime_budget_rule": (
                "candidate runtime <= max(4x iterative baseline, baseline + 2 seconds)"
            ),
        },
        "claim_scope": (
            "physical-read trace reconstruction against explicit reference oligos; "
            "not OligoArk end-to-end archive decoding"
        ),
    }
    not_applicable = [
        {
            "method": "oligoark_fountain_hybrid",
            "reason": (
                "CNR strands were not encoded with OligoArk fountain/hybrid redundancy, so "
                "retroactively applying it would change the physical dataset."
            ),
        },
        {
            "method": "oligoark_combined_robust_optimizer",
            "reason": (
                "The optimizer chooses OligoArk archive/codec parameters before encoding; "
                "CNR is an externally encoded physical dataset, so this is not a fair "
                "same-read reconstruction comparison."
            ),
        },
        {
            "metric": "full_file_sha256_recovery",
            "reason": (
                "CNR publishes reference strands and clustered reads, not an original file "
                "encoded in the OligoArk archive format. Subset reference/reconstruction "
                "SHA-256 values are reported instead."
            ),
        },
    ]
    limitations = [
        "The fast CI profile uses a deterministic high-coverage subset, not all 10000 clusters.",
        (
            "Only read coverages 1, 5 and 10 are tested; these are not claims of true "
            "minimum coverage."
        ),
        (
            "The CNR maintainers note that its 10000 generated centers contain unintended "
            "long-range dependencies and some recovered clusters may be malformed."
        ),
        (
            "BBS is the single external recent baseline in this low-infrastructure profile. "
            "Codec-specific methods such as HEDGES cannot be fairly retrofitted onto these "
            "already synthesized CNR strands without changing the encoded data."
        ),
        (
            "The official BBS implementation can choose different equal-score candidates "
            "across processes because its Rust HashMap iteration order is randomized. "
            "The benchmark preserves the official source and records repeated executions "
            "instead of modifying its tie behavior."
        ),
        (
            "The multistart configuration is selected on one disjoint CNR calibration "
            "subset and may not generalize to other physical DNA-storage channels."
        ),
    ]

    error_profile = _error_profile(selected, max_coverage)
    summary_payload = {
        "metadata": metadata,
        "calibration": {
            "baseline": calibration_baseline,
            "candidates": calibration_candidates,
            "winner": calibration_winner,
            "held_out_overlap_count": 0,
        },
        "raw_read_error_profile_on_selected_max_coverage_reads": error_profile,
        "summaries": summaries,
        "paired_exact_recovery_comparisons_vs_bbs": comparisons,
        "paired_improvements_vs_prior_oligoark": prior_comparisons,
        "bbs_repeat_summaries": bbs_repeat_summaries,
        "not_applicable": not_applicable,
        "limitations": limitations,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary_payload, indent=2), encoding="utf-8"
    )
    (output_dir / "raw-results.json").write_text(
        json.dumps(all_rows, indent=2), encoding="utf-8"
    )
    _write_csv(output_dir / "raw-results.csv", all_rows)
    _write_csv(output_dir / "summary.csv", summaries)
    _write_csv(output_dir / "paired-comparisons.csv", comparisons)
    _write_csv(
        output_dir / "paired-improvements-vs-prior.csv",
        prior_comparisons,
    )
    _write_csv(output_dir / "bbs-repeat-summaries.csv", bbs_repeat_summaries)
    _write_csv(output_dir / "calibration-candidates.csv", calibration_candidates)
    _write_csv(
        output_dir / "calibration-summary.csv",
        [calibration_baseline, calibration_winner],
    )
    (output_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary_payload, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--centers")
    parser.add_argument("--clusters")
    parser.add_argument("--bbs-bin")
    parser.add_argument("--output-dir", default="external-cnr-results")
    parser.add_argument("--subset-size", type=int, default=DEFAULT_SUBSET_SIZE)
    parser.add_argument("--coverages", type=int, nargs="+", default=list(DEFAULT_COVERAGES))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--calibration-seed",
        type=int,
        default=DEFAULT_CALIBRATION_SEED,
    )
    parser.add_argument(
        "--calibration-size",
        type=int,
        default=DEFAULT_CALIBRATION_SIZE,
    )
    parser.add_argument("--bbs-repeats", type=int, default=5)
    parser.add_argument("--worker-method", choices=OLIGOARK_METHODS)
    parser.add_argument("--subset-json")
    parser.add_argument("--worker-output")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.worker_method:
        if not args.subset_json or not args.worker_output:
            raise ValueError("worker mode requires --subset-json and --worker-output")
        run_worker(args.worker_method, Path(args.subset_json), Path(args.worker_output))
        return
    run_benchmark(args)


if __name__ == "__main__":
    main()

"""Cross-dataset physical-read benchmark on the Grass et al. DNA-storage data."""

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
    confidence_fusion_trace_consensus,
    edit_distance,
    global_align,
)

DATASET_NAME = "Grass et al. prepared binned physical DNA-storage reads"
DATASET_RECORD = "https://zenodo.org/records/14296588"
DATASET_DOI = "10.5281/zenodo.14296588"
SOURCE_PAPER_DOI = "10.1002/anie.201411378"
DATASET_FILE = "Grass.txt"
DATASET_MD5 = "b076770ad26e955ff60b94e6344e3a05"
DATASET_LICENSE = "CC BY (prepared benchmark record)"
BBS_REPOSITORY = "GZHoffie/bbs"
BBS_COMMIT = "3e4ab46871929819e4f3e34a831c57cac88bb456"

DEVELOPMENT_SEED = 20261011
HELD_OUT_SEED = 20261012
DEVELOPMENT_SIZE = 48
HELD_OUT_SIZE = 96
COVERAGES = (1, 5, 10)

FROZEN_FUSION_CONFIG: dict[str, object] = {
    "name": "fusion-fast-t2-c8-r0-q025-g005",
    "anchors": 3,
    "rounds": 1,
    "top_positions": 2,
    "max_candidates": 8,
    "trim_farthest": 0,
    "qgram_width": 4,
    "qgram_weight": 0.25,
    "minimum_score_gain": 0.05,
}


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _validate_dna(value: str) -> str:
    sequence = value.strip().upper()
    if not sequence or any(base not in "ACGT" for base in sequence):
        raise ValueError("expected a non-empty DNA sequence over A/C/G/T")
    return sequence


def load_grass_binned(path: Path) -> list[dict[str, Any]]:
    """Parse the prepared Grass format: reference, star marker, reads, blank line(s)."""
    records: list[dict[str, Any]] = []
    reference: str | None = None
    reads: list[str] = []
    awaiting_marker = False

    def finish() -> None:
        nonlocal reference, reads, awaiting_marker
        if reference is not None:
            records.append(
                {
                    "cluster_index": len(records) + 1,
                    "reference": reference,
                    "reads": reads,
                }
            )
        reference = None
        reads = []
        awaiting_marker = False

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if reference is None:
                if not line:
                    continue
                reference = _validate_dna(line)
                awaiting_marker = True
                continue
            if awaiting_marker:
                if not line:
                    continue
                if set(line) != {"*"}:
                    raise ValueError(
                        f"expected star marker after cluster {len(records) + 1}"
                    )
                awaiting_marker = False
                continue
            if not line:
                finish()
                continue
            reads.append(_validate_dna(line))
    finish()

    if not records:
        raise ValueError("no Grass clusters parsed")
    lengths = {len(str(record["reference"])) for record in records}
    if len(lengths) != 1:
        raise ValueError(f"expected one reference length, found {sorted(lengths)}")
    return records


def _rank_key(seed: int, *parts: object) -> bytes:
    material = ":".join([str(seed), *(str(part) for part in parts)])
    return hashlib.sha256(material.encode("utf-8")).digest()


def select_split(
    records: list[dict[str, Any]],
    *,
    size: int,
    seed: int,
    max_coverage: int,
    excluded_indices: set[int] | None = None,
) -> list[dict[str, Any]]:
    excluded = excluded_indices or set()
    eligible = [
        record
        for record in records
        if len(record["reads"]) >= max_coverage
        and int(record["cluster_index"]) not in excluded
    ]
    if len(eligible) < size:
        raise ValueError(
            f"only {len(eligible)} eligible clusters; cannot select {size}"
        )
    selected = sorted(
        eligible,
        key=lambda record: _rank_key(seed, int(record["cluster_index"])),
    )[:size]
    output: list[dict[str, Any]] = []
    for record in selected:
        cluster_index = int(record["cluster_index"])
        ranked_reads = [
            read
            for _, read in sorted(
                enumerate(record["reads"]),
                key=lambda item: _rank_key(
                    seed,
                    cluster_index,
                    item[0],
                    item[1],
                ),
            )
        ]
        output.append(
            {
                "cluster_index": cluster_index,
                "reference": str(record["reference"]),
                "available_reads": len(record["reads"]),
                "ranked_reads": ranked_reads,
            }
        )
    return output


def _coverage_records(
    records: list[dict[str, Any]], coverage: int
) -> list[dict[str, Any]]:
    return [
        {
            "cluster_index": int(record["cluster_index"]),
            "reference": str(record["reference"]),
            "reads": list(record["ranked_reads"][:coverage]),
        }
        for record in records
    ]


def _classify_error(reference: str, reconstruction: str) -> dict[str, Any]:
    if reference == reconstruction:
        return {
            "error_class": "exact",
            "substitutions": 0,
            "insertions": 0,
            "deletions": 0,
        }
    aligned_reference, aligned_reconstruction = global_align(reference, reconstruction)
    substitutions = 0
    insertions = 0
    deletions = 0
    for ref_base, out_base in zip(
        aligned_reference, aligned_reconstruction, strict=True
    ):
        if ref_base == "-":
            insertions += 1
        elif out_base == "-":
            deletions += 1
        elif ref_base != out_base:
            substitutions += 1
    distance = edit_distance(reference, reconstruction)
    if distance == 1 and substitutions == 1 and not insertions and not deletions:
        error_class = "one_edit_substitution"
    elif distance == 1 and insertions == 1 and not substitutions and not deletions:
        error_class = "one_edit_insertion"
    elif distance == 1 and deletions == 1 and not substitutions and not insertions:
        error_class = "one_edit_deletion"
    elif insertions and deletions and len(reference) == len(reconstruction):
        error_class = "alignment_shift"
    elif len(reference) != len(reconstruction):
        error_class = "length_error"
    else:
        error_class = "multi_edit"
    return {
        "error_class": error_class,
        "substitutions": substitutions,
        "insertions": insertions,
        "deletions": deletions,
    }


def _case_record(
    *,
    method: str,
    family: str,
    coverage: int,
    cluster_index: int,
    reference: str,
    reconstruction: str,
) -> dict[str, Any]:
    distance = edit_distance(reference, reconstruction)
    return {
        "method": method,
        "method_family": family,
        "coverage": coverage,
        "cluster_index": cluster_index,
        "exact": reconstruction == reference,
        "edit_distance": distance,
        "normalized_edit_distance": round(distance / len(reference), 8),
        "reference": reference,
        "reconstruction": reconstruction,
        "reconstructed_length": len(reconstruction),
        **_classify_error(reference, reconstruction),
    }


def reconstruct_confidence_fusion(reads: list[str], target_length: int) -> str:
    config = FROZEN_FUSION_CONFIG
    return confidence_fusion_trace_consensus(
        reads,
        target_length=target_length,
        anchors=int(config["anchors"]),
        rounds=int(config["rounds"]),
        top_positions=int(config["top_positions"]),
        max_candidates=int(config["max_candidates"]),
        trim_farthest=int(config["trim_farthest"]),
        qgram_width=int(config["qgram_width"]),
        qgram_weight=float(config["qgram_weight"]),
        minimum_score_gain=float(config["minimum_score_gain"]),
    )


def run_worker(subset_json: Path, worker_output: Path) -> None:
    payload = json.loads(subset_json.read_text(encoding="utf-8"))
    coverage = int(payload["coverage"])
    target_length = int(payload["target_length"])
    rows: list[dict[str, Any]] = []
    for record in payload["records"]:
        reads = [str(read) for read in record["reads"]]
        reference = str(record["reference"])
        reconstruction = reconstruct_confidence_fusion(reads, target_length)
        rows.append(
            _case_record(
                method="confidence_fusion",
                family="oligoark",
                coverage=coverage,
                cluster_index=int(record["cluster_index"]),
                reference=reference,
                reconstruction=reconstruction,
            )
        )
    worker_output.write_text(json.dumps(rows, indent=2), encoding="utf-8")


def write_bbs_clusters(path: Path, records: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write("========================================\n")
            for read in record["reads"]:
                handle.write(f"{read}\n")


def _timed_subprocess(command: list[str], timing_path: Path) -> tuple[float, float | None]:
    gnu_time = Path("/usr/bin/time")
    if gnu_time.exists():
        subprocess.run(
            [
                str(gnu_time),
                "-f",
                "%e,%M",
                "-o",
                str(timing_path),
                *command,
            ],
            check=True,
        )
        elapsed_text, rss_text = (
            timing_path.read_text(encoding="utf-8").strip().split(",")
        )
        return float(elapsed_text), round(float(rss_text) / 1024.0, 4)
    started = time.perf_counter()
    subprocess.run(command, check=True)
    return round(time.perf_counter() - started, 6), None


def _parse_bbs(
    path: Path,
    records: list[dict[str, Any]],
    coverage: int,
) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        parsed = list(csv.DictReader(handle))
    if len(parsed) != len(records):
        raise ValueError(
            f"BBS returned {len(parsed)} rows for {len(records)} input clusters"
        )
    rows: list[dict[str, Any]] = []
    for output_row, record in zip(parsed, records, strict=True):
        reference = str(record["reference"])
        reconstruction = str(output_row["reconstruction_result"]).strip().upper()
        row = _case_record(
            method="external_bbs",
            family="external",
            coverage=coverage,
            cluster_index=int(record["cluster_index"]),
            reference=reference,
            reconstruction=reconstruction,
        )
        row.update(
            {
                "bbs_k": int(output_row["k"]),
                "bbs_path_weight": float(output_row["path_weight"]),
                "bbs_confidence": float(output_row["confidence"]),
            }
        )
        rows.append(row)
    return rows


def wilson_interval(successes: int, trials: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = successes / trials
    denominator = 1.0 + z * z / trials
    center = (p + z * z / (2 * trials)) / denominator
    spread = (
        z
        * math.sqrt(p * (1.0 - p) / trials + z * z / (4 * trials * trials))
        / denominator
    )
    return max(0.0, center - spread), min(1.0, center + spread)


def summarize_rows(
    rows: list[dict[str, Any]],
    elapsed_seconds: float,
    peak_rss_mb: float | None,
) -> dict[str, Any]:
    successes = sum(bool(row["exact"]) for row in rows)
    trials = len(rows)
    low, high = wilson_interval(successes, trials)
    distances = [int(row["edit_distance"]) for row in rows]
    references = "".join(str(row["reference"]) for row in rows)
    reconstructions = "".join(str(row["reconstruction"]) for row in rows)
    failure_classes: dict[str, int] = {}
    for row in rows:
        label = str(row["error_class"])
        failure_classes[label] = failure_classes.get(label, 0) + 1
    return {
        "method": str(rows[0]["method"]),
        "method_family": str(rows[0]["method_family"]),
        "coverage": int(rows[0]["coverage"]),
        "successes": successes,
        "trials": trials,
        "exact_recovery_rate": round(successes / trials, 8),
        "recovery_ci95_low": round(low, 8),
        "recovery_ci95_high": round(high, 8),
        "mean_edit_distance": round(statistics.fmean(distances), 8),
        "median_edit_distance": round(float(statistics.median(distances)), 8),
        "one_edit_failures": sum(distance == 1 for distance in distances),
        "two_edit_failures": sum(distance == 2 for distance in distances),
        "three_plus_edit_failures": sum(distance >= 3 for distance in distances),
        "max_edit_distance": max(distances, default=0),
        "failure_classes": dict(sorted(failure_classes.items())),
        "elapsed_seconds": round(elapsed_seconds, 6),
        "seconds_per_cluster": round(elapsed_seconds / trials, 8),
        "peak_rss_mb": peak_rss_mb,
        "reference_subset_sha256": _sha256_text(references),
        "reconstructed_subset_sha256": _sha256_text(reconstructions),
    }


def mcnemar_exact(
    left: list[dict[str, Any]], right: list[dict[str, Any]]
) -> dict[str, Any]:
    left_by_id = {int(row["cluster_index"]): bool(row["exact"]) for row in left}
    right_by_id = {int(row["cluster_index"]): bool(row["exact"]) for row in right}
    if left_by_id.keys() != right_by_id.keys():
        raise ValueError("paired methods must contain identical cluster IDs")
    left_only = sum(
        left_by_id[index] and not right_by_id[index] for index in left_by_id
    )
    right_only = sum(
        right_by_id[index] and not left_by_id[index] for index in left_by_id
    )
    discordant = left_only + right_only
    if discordant == 0:
        p_value = 1.0
    else:
        tail = sum(
            math.comb(discordant, value)
            for value in range(min(left_only, right_only) + 1)
        ) / (2**discordant)
        p_value = min(1.0, 2.0 * tail)
    return {
        "oligoark_only_successes": left_only,
        "bbs_only_successes": right_only,
        "discordant_pairs": discordant,
        "two_sided_exact_p": p_value,
    }


def _read_error_profile(
    records: list[dict[str, Any]],
    max_coverage: int,
) -> dict[str, Any]:
    substitutions = 0
    insertions = 0
    deletions = 0
    reference_bases = 0
    reads_profiled = 0
    for record in records:
        reference = str(record["reference"])
        for read in record["ranked_reads"][:max_coverage]:
            aligned_reference, aligned_read = global_align(reference, str(read))
            for ref_base, read_base in zip(
                aligned_reference, aligned_read, strict=True
            ):
                if ref_base == "-":
                    insertions += 1
                elif read_base == "-":
                    deletions += 1
                    reference_bases += 1
                else:
                    substitutions += int(ref_base != read_base)
                    reference_bases += 1
            reads_profiled += 1
    denominator = max(1, reference_bases)
    return {
        "reads_profiled": reads_profiled,
        "reference_bases_profiled": reference_bases,
        "substitutions": substitutions,
        "insertions": insertions,
        "deletions": deletions,
        "substitution_rate_per_reference_base": round(
            substitutions / denominator, 8
        ),
        "insertion_rate_per_reference_base": round(insertions / denominator, 8),
        "deletion_rate_per_reference_base": round(deletions / denominator, 8),
    }


def _split_metadata(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "cluster_index": int(record["cluster_index"]),
            "available_reads": int(record["available_reads"]),
            "reference_sha256": _sha256_text(str(record["reference"])),
        }
        for record in records
    ]


def run_benchmark(args: argparse.Namespace) -> None:
    dataset_path = Path(args.dataset)
    bbs_bin = Path(args.bbs_bin)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    timing_dir = output_dir / "timing"
    timing_dir.mkdir(exist_ok=True)

    if _md5_file(dataset_path) != DATASET_MD5:
        raise ValueError("Grass.txt MD5 does not match the pinned Zenodo record")

    records = load_grass_binned(dataset_path)
    target_length = len(str(records[0]["reference"]))
    max_coverage = max(COVERAGES)
    eligible_count = sum(len(record["reads"]) >= max_coverage for record in records)

    development = select_split(
        records,
        size=args.development_size,
        seed=args.development_seed,
        max_coverage=max_coverage,
    )
    development_ids = {int(record["cluster_index"]) for record in development}
    held_out = select_split(
        records,
        size=args.held_out_size,
        seed=args.held_out_seed,
        max_coverage=max_coverage,
        excluded_indices=development_ids,
    )
    held_out_ids = {int(record["cluster_index"]) for record in held_out}
    if development_ids & held_out_ids:
        raise RuntimeError("development and held-out Grass splits overlap")

    (output_dir / "development-clusters.json").write_text(
        json.dumps(_split_metadata(development), indent=2),
        encoding="utf-8",
    )
    (output_dir / "held-out-clusters.json").write_text(
        json.dumps(_split_metadata(held_out), indent=2),
        encoding="utf-8",
    )

    development_summaries: list[dict[str, Any]] = []
    for coverage in (5, 10):
        started = time.perf_counter()
        rows: list[dict[str, Any]] = []
        for record in _coverage_records(development, coverage):
            reference = str(record["reference"])
            reconstruction = reconstruct_confidence_fusion(
                [str(read) for read in record["reads"]],
                target_length,
            )
            rows.append(
                _case_record(
                    method="confidence_fusion",
                    family="oligoark",
                    coverage=coverage,
                    cluster_index=int(record["cluster_index"]),
                    reference=reference,
                    reconstruction=reconstruction,
                )
            )
        development_summaries.append(
            summarize_rows(rows, time.perf_counter() - started, None)
        )

    summaries: list[dict[str, Any]] = []
    raw_rows: list[dict[str, Any]] = []
    bbs_repeat_summaries: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []

    for coverage in COVERAGES:
        coverage_records = _coverage_records(held_out, coverage)
        subset_json = output_dir / f"subset-coverage-{coverage}.json"
        subset_json.write_text(
            json.dumps(
                {
                    "coverage": coverage,
                    "target_length": target_length,
                    "records": coverage_records,
                }
            ),
            encoding="utf-8",
        )

        worker_output = output_dir / f"confidence-fusion-{coverage}.json"
        elapsed, peak_rss = _timed_subprocess(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--worker",
                "--subset-json",
                str(subset_json),
                "--worker-output",
                str(worker_output),
            ],
            timing_dir / f"confidence-fusion-{coverage}.txt",
        )
        oligoark_rows = json.loads(worker_output.read_text(encoding="utf-8"))
        raw_rows.extend(oligoark_rows)
        summaries.append(summarize_rows(oligoark_rows, elapsed, peak_rss))

        bbs_input = output_dir / f"bbs-input-{coverage}.txt"
        write_bbs_clusters(bbs_input, coverage_records)
        first_bbs_rows: list[dict[str, Any]] | None = None
        for repeat in range(1, args.bbs_repeats + 1):
            bbs_output = output_dir / f"bbs-{coverage}-repeat-{repeat}.csv"
            bbs_elapsed, bbs_rss = _timed_subprocess(
                [
                    str(bbs_bin),
                    str(bbs_input),
                    "-l",
                    str(target_length),
                    "-t",
                    "1",
                    "-o",
                    str(bbs_output),
                ],
                timing_dir / f"bbs-{coverage}-repeat-{repeat}.txt",
            )
            bbs_rows = _parse_bbs(bbs_output, coverage_records, coverage)
            bbs_summary = summarize_rows(bbs_rows, bbs_elapsed, bbs_rss)
            bbs_summary["repeat"] = repeat
            bbs_repeat_summaries.append(bbs_summary)
            if repeat == 1:
                first_bbs_rows = bbs_rows
                raw_rows.extend(bbs_rows)
                summaries.append(bbs_summary)

        assert first_bbs_rows is not None
        comparisons.append(
            {
                "coverage": coverage,
                **mcnemar_exact(oligoark_rows, first_bbs_rows),
            }
        )

    metadata = {
        "dataset": {
            "name": DATASET_NAME,
            "record": DATASET_RECORD,
            "doi": DATASET_DOI,
            "source_paper_doi": SOURCE_PAPER_DOI,
            "file": DATASET_FILE,
            "file_md5": _md5_file(dataset_path),
            "file_sha256": _sha256_file(dataset_path),
            "license": DATASET_LICENSE,
            "cluster_count": len(records),
            "eligible_clusters_with_at_least_10_reads": eligible_count,
            "reference_length": target_length,
            "synthesis": "CustomArray electrochemical oligo synthesis",
            "sequencing": "Illumina MiSeq 2x150 TruSeq",
            "prepared_format": (
                "each cluster stores the encoded reference sequence followed by its "
                "associated physical reads"
            ),
        },
        "oligoark": {
            "algorithm": "confidence_fusion_trace_consensus",
            "frozen_config": FROZEN_FUSION_CONFIG,
            "dataset_specific_change": (
                "target_length inferred from Grass references; no scoring/search "
                "parameters changed from the merged CNR configuration"
            ),
        },
        "external_baseline": {
            "name": "Bidirectional Beam Search (BBS)",
            "repository": BBS_REPOSITORY,
            "commit": BBS_COMMIT,
            "threads": 1,
            "repeats_per_coverage": args.bbs_repeats,
        },
        "splits": {
            "development_seed": args.development_seed,
            "development_size": args.development_size,
            "held_out_seed": args.held_out_seed,
            "held_out_size": args.held_out_size,
            "overlap_count": 0,
            "read_ranking": "deterministic SHA-256 rank within each selected cluster",
        },
        "python": platform.python_version(),
        "platform": platform.platform(),
    }
    payload = {
        "metadata": metadata,
        "development_diagnostic_only": development_summaries,
        "held_out_summaries": summaries,
        "paired_exact_recovery_vs_first_bbs_repeat": comparisons,
        "bbs_repeat_summaries": bbs_repeat_summaries,
        "held_out_raw_read_error_profile": _read_error_profile(
            held_out, max(COVERAGES)
        ),
        "claim_boundary": (
            "reference-level physical-read reconstruction on an independently "
            "published Illumina DNA-storage dataset; not end-to-end archive decoding"
        ),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
    (output_dir / "raw-results.json").write_text(
        json.dumps(raw_rows, indent=2),
        encoding="utf-8",
    )
    if raw_rows:
        fieldnames = sorted(
            {key for row in raw_rows for key in row},
            key=lambda key: (
                key
                not in {
                    "method",
                    "method_family",
                    "coverage",
                    "cluster_index",
                    "exact",
                    "edit_distance",
                },
                key,
            ),
        )
        with (output_dir / "raw-results.csv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(raw_rows)
    print(json.dumps(payload, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset")
    parser.add_argument("--bbs-bin")
    parser.add_argument("--output-dir", default="external-grass-results")
    parser.add_argument("--development-seed", type=int, default=DEVELOPMENT_SEED)
    parser.add_argument("--held-out-seed", type=int, default=HELD_OUT_SEED)
    parser.add_argument("--development-size", type=int, default=DEVELOPMENT_SIZE)
    parser.add_argument("--held-out-size", type=int, default=HELD_OUT_SIZE)
    parser.add_argument("--bbs-repeats", type=int, default=5)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--subset-json")
    parser.add_argument("--worker-output")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.worker:
        if not args.subset_json or not args.worker_output:
            raise ValueError("worker mode requires --subset-json and --worker-output")
        run_worker(Path(args.subset_json), Path(args.worker_output))
        return
    if not args.dataset or not args.bbs_bin:
        raise ValueError("--dataset and --bbs-bin are required")
    if args.bbs_repeats < 1:
        raise ValueError("--bbs-repeats must be positive")
    run_benchmark(args)


if __name__ == "__main__":
    main()

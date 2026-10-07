"""Third independent physical-read benchmark on the LCRC HFS-Pool-11.7K dataset."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
import statistics
import subprocess
import sys
import time
from collections import Counter
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from oligoark.reconstruct import (
    confidence_fusion_trace_consensus,
    edit_distance,
    global_align,
)

DATASET_NAME = "LCRC HFS-Pool-11.7K + NGS"
DATASET_REPOSITORY = "dna-storage-lab/DNAStorage_LCRC"
DATASET_COMMIT = "ce35bb2720c49ea6d1e6bc17903655b45f9c6c51"
REFERENCE_PATH = "LCRC_SmallScale/reference/DNA_oligoPool/oligoPool_11.7K.fa"
READS_PATH = "LCRC_SmallScale/fastq/HFS_Pool_11.7K_NGS_PEmerged.fastq.gz"
REFERENCE_GIT_BLOB = "3676efd8535e0b5554408154234ed4e2dd1d18ab"
READS_GIT_BLOB = "c9926ac7e348e0aba8906ccf6dfd02133ad6b433"
SOURCE_PAPER_DOI = "10.1126/sciadv.aec1469"
SRA_ACCESSION = "PRJNA1371011"
DATASET_LICENSE = "MIT software repository; sequencing data cited to SRA PRJNA1371011"
BBS_REPOSITORY = "GZHoffie/bbs"
BBS_COMMIT = "3e4ab46871929819e4f3e34a831c57cac88bb456"

REFERENCE_COUNT = 11_745
TARGET_LENGTH = 200
COVERAGES = (1, 5, 10)
DEVELOPMENT_SEED = 20261013
HELD_OUT_SEED = 20261014
DEVELOPMENT_SIZE = 48
HELD_OUT_SIZE = 96

ANCHOR_K = 13
ANCHOR_STEP = 8
ANCHOR_LEFT = 20
ANCHOR_RIGHT = 20
MIN_ANCHOR_HITS = 3
MIN_ANCHOR_MARGIN = 1
MAX_ASSIGNMENT_DISTANCE = 15

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

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
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


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def load_fasta(path: Path) -> list[str]:
    sequences: list[str] = []
    current: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences.append(_validate_dna("".join(current)))
                    current = []
                continue
            current.append(line)
    if current:
        sequences.append(_validate_dna("".join(current)))
    if not sequences:
        raise ValueError("no FASTA sequences parsed")
    return sequences


def iter_fastq_gz(path: Path) -> Iterable[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip().upper()
            plus = handle.readline()
            quality = handle.readline()
            if not plus or not quality:
                raise ValueError("truncated FASTQ record")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError("malformed FASTQ record")
            if sequence and all(base in "ACGT" for base in sequence):
                yield sequence


def _rank_key(seed: int, *parts: object) -> bytes:
    material = ":".join([str(seed), *(str(part) for part in parts)])
    return hashlib.sha256(material.encode("utf-8")).digest()


def build_unique_anchor_index(references: list[str]) -> dict[str, int]:
    owners: dict[str, int] = {}
    ambiguous: set[str] = set()
    for index, reference in enumerate(references):
        stop = len(reference) - ANCHOR_RIGHT - ANCHOR_K + 1
        for position in range(ANCHOR_LEFT, stop, ANCHOR_STEP):
            kmer = reference[position : position + ANCHOR_K]
            previous = owners.get(kmer)
            if previous is None and kmer not in ambiguous:
                owners[kmer] = index
            elif previous != index:
                owners.pop(kmer, None)
                ambiguous.add(kmer)
    return owners


def _bounded_edit_distance(left: str, right: str, limit: int) -> int:
    if abs(len(left) - len(right)) > limit:
        return limit + 1
    if len(left) == len(right):
        mismatches = sum(a != b for a, b in zip(left, right, strict=True))
        if mismatches <= limit:
            return mismatches

    previous = {column: column for column in range(0, min(len(right), limit) + 1)}
    for row, left_base in enumerate(left, start=1):
        start = max(0, row - limit)
        end = min(len(right), row + limit)
        current: dict[int, int] = {}
        for column in range(start, end + 1):
            if column == 0:
                current[column] = row
                continue
            best = limit + 1
            if column in previous:
                best = min(best, previous[column] + 1)
            if column - 1 in current:
                best = min(best, current[column - 1] + 1)
            if column - 1 in previous:
                best = min(
                    best,
                    previous[column - 1]
                    + int(left_base != right[column - 1]),
                )
            current[column] = best
        if current and min(current.values()) > limit:
            return limit + 1
        previous = current
    return previous.get(len(right), limit + 1)


def assign_read(
    read: str,
    references: list[str],
    anchor_index: dict[str, int],
) -> tuple[int, str, int, int] | None:
    best_assignment: tuple[int, str, int, int, int] | None = None
    for oriented in (read, reverse_complement(read)):
        counts: Counter[int] = Counter()
        seen_kmers: set[str] = set()
        for position in range(0, len(oriented) - ANCHOR_K + 1):
            kmer = oriented[position : position + ANCHOR_K]
            if kmer in seen_kmers:
                continue
            seen_kmers.add(kmer)
            owner = anchor_index.get(kmer)
            if owner is not None:
                counts[owner] += 1
        if not counts:
            continue
        ranked = counts.most_common(2)
        candidate_index, hits = ranked[0]
        second_hits = ranked[1][1] if len(ranked) > 1 else 0
        margin = hits - second_hits
        if hits < MIN_ANCHOR_HITS or margin < MIN_ANCHOR_MARGIN:
            continue
        distance = _bounded_edit_distance(
            references[candidate_index],
            oriented,
            MAX_ASSIGNMENT_DISTANCE,
        )
        if distance > MAX_ASSIGNMENT_DISTANCE:
            continue
        assignment = (distance, -hits, -margin, candidate_index, oriented)
        if best_assignment is None or assignment < best_assignment:
            best_assignment = assignment

    if best_assignment is None:
        return None
    distance, negative_hits, negative_margin, candidate_index, oriented = best_assignment
    return candidate_index, oriented, -negative_hits, distance


def cluster_reads(
    references: list[str],
    reads_path: Path,
) -> tuple[list[list[str]], dict[str, Any]]:
    anchor_index = build_unique_anchor_index(references)
    clusters: list[list[str]] = [[] for _ in references]
    total = 0
    mapped = 0
    rejected = 0
    read_lengths: Counter[int] = Counter()
    anchor_hits: list[int] = []
    assignment_distances: list[int] = []

    for read in iter_fastq_gz(reads_path):
        total += 1
        read_lengths[len(read)] += 1
        assigned = assign_read(read, references, anchor_index)
        if assigned is None:
            rejected += 1
            continue
        index, oriented, hits, distance = assigned
        clusters[index].append(oriented)
        mapped += 1
        anchor_hits.append(hits)
        assignment_distances.append(distance)

    populated = sum(bool(cluster) for cluster in clusters)
    eligible_10 = sum(len(cluster) >= 10 for cluster in clusters)
    return clusters, {
        "total_fastq_reads": total,
        "mapped_reads": mapped,
        "rejected_reads": rejected,
        "mapping_rate": round(mapped / max(1, total), 8),
        "unique_anchor_kmers": len(anchor_index),
        "populated_reference_clusters": populated,
        "clusters_with_at_least_10_reads": eligible_10,
        "read_length_counts": dict(sorted(read_lengths.items())),
        "median_anchor_hits": (
            float(statistics.median(anchor_hits)) if anchor_hits else 0.0
        ),
        "median_assignment_distance": (
            float(statistics.median(assignment_distances))
            if assignment_distances
            else 0.0
        ),
        "assignment_rule": {
            "anchor_k": ANCHOR_K,
            "anchor_step": ANCHOR_STEP,
            "anchor_left_exclusion": ANCHOR_LEFT,
            "anchor_right_exclusion": ANCHOR_RIGHT,
            "minimum_anchor_hits": MIN_ANCHOR_HITS,
            "minimum_anchor_margin": MIN_ANCHOR_MARGIN,
            "maximum_reference_assignment_edit_distance": MAX_ASSIGNMENT_DISTANCE,
            "orientation": "evaluate forward and reverse complement deterministically",
        },
    }


def select_split(
    references: list[str],
    clusters: list[list[str]],
    *,
    size: int,
    seed: int,
    max_coverage: int,
    excluded_indices: set[int] | None = None,
) -> list[dict[str, Any]]:
    excluded = excluded_indices or set()
    eligible = [
        index
        for index, cluster in enumerate(clusters)
        if len(cluster) >= max_coverage and index not in excluded
    ]
    if len(eligible) < size:
        raise ValueError(
            f"only {len(eligible)} clusters have at least {max_coverage} reads; "
            f"cannot select {size}"
        )
    selected = sorted(eligible, key=lambda index: _rank_key(seed, index))[:size]
    records: list[dict[str, Any]] = []
    for index in selected:
        ranked_reads = [
            read
            for _, read in sorted(
                enumerate(clusters[index]),
                key=lambda item: _rank_key(seed, index, item[0], item[1]),
            )
        ]
        records.append(
            {
                "cluster_index": index + 1,
                "reference": references[index],
                "available_reads": len(clusters[index]),
                "ranked_reads": ranked_reads,
            }
        )
    return records


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
        label = "one_edit_substitution"
    elif distance == 1 and insertions == 1 and not substitutions and not deletions:
        label = "one_edit_insertion"
    elif distance == 1 and deletions == 1 and not substitutions and not insertions:
        label = "one_edit_deletion"
    elif insertions and deletions and len(reference) == len(reconstruction):
        label = "alignment_shift"
    elif len(reference) != len(reconstruction):
        label = "length_error"
    else:
        label = "multi_edit"
    return {
        "error_class": label,
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
    classes: Counter[str] = Counter(str(row["error_class"]) for row in rows)
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
        "failure_classes": dict(sorted(classes.items())),
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
    records: list[dict[str, Any]], max_coverage: int
) -> dict[str, Any]:
    substitutions = 0
    insertions = 0
    deletions = 0
    bases = 0
    reads = 0
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
                    bases += 1
                else:
                    substitutions += int(ref_base != read_base)
                    bases += 1
            reads += 1
    denominator = max(1, bases)
    return {
        "reads_profiled": reads,
        "reference_bases_profiled": bases,
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
    reference_path = Path(args.references)
    reads_path = Path(args.reads)
    bbs_bin = Path(args.bbs_bin)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    timing_dir = output_dir / "timing"
    timing_dir.mkdir(exist_ok=True)

    references = load_fasta(reference_path)
    if len(references) != REFERENCE_COUNT:
        raise ValueError(
            f"expected {REFERENCE_COUNT} references, got {len(references)}"
        )
    lengths = {len(reference) for reference in references}
    if lengths != {TARGET_LENGTH}:
        raise ValueError(f"expected only {TARGET_LENGTH}-nt references, got {lengths}")

    mapping_started = time.perf_counter()
    clusters, mapping = cluster_reads(references, reads_path)
    mapping_seconds = time.perf_counter() - mapping_started

    development = select_split(
        references,
        clusters,
        size=args.development_size,
        seed=args.development_seed,
        max_coverage=max(COVERAGES),
    )
    development_ids = {int(row["cluster_index"]) - 1 for row in development}
    held_out = select_split(
        references,
        clusters,
        size=args.held_out_size,
        seed=args.held_out_seed,
        max_coverage=max(COVERAGES),
        excluded_indices=development_ids,
    )
    held_out_ids = {int(row["cluster_index"]) - 1 for row in held_out}
    if development_ids & held_out_ids:
        raise RuntimeError("development and held-out LCRC splits overlap")

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
                TARGET_LENGTH,
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
    comparisons: list[dict[str, Any]] = []
    bbs_repeat_summaries: list[dict[str, Any]] = []

    for coverage in COVERAGES:
        coverage_records = _coverage_records(held_out, coverage)
        subset_json = output_dir / f"subset-coverage-{coverage}.json"
        subset_json.write_text(
            json.dumps(
                {
                    "coverage": coverage,
                    "target_length": TARGET_LENGTH,
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
                    str(TARGET_LENGTH),
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
            "repository": DATASET_REPOSITORY,
            "commit": DATASET_COMMIT,
            "source_paper_doi": SOURCE_PAPER_DOI,
            "sra_accession": SRA_ACCESSION,
            "license_note": DATASET_LICENSE,
            "reference_path": REFERENCE_PATH,
            "reads_path": READS_PATH,
            "reference_git_blob": REFERENCE_GIT_BLOB,
            "reads_git_blob": READS_GIT_BLOB,
            "reference_sha256": _sha256_file(reference_path),
            "reads_sha256": _sha256_file(reads_path),
            "reference_count": len(references),
            "reference_length": TARGET_LENGTH,
            "synthesis": "Twist Bioscience high-fidelity HFS-Pool-11.7K",
            "sequencing": "Illumina/NGS PE150, repository-provided paired-end merged FASTQ",
        },
        "mapping": {
            **mapping,
            "elapsed_seconds": round(mapping_seconds, 6),
            "boundary": (
                "reference library used only for deterministic read-to-strand "
                "association; reconstruction receives observed clustered reads only"
            ),
        },
        "oligoark": {
            "algorithm": "confidence_fusion_trace_consensus",
            "frozen_config": FROZEN_FUSION_CONFIG,
            "dataset_specific_change": (
                "target length set mechanically to the published 200-nt design; "
                "all scoring/search parameters unchanged from CNR and Grass"
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
            "read_ranking": "deterministic SHA-256 rank within mapped cluster",
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
            "reference-level reconstruction on a third independent physical DNA-storage "
            "experiment; not end-to-end decoding through the authors' LCRC/ECC codec"
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
        fieldnames = sorted({key for row in raw_rows for key in row})
        with (output_dir / "raw-results.csv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(raw_rows)

    print(json.dumps(payload, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--references")
    parser.add_argument("--reads")
    parser.add_argument("--bbs-bin")
    parser.add_argument("--output-dir", default="external-lcrc-results")
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
    if not args.references or not args.reads or not args.bbs_bin:
        raise ValueError("--references, --reads and --bbs-bin are required")
    if args.bbs_repeats < 1:
        raise ValueError("--bbs-repeats must be positive")
    run_benchmark(args)


if __name__ == "__main__":
    main()

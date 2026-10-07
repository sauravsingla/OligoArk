"""Third independent physical-read benchmark on the LCRC HFS 11.7K ONT pool."""

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
import tarfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, BinaryIO, Iterable

from oligoark.reconstruct import confidence_fusion_trace_consensus, edit_distance, global_align

DATASET_NAME = "LCRC HFS-Pool-11.7K ONT FAST real-time physical reads"
DATASET_REPOSITORY = "dna-storage-lab/DNAStorage_LCRC"
DATASET_COMMIT = "ce35bb2720c49ea6d1e6bc17903655b45f9c6c51"
SOURCE_PAPER_DOI = "10.1126/sciadv.aec1469"
SRA_ACCESSION = "PRJNA1371011"
REFERENCE_PATH = "LCRC_SmallScale/reference/DNA_oligoPool/oligoPool_11.7K.fa"
READ_ARCHIVE_PATH = "LCRC_SmallScale/fastq/HFS_Pool_11.7K_ONT_FAST_RT.tar.gz"
REFERENCE_GIT_BLOB = "3676efd8535e0b5554408154234ed4e2dd1d18ab"
READ_ARCHIVE_GIT_BLOB = "e1529d774ef1b327d0f601c4890c80db9cde90e6"
REPOSITORY_LICENSE = "MIT"
BBS_REPOSITORY = "GZHoffie/bbs"
BBS_COMMIT = "3e4ab46871929819e4f3e34a831c57cac88bb456"

REFERENCE_LENGTH = 200
FORWARD_PRIMER = "ATAATTGGCTCCTGCTTGCA"
REVERSE_PRIMER = "AATGTAGGCGGAAAGTGCAA"
PRIMER_ANCHOR_LENGTH = 12
PRIMER_MAX_MISMATCHES = 2
PRIMER_SEARCH_SPAN = 70
MAPPING_K = 17
MAPPING_MIN_VOTES = 4
MAPPING_MIN_MARGIN = 2

DEVELOPMENT_SEED = 20261015
HELD_OUT_SEED = 20261016
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


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _validate_dna(value: str, *, expected_length: int | None = None) -> str:
    sequence = value.strip().upper()
    if not sequence or any(base not in "ACGT" for base in sequence):
        raise ValueError("expected a non-empty DNA sequence over A/C/G/T")
    if expected_length is not None and len(sequence) != expected_length:
        raise ValueError(
            f"expected sequence length {expected_length}, got {len(sequence)}"
        )
    return sequence


def load_references(path: Path) -> list[str]:
    references: list[str] = []
    current: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    references.append(
                        _validate_dna("".join(current), expected_length=REFERENCE_LENGTH)
                    )
                    current = []
                continue
            current.append(line)
    if current:
        references.append(
            _validate_dna("".join(current), expected_length=REFERENCE_LENGTH)
        )
    if len(references) != 11_745:
        raise ValueError(f"expected 11745 LCRC references, got {len(references)}")
    if not all(
        reference.startswith(FORWARD_PRIMER) and reference.endswith(REVERSE_PRIMER)
        for reference in references
    ):
        raise ValueError("LCRC references do not match the pinned universal primers")
    return references


_BASE_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}


def _encode_kmer(sequence: str) -> int | None:
    value = 0
    for base in sequence:
        bits = _BASE_BITS.get(base)
        if bits is None:
            return None
        value = (value << 2) | bits
    return value


def build_unique_kmer_index(
    references: list[str],
    *,
    k: int = MAPPING_K,
) -> dict[int, int]:
    """Map exact payload k-mers to one unique reference index, or -1 if ambiguous."""
    index: dict[int, int] = {}
    for reference_index, reference in enumerate(references):
        payload = reference[len(FORWARD_PRIMER) : -len(REVERSE_PRIMER)]
        seen_in_reference: set[int] = set()
        for position in range(len(payload) - k + 1):
            encoded = _encode_kmer(payload[position : position + k])
            if encoded is None or encoded in seen_in_reference:
                continue
            seen_in_reference.add(encoded)
            previous = index.get(encoded)
            if previous is None:
                index[encoded] = reference_index
            elif previous != reference_index:
                index[encoded] = -1
    return index


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _vote_reference(
    sequence: str,
    kmer_index: dict[int, int],
    *,
    k: int = MAPPING_K,
) -> tuple[int | None, int, int]:
    votes: Counter[int] = Counter()
    seen_kmers: set[int] = set()
    for position in range(len(sequence) - k + 1):
        encoded = _encode_kmer(sequence[position : position + k])
        if encoded is None or encoded in seen_kmers:
            continue
        seen_kmers.add(encoded)
        reference_index = kmer_index.get(encoded)
        if reference_index is not None and reference_index >= 0:
            votes[reference_index] += 1
    if not votes:
        return None, 0, 0
    ranked = votes.most_common(2)
    best_reference, best_votes = ranked[0]
    second_votes = ranked[1][1] if len(ranked) > 1 else 0
    return best_reference, best_votes, second_votes


def orient_and_map_read(
    sequence: str,
    kmer_index: dict[int, int],
) -> tuple[int, str, int] | None:
    normalized = sequence.strip().upper()
    if len(normalized) < 60:
        return None
    forward = _vote_reference(normalized, kmer_index)
    reverse_sequence = reverse_complement(normalized)
    reverse = _vote_reference(reverse_sequence, kmer_index)

    choices = [
        (forward[1], forward[1] - forward[2], 0, forward[0], normalized),
        (reverse[1], reverse[1] - reverse[2], 1, reverse[0], reverse_sequence),
    ]
    best_votes, margin, orientation_rank, reference_index, oriented = max(choices)
    del orientation_rank
    if (
        reference_index is None
        or best_votes < MAPPING_MIN_VOTES
        or margin < MAPPING_MIN_MARGIN
    ):
        return None
    return reference_index, oriented, best_votes


def _best_hamming_anchor(
    sequence: str,
    anchor: str,
    *,
    start: int,
    stop: int,
) -> tuple[int, int] | None:
    if stop - start < len(anchor):
        return None
    best: tuple[int, int] | None = None
    for position in range(start, stop - len(anchor) + 1):
        window = sequence[position : position + len(anchor)]
        mismatches = sum(a != b for a, b in zip(window, anchor, strict=True))
        candidate = (mismatches, position)
        if best is None or candidate < best:
            best = candidate
    return best


def trim_universal_primers(oriented_read: str) -> str | None:
    """Crop one oriented read using only the two universal primer sequences."""
    prefix_anchor = FORWARD_PRIMER[:PRIMER_ANCHOR_LENGTH]
    suffix_anchor = REVERSE_PRIMER[-PRIMER_ANCHOR_LENGTH:]
    prefix = _best_hamming_anchor(
        oriented_read,
        prefix_anchor,
        start=0,
        stop=min(len(oriented_read), PRIMER_SEARCH_SPAN),
    )
    suffix_start = max(0, len(oriented_read) - PRIMER_SEARCH_SPAN)
    suffix = _best_hamming_anchor(
        oriented_read,
        suffix_anchor,
        start=suffix_start,
        stop=len(oriented_read),
    )
    if prefix is None or suffix is None:
        return None
    prefix_mismatches, start = prefix
    suffix_mismatches, suffix_anchor_start = suffix
    if (
        prefix_mismatches > PRIMER_MAX_MISMATCHES
        or suffix_mismatches > PRIMER_MAX_MISMATCHES
    ):
        return None
    end = suffix_anchor_start + PRIMER_ANCHOR_LENGTH
    if end <= start or end - start < 120 or end - start > 260:
        return None
    cropped = oriented_read[start:end]
    if any(base not in "ACGT" for base in cropped):
        return None
    return cropped


def _iter_fastq_records(handle: BinaryIO) -> Iterable[tuple[str, str]]:
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        if not sequence or not plus or not quality:
            raise ValueError("truncated FASTQ record in LCRC archive")
        header_text = header.decode("utf-8", errors="strict").strip()
        if not header_text.startswith("@"):
            raise ValueError("invalid FASTQ header")
        yield header_text[1:].split()[0], sequence.decode("utf-8").strip()


def bin_ont_reads(
    archive_path: Path,
    references: list[str],
) -> tuple[list[list[str]], dict[str, Any]]:
    kmer_index = build_unique_kmer_index(references)
    clusters: list[list[str]] = [[] for _ in references]
    seen_ids: set[str] = set()
    stats: Counter[str] = Counter()

    with tarfile.open(archive_path, "r:gz") as archive:
        members = sorted(
            (
                member
                for member in archive.getmembers()
                if member.isfile()
                and member.name.lower().endswith((".fastq", ".fq"))
            ),
            key=lambda member: member.name,
        )
        stats["fastq_members"] = len(members)
        if not members:
            raise ValueError("LCRC ONT archive contains no FASTQ members")

        for member in members:
            extracted = archive.extractfile(member)
            if extracted is None:
                continue
            for read_id, sequence in _iter_fastq_records(extracted):
                stats["records_seen"] += 1
                if read_id in seen_ids:
                    stats["duplicate_read_ids"] += 1
                    continue
                seen_ids.add(read_id)
                if any(base not in "ACGTNacgtn" for base in sequence):
                    stats["non_dna_reads"] += 1
                    continue
                mapped = orient_and_map_read(sequence, kmer_index)
                if mapped is None:
                    stats["unmapped_or_ambiguous"] += 1
                    continue
                reference_index, oriented, votes = mapped
                stats["mapped_by_unique_kmers"] += 1
                stats["mapping_votes_total"] += votes
                cropped = trim_universal_primers(oriented)
                if cropped is None:
                    stats["primer_trim_rejected"] += 1
                    continue
                clusters[reference_index].append(cropped)
                stats["accepted_reads"] += 1

    stats["references_with_reads"] = sum(bool(cluster) for cluster in clusters)
    stats["references_with_5_reads"] = sum(len(cluster) >= 5 for cluster in clusters)
    stats["references_with_10_reads"] = sum(len(cluster) >= 10 for cluster in clusters)
    return clusters, dict(stats)


def _rank_key(seed: int, *parts: object) -> bytes:
    material = ":".join([str(seed), *(str(part) for part in parts)])
    return hashlib.sha256(material.encode("utf-8")).digest()


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
            f"only {len(eligible)} references have at least {max_coverage} accepted reads; "
            f"cannot select {size}"
        )
    selected = sorted(eligible, key=lambda index: _rank_key(seed, index))[:size]
    records: list[dict[str, Any]] = []
    for index in selected:
        ranked_reads = sorted(
            enumerate(clusters[index]),
            key=lambda item: _rank_key(seed, index, item[0], item[1]),
        )
        records.append(
            {
                "cluster_index": index + 1,
                "reference": references[index],
                "available_reads": len(clusters[index]),
                "ranked_reads": [read for _, read in ranked_reads],
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


def reconstruct_confidence_fusion(reads: list[str]) -> str:
    config = FROZEN_FUSION_CONFIG
    return confidence_fusion_trace_consensus(
        reads,
        target_length=REFERENCE_LENGTH,
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


def run_worker(subset_json: Path, worker_output: Path) -> None:
    payload = json.loads(subset_json.read_text(encoding="utf-8"))
    coverage = int(payload["coverage"])
    rows: list[dict[str, Any]] = []
    for record in payload["records"]:
        reference = str(record["reference"])
        reconstruction = reconstruct_confidence_fusion(
            [str(read) for read in record["reads"]]
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
    failure_classes: dict[str, int] = {}
    for row in rows:
        label = str(row["error_class"])
        failure_classes[label] = failure_classes.get(label, 0) + 1
    references = "".join(str(row["reference"]) for row in rows)
    reconstructions = "".join(str(row["reconstruction"]) for row in rows)
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
    reference_path = Path(args.references)
    archive_path = Path(args.read_archive)
    bbs_bin = Path(args.bbs_bin)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    timing_dir = output_dir / "timing"
    timing_dir.mkdir(exist_ok=True)

    references = load_references(reference_path)
    started = time.perf_counter()
    clusters, binning_stats = bin_ont_reads(archive_path, references)
    binning_seconds = time.perf_counter() - started

    development = select_split(
        references,
        clusters,
        size=args.development_size,
        seed=args.development_seed,
        max_coverage=max(COVERAGES),
    )
    development_ids = {int(record["cluster_index"]) - 1 for record in development}
    held_out = select_split(
        references,
        clusters,
        size=args.held_out_size,
        seed=args.held_out_seed,
        max_coverage=max(COVERAGES),
        excluded_indices=development_ids,
    )
    held_out_ids = {int(record["cluster_index"]) for record in held_out}
    if {
        int(record["cluster_index"]) for record in development
    } & held_out_ids:
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
        dev_started = time.perf_counter()
        rows: list[dict[str, Any]] = []
        for record in _coverage_records(development, coverage):
            reference = str(record["reference"])
            reconstruction = reconstruct_confidence_fusion(
                [str(read) for read in record["reads"]]
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
            summarize_rows(rows, time.perf_counter() - dev_started, None)
        )

    summaries: list[dict[str, Any]] = []
    raw_rows: list[dict[str, Any]] = []
    bbs_repeat_summaries: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []

    for coverage in COVERAGES:
        coverage_records = _coverage_records(held_out, coverage)
        subset_json = output_dir / f"subset-coverage-{coverage}.json"
        subset_json.write_text(
            json.dumps({"coverage": coverage, "records": coverage_records}),
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
                    str(REFERENCE_LENGTH),
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
            "paper_doi": SOURCE_PAPER_DOI,
            "sra_accession": SRA_ACCESSION,
            "repository_license": REPOSITORY_LICENSE,
            "reference_path": REFERENCE_PATH,
            "reference_git_blob": REFERENCE_GIT_BLOB,
            "reference_sha256": _sha256_file(reference_path),
            "read_archive_path": READ_ARCHIVE_PATH,
            "read_archive_git_blob": READ_ARCHIVE_GIT_BLOB,
            "read_archive_sha256": _sha256_file(archive_path),
            "reference_count": len(references),
            "reference_length": REFERENCE_LENGTH,
            "payload_length": REFERENCE_LENGTH
            - len(FORWARD_PRIMER)
            - len(REVERSE_PRIMER),
            "synthesis": "Twist Bioscience high-fidelity HFS-Pool-11.7K",
            "sequencing": (
                "Oxford Nanopore MinION R10.4.1, Ligation Sequencing Kit V14 "
                "(SQK-LSK114), FAST real-time read archive"
            ),
        },
        "association": {
            "reference_assisted_for_binning_only": True,
            "mapping_rule": (
                f"unique exact {MAPPING_K}-mer votes from the 160-nt payload; "
                f"minimum {MAPPING_MIN_VOTES} votes and vote margin "
                f"{MAPPING_MIN_MARGIN}"
            ),
            "orientation": "highest unique-kmer vote between read and reverse complement",
            "cropping": (
                f"universal-primer anchors only: {PRIMER_ANCHOR_LENGTH}-nt anchors, "
                f"at most {PRIMER_MAX_MISMATCHES} substitutions per anchor"
            ),
            "reference_used_during_reconstruction": False,
            "binning_seconds": round(binning_seconds, 6),
            "binning_stats": binning_stats,
        },
        "oligoark": {
            "algorithm": "confidence_fusion_trace_consensus",
            "frozen_config": FROZEN_FUSION_CONFIG,
            "parameter_changes_from_cnr_grass": "none",
            "mechanical_target_length": REFERENCE_LENGTH,
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
            "read_ranking": "deterministic SHA-256 ranking inside each mapped cluster",
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
            "reference-level reconstruction conditional on transparent reference-assisted "
            "read binning; not end-to-end archive decoding"
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
    parser.add_argument("--read-archive")
    parser.add_argument("--bbs-bin")
    parser.add_argument("--output-dir", default="external-lcrc-ont-results")
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
    if not args.references or not args.read_archive or not args.bbs_bin:
        raise ValueError("--references, --read-archive and --bbs-bin are required")
    if args.bbs_repeats < 1:
        raise ValueError("--bbs-repeats must be positive")
    run_benchmark(args)


if __name__ == "__main__":
    main()

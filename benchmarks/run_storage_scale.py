"""Scalable OligoArk archive benchmark.

The benchmark distinguishes software archive scaling from physical-read reconstruction.
The 100 MiB acceptance profile emits the milestone sentence only after SHA-256 verified
clean and controlled-loss recovery both succeed.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import platform
import random
import resource
import shutil
import subprocess
import sys
import time
import zipfile
import zlib
from pathlib import Path
from typing import Any

from oligoark import __version__
from oligoark.archive import ArchiveConfig
from oligoark.profiles import physical_strand_profile
from oligoark.streaming import (
    StreamingFaultProfile,
    archive_file_streaming,
    recover_file_streaming,
)

KIB = 1024
MIB = 1024 * 1024
GIB = 1024 * MIB
SCALE_SIZES = (1 * KIB, 64 * KIB, 1 * MIB, 10 * MIB, 100 * MIB, 1 * GIB)
SCHEMES = ("none", "xor", "fountain", "hybrid")
SCALE_PROFILE = "scale-1024"
PHYSICAL_PROFILE = "oligoark-248"
PHYSICAL_MATRIX_SIZE = 64 * KIB
MILESTONE_SIZE = 1 * GIB
MILESTONE_STATEMENT = (
    "OligoArk successfully archives and SHA-256 recovers a 1 GiB heterogeneous "
    "dataset under clean, 1% and 5% controlled strand loss, while preserving "
    "bounded-memory streaming measurements and exact-recovery evidence."
)

FAULTS: dict[str, StreamingFaultProfile] = {
    "clean": StreamingFaultProfile(seed=2026),
    "dropout-1": StreamingFaultProfile(dropout_rate=0.01, seed=2026),
    "dropout-5": StreamingFaultProfile(dropout_rate=0.05, seed=2026),
    "substitution-low": StreamingFaultProfile(substitution_rate=0.0005, seed=2026),
    "indel-low": StreamingFaultProfile(
        insertion_rate=0.0001,
        deletion_rate=0.0001,
        seed=2026,
    ),
    "mixed": StreamingFaultProfile(
        dropout_rate=0.01,
        substitution_rate=0.0002,
        insertion_rate=0.00005,
        deletion_rate=0.00005,
        seed=2026,
    ),
}


def _benchmark_config(strand_profile: str, scheme: str) -> ArchiveConfig:
    """Resolve either the fast scale profile or an RS-enabled physical profile."""
    if strand_profile == SCALE_PROFILE:
        return ArchiveConfig(
            chunk_size=237,
            rs_nsym=0,
            parity_group_size=8,
            adaptive_masks=False,
            redundancy_scheme=scheme,
            fountain_redundancy=0.25,
            min_gc_fraction=0.0,
            max_gc_fraction=1.0,
            max_homopolymer=1024,
            mask_search_limit=1,
        )
    return physical_strand_profile(strand_profile).with_scheme(scheme).to_archive_config(
        fountain_redundancy=0.25,
    )


def _heterogeneous_block() -> bytes:
    rng = random.Random(2026)
    text = (
        b"OligoArk heterogeneous storage benchmark\n"
        b"UTF-8 text, structured records, source code, image bytes, compressed bytes.\n"
    ) * 128
    jsonl = b"".join(
        json.dumps(
            {
                "id": index,
                "kind": "archive-record",
                "active": index % 3 != 0,
                "score": round(math.sin(index) + 1.0, 6),
            },
            sort_keys=True,
        ).encode()
        + b"\n"
        for index in range(256)
    )
    csv_data = b"id,category,value\n" + b"".join(
        f"{index},sample,{index * 17}\n".encode()
        for index in range(256)
    )
    source = (
        b"def recover(payload: bytes) -> bytes:\n"
        b"    # deterministic source-code fixture\n"
        b"    return payload\n"
    ) * 256
    image_pixels = bytes(rng.randrange(256) for _ in range(128 * 128 * 3))
    ppm = b"P6\n128 128\n255\n" + image_pixels
    binary = bytes(rng.randrange(256) for _ in range(96 * KIB))
    compressed = zlib.compress(binary, level=9)

    mixed_archive = io.BytesIO()
    with zipfile.ZipFile(mixed_archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        def add_member(name: str, payload: bytes) -> None:
            # ZipFile.writestr(name, data) otherwise embeds the current local time, which
            # makes the supposedly deterministic heterogeneous fixture differ per worker.
            info = zipfile.ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            bundle.writestr(info, payload)

        add_member("notes.txt", text)
        add_member("records.jsonl", jsonl)
        add_member("module.py", source)
        add_member("image.ppm", ppm)
        add_member("payload.bin", binary)
    return (
        text
        + jsonl
        + csv_data
        + source
        + ppm
        + binary
        + compressed
        + mixed_archive.getvalue()
    )


def create_heterogeneous_payload(path: Path, size: int) -> None:
    block = _heterogeneous_block()
    remaining = size
    with path.open("wb") as handle:
        while remaining:
            part = block[: min(remaining, len(block))]
            handle.write(part)
            remaining -= len(part)


def _peak_rss_mib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return value / (1024 * 1024)
    return value / 1024


def _fault_to_dict(fault: StreamingFaultProfile) -> dict[str, object]:
    return {
        "dropout_rate": fault.dropout_rate,
        "substitution_rate": fault.substitution_rate,
        "insertion_rate": fault.insertion_rate,
        "deletion_rate": fault.deletion_rate,
        "seed": fault.seed,
        "controlled_dropout": fault.controlled_dropout,
    }


def _worker(
    size: int,
    scheme: str,
    fault_name: str,
    strand_profile: str,
    workdir: Path,
) -> dict[str, object]:
    baseline_rss = _peak_rss_mib()
    workdir.mkdir(parents=True, exist_ok=True)
    source = workdir / "heterogeneous.bin"
    archive = workdir / "archive.oab"
    recovered = workdir / "recovered.bin"
    create_heterogeneous_payload(source, size)

    config = _benchmark_config(strand_profile, scheme)
    fault = FAULTS[fault_name]
    encode_started = time.perf_counter()
    stats = archive_file_streaming(source, archive, config)
    encode_seconds = time.perf_counter() - encode_started

    decode_started = time.perf_counter()
    error: str | None = None
    try:
        report = recover_file_streaming(archive, recovered, fault=fault, strict=False)
    except (OSError, ValueError) as exc:
        report = None
        error = f"{type(exc).__name__}: {exc}"
    decode_seconds = time.perf_counter() - decode_started
    source_mib = size / MIB
    recovered_ok = bool(report and report.verified_sha256)
    encoded_ideal = max(1, size * 4)

    row: dict[str, object] = {
        "size_bytes": size,
        "size_mib": round(source_mib, 6),
        "payload_kind": "deterministic-heterogeneous-mixed",
        "scheme": scheme,
        "fault": fault_name,
        "strand_profile": strand_profile,
        "chunk_size": config.chunk_size,
        "rs_nsym": config.rs_nsym,
        "fault_config": _fault_to_dict(fault),
        "sha256_verified": recovered_ok,
        "exact_recovery_rate": 1.0 if recovered_ok else 0.0,
        "encode_seconds": round(encode_seconds, 6),
        "decode_seconds": round(decode_seconds, 6),
        "total_seconds": round(encode_seconds + decode_seconds, 6),
        "encode_throughput_mib_s": round(
            source_mib / encode_seconds if encode_seconds else 0.0,
            6,
        ),
        "decode_throughput_mib_s": round(
            source_mib / decode_seconds if decode_seconds else 0.0,
            6,
        ),
        "peak_rss_mib": round(_peak_rss_mib(), 3),
        "baseline_rss_mib": round(baseline_rss, 3),
        "rss_growth_mib": round(max(0.0, _peak_rss_mib() - baseline_rss), 3),
        "encoded_nucleotides": stats.encoded_nucleotides,
        "logical_bits_per_nucleotide": stats.logical_bits_per_nucleotide,
        "nucleotide_overhead_vs_2bit_ideal": round(
            stats.encoded_nucleotides / encoded_ideal,
            6,
        ),
        "strand_count": stats.strand_count,
        "data_strands": stats.data_strands,
        "parity_strands": stats.parity_strands,
        "fountain_strands": stats.fountain_strands,
        "strand_redundancy_ratio": stats.strand_redundancy_ratio,
        "archive_size_bytes": stats.archive_size_bytes,
        "archive_size_overhead_ratio": stats.archive_size_overhead_ratio,
        "constraint_pass_rate_sampled": stats.constraint_pass_rate,
        "constraint_sample_size": stats.constraint_sample_size,
        "error": error,
    }
    if report:
        row.update(
            {
                "missing_data_strands": report.missing_data_strands,
                "dropped_records": report.dropped_records,
                "undecodable_records": report.undecodable_records,
                "xor_recovered_strands": report.xor_recovered_strands,
                "fountain_recovered_strands": report.fountain_recovered_strands,
                "output_sha256": report.output_sha256,
            }
        )
    return row


def _cases(profile: str) -> list[tuple[int, str, str, str]]:
    if profile == "ci":
        return [
            (1 * KIB, scheme, fault, SCALE_PROFILE)
            for scheme in SCHEMES
            for fault in ("clean", "dropout-5")
        ] + [
            (64 * KIB, "xor", "clean", SCALE_PROFILE),
            (64 * KIB, "xor", "dropout-5", SCALE_PROFILE),
        ]

    acceptance = [
        (size, "xor", "clean", SCALE_PROFILE)
        for size in SCALE_SIZES
    ] + [
        (MILESTONE_SIZE, "xor", "dropout-1", SCALE_PROFILE),
        (MILESTONE_SIZE, "xor", "dropout-5", SCALE_PROFILE),
    ]
    physical_matrix = [
        (PHYSICAL_MATRIX_SIZE, scheme, fault, PHYSICAL_PROFILE)
        for scheme in SCHEMES
        for fault in FAULTS
    ]

    if profile == "acceptance":
        return acceptance
    if profile == "physical":
        return physical_matrix
    if profile == "storage":
        return acceptance + physical_matrix

    # Full is intentionally exhaustive and may be expensive. The 1 GiB point is kept
    # to the mandatory XOR acceptance cases; all-scheme/all-fault sweeps stop at 100 MiB.
    exhaustive = [
        (size, scheme, fault, SCALE_PROFILE)
        for size in SCALE_SIZES[:-1]
        for scheme in SCHEMES
        for fault in FAULTS
    ]
    return acceptance + physical_matrix + exhaustive


def _log_slope(rows: list[dict[str, object]], metric: str) -> float | None:
    pairs = [
        (float(row["size_bytes"]), float(row[metric]))
        for row in rows
        if bool(row["sha256_verified"])
        and metric in row
        and float(row[metric]) > 0
    ]
    unique = {}
    for size, value in pairs:
        unique[size] = value
    if len(unique) < 2:
        return None
    xs = [math.log(size) for size in sorted(unique)]
    ys = [math.log(unique[size]) for size in sorted(unique)]
    xmean = sum(xs) / len(xs)
    ymean = sum(ys) / len(ys)
    denom = sum((x - xmean) ** 2 for x in xs)
    if denom == 0:
        return None
    return sum((x - xmean) * (y - ymean) for x, y in zip(xs, ys, strict=True)) / denom


def _classify_memory(slope: float | None) -> str:
    if slope is None:
        return "insufficient-data"
    if slope < 0.25:
        return "bounded"
    if slope <= 1.25:
        return "approximately-linear"
    return "super-linear"


def _classify_runtime(slope: float | None) -> str:
    if slope is None:
        return "insufficient-data"
    if slope <= 1.25:
        return "linear-or-better"
    return "super-linear"


def _summary(rows: list[dict[str, object]]) -> dict[str, object]:
    clean_xor = [
        row
        for row in rows
        if row["scheme"] == "xor"
        and row["fault"] == "clean"
        and row.get("strand_profile") == SCALE_PROFILE
    ]
    memory_slope = _log_slope(clean_xor, "peak_rss_mib")
    runtime_slope = _log_slope(clean_xor, "total_seconds")
    milestone_rows = [
        row
        for row in rows
        if int(row["size_bytes"]) == MILESTONE_SIZE
        and row["scheme"] == "xor"
        and row["fault"] in {"clean", "dropout-1", "dropout-5"}
        and row.get("strand_profile") == SCALE_PROFILE
    ]
    milestone = (
        len(milestone_rows) == 3
        and all(bool(row["sha256_verified"]) for row in milestone_rows)
    )
    failures = [
        {
            "size_bytes": row["size_bytes"],
            "scheme": row["scheme"],
            "fault": row["fault"],
            "error": row["error"],
            "missing_data_strands": row.get("missing_data_strands"),
        }
        for row in rows
        if not bool(row["sha256_verified"])
    ]
    return {
        "milestone_achieved": milestone,
        "milestone_statement": MILESTONE_STATEMENT if milestone else None,
        "milestone_failure": None
        if milestone
        else (
            "1 GiB clean, controlled 1% dropout and controlled 5% dropout have not "
            "all passed the SHA-256 gate in this result set."
        ),
        "memory_log_log_slope_xor_clean": memory_slope,
        "memory_scaling_xor_clean": _classify_memory(memory_slope),
        "runtime_log_log_slope_xor_clean": runtime_slope,
        "runtime_scaling_xor_clean": _classify_runtime(runtime_slope),
        "failure_count": len(failures),
        "failures": failures,
        "physical_matrix_rows": sum(
            row.get("strand_profile") == PHYSICAL_PROFILE for row in rows
        ),
        "physical_matrix_failures": sum(
            row.get("strand_profile") == PHYSICAL_PROFILE
            and not bool(row.get("sha256_verified"))
            for row in rows
        ),
    }


def _write_csv(rows: list[dict[str, object]], path: Path) -> None:
    if not rows:
        return
    flattened: list[dict[str, object]] = []
    for row in rows:
        flattened.append(
            {
                key: json.dumps(value, sort_keys=True)
                if isinstance(value, (dict, list))
                else value
                for key, value in row.items()
            }
        )
    fieldnames = sorted({key for row in flattened for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(flattened)


def _write_plots(rows: list[dict[str, object]], output: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    clean_xor = [
        row
        for row in rows
        if row["scheme"] == "xor"
        and row["fault"] == "clean"
        and row.get("strand_profile") == SCALE_PROFILE
    ]
    if len(clean_xor) < 2:
        return
    clean_xor.sort(key=lambda row: int(row["size_bytes"]))
    sizes = [float(row["size_mib"]) for row in clean_xor]

    plt.figure(figsize=(8, 4))
    plt.plot(sizes, [float(row["total_seconds"]) for row in clean_xor], marker="o")
    plt.xlabel("Payload size (MiB)")
    plt.ylabel("Encode + recover runtime (s)")
    plt.tight_layout()
    plt.savefig(output / "scale_runtime.png", dpi=160)
    plt.close()

    plt.figure(figsize=(8, 4))
    plt.plot(sizes, [float(row["peak_rss_mib"]) for row in clean_xor], marker="o")
    plt.xlabel("Payload size (MiB)")
    plt.ylabel("Peak RSS (MiB)")
    plt.tight_layout()
    plt.savefig(output / "scale_peak_memory.png", dpi=160)
    plt.close()

    plt.figure(figsize=(8, 4))
    plt.plot(
        sizes,
        [float(row["logical_bits_per_nucleotide"]) for row in clean_xor],
        marker="o",
    )
    plt.xlabel("Payload size (MiB)")
    plt.ylabel("Logical bits per nucleotide")
    plt.tight_layout()
    plt.savefig(output / "scale_density.png", dpi=160)
    plt.close()


def _run_isolated(
    size: int,
    scheme: str,
    fault: str,
    strand_profile: str,
    output: Path,
    ordinal: int,
) -> dict[str, object]:
    case_dir = output / "work" / (
        f"{ordinal:03d}-{size}-{scheme}-{fault}-{strand_profile}"
    )
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--size",
        str(size),
        "--scheme",
        scheme,
        "--fault",
        fault,
        "--strand-profile",
        strand_profile,
        "--workdir",
        str(case_dir),
    ]
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    try:
        if completed.returncode != 0:
            return {
                "size_bytes": size,
                "size_mib": round(size / MIB, 6),
                "scheme": scheme,
                "fault": fault,
                "strand_profile": strand_profile,
                "sha256_verified": False,
                "error": completed.stderr.strip() or completed.stdout.strip(),
            }
        return json.loads(completed.stdout)
    finally:
        shutil.rmtree(case_dir, ignore_errors=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--profile",
        choices=("ci", "acceptance", "physical", "storage", "full"),
        default="ci",
    )
    parser.add_argument("--output", type=Path, default=Path("storage-scale-results"))
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--size", type=int)
    parser.add_argument("--scheme", choices=SCHEMES)
    parser.add_argument("--fault", choices=tuple(FAULTS))
    parser.add_argument(
        "--strand-profile",
        choices=(SCALE_PROFILE, "oligoark-152", "oligoark-200", PHYSICAL_PROFILE),
    )
    parser.add_argument("--workdir", type=Path)
    args = parser.parse_args()

    if args.worker:
        if (
            args.size is None
            or args.scheme is None
            or args.fault is None
            or args.strand_profile is None
            or args.workdir is None
        ):
            raise SystemExit(
                "worker requires --size, --scheme, --fault, --strand-profile and --workdir"
            )
        print(
            json.dumps(
                _worker(
                    args.size,
                    args.scheme,
                    args.fault,
                    args.strand_profile,
                    args.workdir,
                )
            )
        )
        return

    args.output.mkdir(parents=True, exist_ok=True)
    rows = [
        _run_isolated(size, scheme, fault, strand_profile, args.output, ordinal)
        for ordinal, (size, scheme, fault, strand_profile) in enumerate(_cases(args.profile))
    ]
    summary = _summary(rows)
    metadata: dict[str, Any] = {
        "oligoark_version": __version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "profile": args.profile,
        "sizes_bytes": sorted({int(row["size_bytes"]) for row in rows}),
        "strand_profiles": sorted(
            {str(row.get("strand_profile", "")) for row in rows if row.get("strand_profile")}
        ),
        "claim_scope": (
            "software archive scaling plus controlled software-channel faults; "
            "the scale-1024 profile is a streaming systems test, while oligoark-248 "
            "is an RS-enabled realistic-strand software test; physical-read reconstruction "
            "and wet-lab archive storage remain separate evidence classes"
        ),
        "payload_description": (
            "deterministic mixture of UTF-8 text, JSONL, Python source, PPM image bytes, "
            "CSV, random binary, compressed binary and a ZIP mixed-file archive"
        ),
    }
    (args.output / "results.json").write_text(
        json.dumps(rows, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _write_csv(rows, args.output / "results.csv")
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (args.output / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _write_plots(rows, args.output)
    print(json.dumps({"metadata": metadata, "summary": summary}, indent=2))


if __name__ == "__main__":
    main()

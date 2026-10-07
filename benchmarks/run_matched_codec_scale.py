"""Matched 152-nt codec benchmark with reusable immutable encoded artifacts.

A method/size pair is encoded exactly once. Each channel condition is then evaluated in an
isolated worker that loads the same read-only encoded artifact. This prevents repeated encode
work, keeps failures/timeouts visible, and allows a slow noisy condition to time out without
erasing completed clean/dropout evidence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import pickle
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from run_dna_fountain_baseline import (
    CONDITIONS,
    KIB,
    MIB,
    _channel_sequences,
    _payload,
    _write_plots,
)

from oligoark import __version__
from oligoark.archive import archive_bytes, archive_statistics, recover_bytes
from oligoark.baselines import (
    DnaFountainBaselineArchive,
    DnaFountainBaselineConfig,
    RotatingTernaryBaselineArchive,
    RotatingTernaryBaselineConfig,
    decode_dna_fountain_baseline,
    decode_rotating_ternary_baseline,
    encode_dna_fountain_baseline,
    encode_rotating_ternary_baseline,
)
from oligoark.experiments import wilson_interval
from oligoark.profiles import physical_strand_profile

METHODS = (
    "oligoark-compact-hybrid",
    "oligoark-efficient-hybrid-v1",
    "dna-fountain-cleanroom",
    "goldman-rotating-xor",
)
FULL_SIZES = (1 * KIB, 64 * KIB, 1 * MIB, 10 * MIB)
SCALE_SIZES = (1 * KIB, 64 * KIB, 1 * MIB, 10 * MIB, 100 * MIB)
RESEARCH_SIZES = (10 * MIB, 100 * MIB)


def _peak_rss_mib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return value / (1024 * 1024)
    return value / 1024


def _artifact_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * MIB), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json_atomic(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def _encode_artifact(
    method: str,
    size: int,
    redundancy: float,
    artifact_path: Path,
) -> dict[str, object]:
    payload_started = time.perf_counter()
    payload = _payload(size)
    payload_seconds = time.perf_counter() - payload_started
    payload_sha256 = hashlib.sha256(payload).hexdigest()

    encode_started = time.perf_counter()
    extra: dict[str, object] = {}
    if method in {"oligoark-compact-hybrid", "oligoark-efficient-hybrid-v1"}:
        profile_name = (
            "oligoark-152-efficient-v1"
            if method == "oligoark-efficient-hybrid-v1"
            else "oligoark-152-compact"
        )
        profile = physical_strand_profile(profile_name).with_scheme("hybrid")
        xor_share = 1.0 / profile.parity_group_size
        fountain_share = max(0.0, redundancy - xor_share)
        config = profile.to_archive_config(fountain_redundancy=fountain_share)
        archive = archive_bytes(payload, config)
        stats = archive_statistics(archive)
        sequences = archive.strands
        data_units = int(archive.metadata["data_strands"])
        encoded_nucleotides = stats.encoded_nucleotides
        density = stats.logical_bits_per_nucleotide
        extra = {
            "redundancy_allocation": {
                "xor": round(xor_share, 6),
                "fountain": round(fountain_share, 6),
            },
            "profile": profile_name,
            "compact_framing": config.compact_framing,
            "compact_index_bytes": config.compact_index_bytes,
            "compact_typed_index": config.compact_typed_index,
            "rs_nsym": config.rs_nsym,
            "chunk_size": config.chunk_size,
            "fountain_max_degree": config.fountain_max_degree,
        }
    elif method == "dna-fountain-cleanroom":
        config = DnaFountainBaselineConfig(redundancy=redundancy)
        archive = encode_dna_fountain_baseline(payload, config)
        sequences = archive.sequences
        data_units = archive.chunk_count
        encoded_nucleotides = archive.encoded_nucleotides
        density = archive.logical_bits_per_nucleotide
        extra = {
            "accepted_droplet_attempts": archive.attempts,
            "chunk_size": config.chunk_size,
            "rs_nsym": config.rs_nsym,
        }
    elif method == "goldman-rotating-xor":
        config = RotatingTernaryBaselineConfig(redundancy=redundancy)
        archive = encode_rotating_ternary_baseline(payload, config)
        sequences = archive.sequences
        data_units = archive.data_count
        encoded_nucleotides = archive.encoded_nucleotides
        density = archive.logical_bits_per_nucleotide
        extra = {
            "chunk_size": config.chunk_size,
            "parity_group_size": config.parity_group_size,
        }
    else:
        raise ValueError(f"unknown method: {method}")
    encode_seconds = time.perf_counter() - encode_started

    encode_peak_rss_mib = _peak_rss_mib()
    record = {
        "method": method,
        "size_bytes": size,
        "redundancy_budget": redundancy,
        "payload_sha256": payload_sha256,
        "payload_generation_seconds": payload_seconds,
        "encode_seconds": encode_seconds,
        "encode_peak_rss_mib": encode_peak_rss_mib,
        "target_max_strand_nt": 152,
        "max_strand_nt": max(map(len, sequences)),
        "strand_count": len(sequences),
        "data_units": data_units,
        "measured_strand_redundancy_ratio": (
            (len(sequences) - data_units) / max(1, data_units)
        ),
        "encoded_nucleotides": encoded_nucleotides,
        "logical_bits_per_nucleotide": density,
        "archive": archive,
        "extra": extra,
    }

    artifact_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = artifact_path.with_suffix(artifact_path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        pickle.dump(record, handle, protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(artifact_path)
    os.chmod(artifact_path, 0o444)

    return {
        key: value
        for key, value in record.items()
        if key not in {"archive"}
    } | {
        "artifact_size_bytes": artifact_path.stat().st_size,
        "artifact_sha256": _artifact_sha256(artifact_path),
    }


def _load_artifact(path: Path) -> dict[str, object]:
    with path.open("rb") as handle:
        value = pickle.load(handle)
    if not isinstance(value, dict):
        raise ValueError("encoded artifact is not a mapping")
    return value


def _decode_adapter(
    method: str,
    record: dict[str, object],
) -> tuple[tuple[str, ...] | list[str], Any]:
    archive = record["archive"]
    if method in {"oligoark-compact-hybrid", "oligoark-efficient-hybrid-v1"}:
        return archive.strands, lambda reads: recover_bytes(archive, reads)
    if method == "dna-fountain-cleanroom":
        if not isinstance(archive, DnaFountainBaselineArchive):
            raise ValueError("DNA Fountain artifact type mismatch")
        return archive.sequences, lambda reads: decode_dna_fountain_baseline(archive, reads)
    if method == "goldman-rotating-xor":
        if not isinstance(archive, RotatingTernaryBaselineArchive):
            raise ValueError("Goldman artifact type mismatch")
        return archive.sequences, lambda reads: decode_rotating_ternary_baseline(archive, reads)
    raise ValueError(f"unknown method: {method}")


def _aggregate_condition(
    record: dict[str, object],
    condition_name: str,
    requested_trials: int,
    trial_results: list[dict[str, object]],
    *,
    artifact_size_bytes: int,
    artifact_sha256: str,
    timed_out: bool = False,
    timeout_seconds: int | None = None,
    error: str | None = None,
) -> dict[str, object]:
    completed_trials = len(trial_results)
    successes = sum(bool(item.get("sha256_verified")) for item in trial_results)
    failures = completed_trials - successes
    if completed_trials:
        low, high = wilson_interval(successes, completed_trials)
        recovery_rate: float | None = successes / completed_trials
        ci_low: float | None = low
        ci_high: float | None = high
    else:
        recovery_rate = None
        ci_low = None
        ci_high = None

    decode_seconds = sum(float(item.get("decode_seconds", 0.0)) for item in trial_results)
    channel_seconds = sum(float(item.get("channel_seconds", 0.0)) for item in trial_results)
    payload_size = int(record["size_bytes"])
    payload_mib = payload_size / MIB
    mean_decode = decode_seconds / completed_trials if completed_trials else 0.0
    encode_seconds = float(record["encode_seconds"])
    peak_rss = max(float(record["encode_peak_rss_mib"]), _peak_rss_mib())

    row: dict[str, object] = {
        "method": record["method"],
        "size_bytes": payload_size,
        "payload_sha256": record["payload_sha256"],
        "target_max_strand_nt": record["target_max_strand_nt"],
        "max_strand_nt": record["max_strand_nt"],
        "redundancy_budget": record["redundancy_budget"],
        "condition": condition_name,
        **CONDITIONS[condition_name],
        "trials": requested_trials,
        "completed_trials": completed_trials,
        "benchmark_complete": (
            completed_trials == requested_trials
            and not timed_out
            and error is None
        ),
        "successes": successes,
        "failures": failures,
        "recovery_rate": None if recovery_rate is None else round(recovery_rate, 6),
        "recovery_ci95_low": None if ci_low is None else round(ci_low, 6),
        "recovery_ci95_high": None if ci_high is None else round(ci_high, 6),
        "sha256_verified_success_definition": True,
        "strand_count": record["strand_count"],
        "data_units": record["data_units"],
        "measured_strand_redundancy_ratio": round(
            float(record["measured_strand_redundancy_ratio"]),
            6,
        ),
        "encoded_nucleotides": record["encoded_nucleotides"],
        "logical_bits_per_nucleotide": round(
            float(record["logical_bits_per_nucleotide"]),
            6,
        ),
        "payload_generation_seconds": round(
            float(record["payload_generation_seconds"]),
            6,
        ),
        "encode_seconds": round(encode_seconds, 6),
        "channel_seconds": round(channel_seconds, 6),
        "decode_seconds_total": round(decode_seconds, 6),
        "mean_decode_seconds": round(mean_decode, 6),
        "encode_throughput_mib_s": round(
            payload_mib / encode_seconds if encode_seconds else 0.0,
            6,
        ),
        "decode_throughput_mib_s": round(
            payload_mib / mean_decode if mean_decode else 0.0,
            6,
        ),
        "encode_peak_rss_mib": round(float(record["encode_peak_rss_mib"]), 3),
        "peak_rss_mib": round(peak_rss, 3),
        "encoded_artifact_size_bytes": artifact_size_bytes,
        "encoded_artifact_sha256": artifact_sha256,
        "trial_results": trial_results,
        "failed_trial_indexes": [
            int(item["trial"])
            for item in trial_results
            if not bool(item.get("sha256_verified"))
        ],
        **dict(record.get("extra", {})),
    }
    if timed_out:
        row["timed_out"] = True
        row["timeout_seconds"] = timeout_seconds
    if error is not None:
        row["error"] = error
    return row


def _condition_worker(
    artifact_path: Path,
    condition_name: str,
    trials: int,
    checkpoint: Path,
) -> dict[str, object]:
    record = _load_artifact(artifact_path)
    method = str(record["method"])
    sequences, decode = _decode_adapter(method, record)
    expected_sha256 = str(record["payload_sha256"])
    trial_results: list[dict[str, object]] = []

    for trial in range(trials):
        seed = 20_260_000 + trial
        channel_started = time.perf_counter()
        reads = _channel_sequences(sequences, condition_name, seed)
        channel_seconds = time.perf_counter() - channel_started

        decode_started = time.perf_counter()
        failure: str | None = None
        verified = False
        try:
            recovered = decode(reads)
            verified = hashlib.sha256(recovered).hexdigest() == expected_sha256
            if not verified:
                failure = "sha256-mismatch"
        except ValueError as exc:
            failure = f"decode-error: {exc}"
        decode_seconds = time.perf_counter() - decode_started
        trial_results.append(
            {
                "trial": trial,
                "seed": seed,
                "sha256_verified": verified,
                "failure": failure,
                "channel_seconds": round(channel_seconds, 6),
                "decode_seconds": round(decode_seconds, 6),
            }
        )
        _write_json_atomic(
            checkpoint,
            {
                "trial_results": trial_results,
                "peak_rss_mib": round(_peak_rss_mib(), 3),
            },
        )

    return _aggregate_condition(
        record,
        condition_name,
        trials,
        trial_results,
        artifact_size_bytes=artifact_path.stat().st_size,
        artifact_sha256=_artifact_sha256(artifact_path),
    )


def _run_command(
    command: list[str],
    timeout_seconds: int,
) -> subprocess.CompletedProcess[str] | None:
    try:
        return subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired:
        return None


def _prepare_artifact(
    method: str,
    size: int,
    redundancy: float,
    artifact_path: Path,
    timeout_seconds: int,
) -> tuple[dict[str, object] | None, str | None, bool]:
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--encode-worker",
        "--method",
        method,
        "--size",
        str(size),
        "--redundancy",
        str(redundancy),
        "--artifact",
        str(artifact_path),
    ]
    completed = _run_command(command, timeout_seconds)
    if completed is None:
        return None, f"encode worker exceeded {timeout_seconds}s timeout", True
    if completed.returncode:
        error = completed.stderr.strip() or completed.stdout.strip()
        return None, error, False
    value = json.loads(completed.stdout)
    if not isinstance(value, dict):
        raise ValueError("encode worker did not return metadata")
    return value, None, False


def _encode_failure_rows(
    method: str,
    size: int,
    redundancy: float,
    conditions: tuple[str, ...],
    trials: int,
    error: str,
    timed_out: bool,
    timeout_seconds: int,
) -> list[dict[str, object]]:
    return [
        {
            "method": method,
            "size_bytes": size,
            "redundancy_budget": redundancy,
            "condition": condition,
            **CONDITIONS[condition],
            "trials": trials,
            "completed_trials": 0,
            "successes": 0,
            "failures": 0,
            "benchmark_complete": False,
            "stage": "encode",
            "timed_out": timed_out,
            "timeout_seconds": timeout_seconds if timed_out else None,
            "error": error,
        }
        for condition in conditions
    ]


def _run_condition_isolated(
    *,
    artifact_path: Path,
    encode_meta: dict[str, object],
    condition_name: str,
    trials: int,
    timeout_seconds: int,
    checkpoint_path: Path,
) -> dict[str, object]:
    checkpoint_path.unlink(missing_ok=True)
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--condition-worker",
        "--artifact",
        str(artifact_path),
        "--condition",
        condition_name,
        "--trials",
        str(trials),
        "--checkpoint",
        str(checkpoint_path),
    ]
    completed = _run_command(command, timeout_seconds)
    if completed is not None and completed.returncode == 0:
        value = json.loads(completed.stdout)
        if not isinstance(value, dict):
            raise ValueError("condition worker did not return a row")
        return value

    partial: list[dict[str, object]] = []
    if checkpoint_path.exists():
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        raw_trials = checkpoint.get("trial_results", [])
        if isinstance(raw_trials, list):
            partial = [item for item in raw_trials if isinstance(item, dict)]

    record = dict(encode_meta)
    record["extra"] = encode_meta.get("extra", {})
    if completed is None:
        return _aggregate_condition(
            record,
            condition_name,
            trials,
            partial,
            artifact_size_bytes=int(encode_meta["artifact_size_bytes"]),
            artifact_sha256=str(encode_meta["artifact_sha256"]),
            timed_out=True,
            timeout_seconds=timeout_seconds,
            error=(
                f"condition worker exceeded {timeout_seconds}s after "
                f"{len(partial)}/{trials} completed trials"
            ),
        )

    error = completed.stderr.strip() or completed.stdout.strip()
    return _aggregate_condition(
        record,
        condition_name,
        trials,
        partial,
        artifact_size_bytes=int(encode_meta["artifact_size_bytes"]),
        artifact_sha256=str(encode_meta["artifact_sha256"]),
        error=error,
    )


def _write_csv(rows: list[dict[str, object]], path: Path) -> None:
    flattened = [
        {
            key: json.dumps(value, sort_keys=True)
            if isinstance(value, (dict, list))
            else value
            for key, value in row.items()
        }
        for row in rows
    ]
    fields = sorted({key for row in flattened for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flattened)


def _summary(rows: list[dict[str, object]], conditions: tuple[str, ...]) -> dict[str, object]:
    sizes = sorted({int(row["size_bytes"]) for row in rows})
    common_sizes: list[int] = []
    for size in sizes:
        complete = True
        for method in METHODS:
            selected = [
                row
                for row in rows
                if row.get("method") == method and int(row["size_bytes"]) == size
            ]
            if (
                len(selected) != len(conditions)
                or any(not bool(row.get("benchmark_complete")) for row in selected)
                or any("error" in row for row in selected)
            ):
                complete = False
                break
        if complete:
            common_sizes.append(size)

    largest = common_sizes[-1] if common_sizes else None
    comparison: dict[str, object] = {}
    if largest is not None:
        clean = [
            row
            for row in rows
            if int(row["size_bytes"]) == largest
            and row.get("condition") == "clean"
            and "logical_bits_per_nucleotide" in row
        ]
        dropout = [
            row
            for row in rows
            if int(row["size_bytes"]) == largest
            and row.get("condition") == "dropout-5"
            and row.get("recovery_rate") is not None
        ]
        if clean:
            winner = max(clean, key=lambda row: float(row["logical_bits_per_nucleotide"]))
            comparison["best_clean_density_method"] = winner["method"]
            comparison["best_clean_density_bits_per_nt"] = winner[
                "logical_bits_per_nucleotide"
            ]
        if dropout:
            best_rate = max(float(row["recovery_rate"]) for row in dropout)
            comparison["best_5pct_dropout_recovery_rate"] = best_rate
            comparison["best_5pct_dropout_methods"] = [
                row["method"] for row in dropout if float(row["recovery_rate"]) == best_rate
            ]

    targets: dict[str, object] = {}
    for target in RESEARCH_SIZES:
        target_rows = [row for row in rows if int(row["size_bytes"]) == target]
        if not target_rows:
            continue
        failures = [
            {
                "method": row.get("method"),
                "condition": row.get("condition"),
                "stage": row.get("stage", "condition"),
                "completed_trials": row.get("completed_trials", 0),
                "timed_out": bool(row.get("timed_out")),
                "error": row.get("error"),
            }
            for row in target_rows
            if "error" in row or not bool(row.get("benchmark_complete"))
        ]
        targets[str(target)] = {
            "common_completed": target in common_sizes,
            "failure_rows": len(failures),
            "failures": failures,
        }

    return {
        "largest_common_completed_size_bytes": largest,
        "common_completed_sizes_bytes": common_sizes,
        "failure_rows": sum("error" in row for row in rows),
        "timed_out_rows": sum(bool(row.get("timed_out")) for row in rows),
        "comparison_at_largest_common_size": comparison,
        "research_target_status": targets,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--profile",
        choices=("ci", "full", "scale", "research", "research-10", "research-100"),
        default="full",
    )
    parser.add_argument("--output", type=Path, default=Path("matched-codec-results"))
    parser.add_argument("--redundancy", type=float, default=0.25)
    parser.add_argument("--trials", type=int)
    parser.add_argument("--timeout-seconds", type=int)
    parser.add_argument("--encode-timeout-seconds", type=int, default=1200)
    parser.add_argument("--condition-timeout-seconds", type=int, default=600)
    parser.add_argument("--keep-encoded-artifacts", action="store_true")
    parser.add_argument("--encode-worker", action="store_true")
    parser.add_argument("--condition-worker", action="store_true")
    parser.add_argument("--method", choices=METHODS)
    parser.add_argument("--size", type=int)
    parser.add_argument("--artifact", type=Path)
    parser.add_argument("--condition", choices=tuple(CONDITIONS))
    parser.add_argument("--checkpoint", type=Path)
    args = parser.parse_args()

    if args.timeout_seconds is not None:
        args.encode_timeout_seconds = args.timeout_seconds
        args.condition_timeout_seconds = args.timeout_seconds

    if args.profile == "ci":
        conditions = ("clean", "dropout-5", "indel-low", "mixed")
        sizes = (1 * KIB,)
        default_trials = 3
    elif args.profile == "full":
        conditions = tuple(CONDITIONS)
        sizes = FULL_SIZES
        default_trials = 10
    elif args.profile == "research":
        conditions = tuple(CONDITIONS)
        sizes = RESEARCH_SIZES
        default_trials = 10
    elif args.profile == "research-10":
        conditions = tuple(CONDITIONS)
        sizes = (10 * MIB,)
        default_trials = 10
    elif args.profile == "research-100":
        conditions = tuple(CONDITIONS)
        sizes = (100 * MIB,)
        default_trials = 10
    else:
        conditions = tuple(CONDITIONS)
        sizes = SCALE_SIZES
        default_trials = 10
    trials = args.trials if args.trials is not None else default_trials

    if args.encode_worker:
        if args.method is None or args.size is None or args.artifact is None:
            raise SystemExit("encode worker requires --method, --size and --artifact")
        print(
            json.dumps(
                _encode_artifact(
                    args.method,
                    args.size,
                    args.redundancy,
                    args.artifact,
                )
            )
        )
        return

    if args.condition_worker:
        if args.artifact is None or args.condition is None or args.checkpoint is None:
            raise SystemExit(
                "condition worker requires --artifact, --condition and --checkpoint"
            )
        print(
            json.dumps(
                _condition_worker(
                    args.artifact,
                    args.condition,
                    trials,
                    args.checkpoint,
                )
            )
        )
        return

    args.output.mkdir(parents=True, exist_ok=True)
    artifact_dir = args.output / "_encoded"
    checkpoint_dir = args.output / "_checkpoints"
    artifact_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    encode_records: list[dict[str, object]] = []
    for size in sizes:
        for method in METHODS:
            artifact_path = artifact_dir / f"{method}-{size}.pickle"
            encode_meta, encode_error, encode_timed_out = _prepare_artifact(
                method,
                size,
                args.redundancy,
                artifact_path,
                args.encode_timeout_seconds,
            )
            if encode_meta is None:
                rows.extend(
                    _encode_failure_rows(
                        method,
                        size,
                        args.redundancy,
                        conditions,
                        trials,
                        encode_error or "unknown encode failure",
                        encode_timed_out,
                        args.encode_timeout_seconds,
                    )
                )
                continue
            encode_records.append(encode_meta)

            for condition in conditions:
                checkpoint = checkpoint_dir / f"{method}-{size}-{condition}.json"
                rows.append(
                    _run_condition_isolated(
                        artifact_path=artifact_path,
                        encode_meta=encode_meta,
                        condition_name=condition,
                        trials=trials,
                        timeout_seconds=args.condition_timeout_seconds,
                        checkpoint_path=checkpoint,
                    )
                )

            if not args.keep_encoded_artifacts:
                os.chmod(artifact_path, 0o644)
                artifact_path.unlink(missing_ok=True)

    metadata: dict[str, Any] = {
        "oligoark_version": __version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "profile": args.profile,
        "payload_sizes_bytes": list(sizes),
        "redundancy_budget": args.redundancy,
        "target_max_strand_nt": 152,
        "trials_per_condition": trials,
        "conditions": list(conditions),
        "methods": list(METHODS),
        "encode_timeout_seconds": args.encode_timeout_seconds,
        "condition_timeout_seconds": args.condition_timeout_seconds,
        "encoded_artifacts_kept": args.keep_encoded_artifacts,
        "orchestration": (
            "one immutable pickle artifact per method x size; every condition loads the same "
            "artifact in an isolated worker; condition checkpoints preserve completed trials "
            "when a worker times out"
        ),
        "fairness": (
            "same deterministic payload bytes, 152-nt ceiling, nominal 25% redundancy budget, "
            "channel rates, trial seeds and SHA-256 exact-recovery definition; both OligoArk "
            "profiles split the budget between XOR parity and fountain symbols"
        ),
        "claim_scope": (
            "software codec comparison only; no wet-lab performance or historical "
            "bit-compatibility is claimed. 100 MiB is a target size: explicit timeout, "
            "resource failure or incomplete trials remain reportable negative results."
        ),
        "encode_records": encode_records,
    }
    summary = _summary(rows, conditions)
    (args.output / "results.json").write_text(
        json.dumps(rows, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (args.output / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _write_csv(rows, args.output / "results.csv")
    _write_plots(rows, args.output)
    print(json.dumps({"metadata": metadata, "summary": summary}, indent=2))


if __name__ == "__main__":
    main()

"""Efficient matched codec scale benchmark.

Each method/size pair is encoded once, then the same encoded archive is evaluated under every
channel condition. This keeps large comparisons practical while preserving identical payload,
strand-length, redundancy, seed and SHA-256 success constraints.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import resource
import signal
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
    DnaFountainBaselineConfig,
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
    "oligoark-compact-v3-hybrid",
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


def _repository_commit() -> str:
    root = Path(__file__).resolve().parents[1]
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.SubprocessError):
        return os.environ.get("GITHUB_SHA", "unknown")


def _metrics(
    *,
    payload_size: int,
    encoded_nucleotides: int,
    encode_seconds: float,
    decode_seconds: float,
    trials: int,
) -> dict[str, object]:
    payload_mib = payload_size / MIB
    mean_decode = decode_seconds / trials
    return {
        "nucleotide_overhead_vs_2bit_ideal": round(
            encoded_nucleotides / max(1, payload_size * 4),
            6,
        ),
        "encode_seconds": round(encode_seconds, 6),
        "mean_decode_seconds": round(mean_decode, 6),
        "encode_throughput_mib_s": round(
            payload_mib / encode_seconds if encode_seconds else 0.0,
            6,
        ),
        "decode_throughput_mib_s": round(
            payload_mib / mean_decode if mean_decode else 0.0,
            6,
        ),
        "peak_rss_mib": round(_peak_rss_mib(), 3),
    }


class _TrialDeadline(TimeoutError):
    """Raised when one matched benchmark trial exceeds its measured budget."""


def _alarm_handler(signum: int, frame: object) -> None:
    del signum, frame
    raise _TrialDeadline


def _condition_row(
    *,
    method: str,
    payload: bytes,
    sequences: tuple[str, ...] | list[str],
    decode: Any,
    condition_name: str,
    trials: int,
    redundancy: float,
    target_max_strand_nt: int,
    data_units: int,
    encoded_nucleotides: int,
    logical_bits_per_nucleotide: float,
    encode_seconds: float,
    condition_timeout_seconds: int,
    extra: dict[str, object] | None = None,
) -> dict[str, object]:
    """Run all requested trials; a timeout is a retained failed trial, never an omitted run."""
    successes = 0
    completed_trials = 0
    timeout_trials = 0
    decode_seconds = 0.0
    channel_seconds = 0.0
    trial_results: list[dict[str, object]] = []
    condition_started = time.perf_counter()
    previous_handler: Any = None
    deadline_enabled = (
        condition_timeout_seconds > 0
        and hasattr(signal, "SIGALRM")
        and hasattr(signal, "setitimer")
    )
    if deadline_enabled:
        previous_handler = signal.getsignal(signal.SIGALRM)
        signal.signal(signal.SIGALRM, _alarm_handler)

    try:
        for trial in range(trials):
            trial_seed = 20_260_000 + trial
            current_stage = "channel"
            trial_started = time.perf_counter()
            decode_started: float | None = None
            if deadline_enabled:
                signal.setitimer(signal.ITIMER_REAL, condition_timeout_seconds)
            try:
                channel_started = time.perf_counter()
                reads = _channel_sequences(sequences, condition_name, trial_seed)
                trial_channel_seconds = time.perf_counter() - channel_started
                channel_seconds += trial_channel_seconds

                current_stage = "decode"
                decode_started = time.perf_counter()
                error: str | None = None
                success = False
                try:
                    recovered = decode(reads)
                    success = recovered == payload
                    successes += int(success)
                except ValueError as exc:
                    error = str(exc)
                trial_decode_seconds = time.perf_counter() - decode_started
                decode_seconds += trial_decode_seconds
                completed_trials += 1
                trial_results.append(
                    {
                        "trial": trial,
                        "seed": trial_seed,
                        "success": success,
                        "sha256_verified": success,
                        "channel_seconds": round(trial_channel_seconds, 6),
                        "decode_seconds": round(trial_decode_seconds, 6),
                        "peak_rss_mib": round(_peak_rss_mib(), 3),
                        "error": error,
                    }
                )
            except _TrialDeadline:
                timeout_trials += 1
                elapsed = time.perf_counter() - trial_started
                if current_stage == "decode" and decode_started is not None:
                    observed_decode = time.perf_counter() - decode_started
                else:
                    observed_decode = 0.0
                trial_results.append(
                    {
                        "trial": trial,
                        "seed": trial_seed,
                        "success": False,
                        "sha256_verified": False,
                        "timed_out": True,
                        "timeout_stage": current_stage,
                        "timeout_seconds": condition_timeout_seconds,
                        "elapsed_trial_seconds": round(elapsed, 6),
                        "observed_decode_seconds_before_timeout": round(
                            observed_decode,
                            6,
                        ),
                        "peak_rss_mib": round(_peak_rss_mib(), 3),
                    }
                )
            finally:
                if deadline_enabled:
                    signal.setitimer(signal.ITIMER_REAL, 0.0)
    finally:
        if deadline_enabled:
            signal.signal(signal.SIGALRM, previous_handler)

    recovery_rate = round(successes / trials, 6)
    wilson_low, wilson_high = wilson_interval(successes, trials)
    strand_count = len(sequences)
    row: dict[str, object] = {
        "method": method,
        "size_bytes": len(payload),
        "target_max_strand_nt": target_max_strand_nt,
        "max_strand_nt": max(map(len, sequences)),
        "redundancy_budget": redundancy,
        "condition": condition_name,
        **CONDITIONS[condition_name],
        "trials": trials,
        "completed_trials": completed_trials,
        "successes": successes,
        "failed_trials": trials - successes,
        "timeout_trials": timeout_trials,
        "recovery_rate": recovery_rate,
        "recovery_ci95_low": round(wilson_low, 6),
        "recovery_ci95_high": round(wilson_high, 6),
        "sha256_verified_success_definition": True,
        "strand_count": strand_count,
        "data_units": data_units,
        "measured_strand_redundancy_ratio": round(
            (strand_count - data_units) / max(1, data_units),
            6,
        ),
        "encoded_nucleotides": encoded_nucleotides,
        "logical_bits_per_nucleotide": round(logical_bits_per_nucleotide, 6),
        "condition_elapsed_seconds": round(
            time.perf_counter() - condition_started,
            6,
        ),
        "channel_seconds": round(channel_seconds, 6),
        "mean_channel_seconds": round(
            channel_seconds / max(1, completed_trials),
            6,
        ),
        "trial_results": trial_results,
        **_metrics(
            payload_size=len(payload),
            encoded_nucleotides=encoded_nucleotides,
            encode_seconds=encode_seconds,
            decode_seconds=decode_seconds,
            trials=max(1, completed_trials),
        ),
    }
    if timeout_trials:
        timeout_stages = sorted(
            {
                str(trial["timeout_stage"])
                for trial in trial_results
                if trial.get("timed_out")
            }
        )
        row.update(
            {
                "timed_out": True,
                "timeout_seconds": condition_timeout_seconds,
                "timeout_scope": "per-trial",
                "timeout_stage": ",".join(timeout_stages),
                "error": (
                    f"{timeout_trials}/{trials} trials exceeded the "
                    f"{condition_timeout_seconds}s per-trial deadline"
                ),
            }
        )
    if extra:
        row.update(extra)
    return row

def _worker(
    method: str,
    size: int,
    redundancy: float,
    trials: int,
    conditions: tuple[str, ...],
    condition_timeout_seconds: int,
) -> list[dict[str, object]]:
    payload = _payload(size)
    payload_sha256 = hashlib.sha256(payload).hexdigest()
    started = time.perf_counter()

    if method in {"oligoark-compact-hybrid", "oligoark-compact-v3-hybrid"}:
        profile_name = (
            "oligoark-152-compact-v3"
            if method == "oligoark-compact-v3-hybrid"
            else "oligoark-152-compact"
        )
        profile = physical_strand_profile(profile_name).with_scheme("hybrid")
        xor_share = 1.0 / profile.parity_group_size
        fountain_share = max(0.0, redundancy - xor_share)
        config = profile.to_archive_config(fountain_redundancy=fountain_share)
        archive = archive_bytes(payload, config)
        encode_seconds = time.perf_counter() - started
        stats = archive_statistics(archive)
        rows = [
            _condition_row(
                method=method,
                payload=payload,
                sequences=archive.strands,
                decode=lambda reads: recover_bytes(archive, reads),
                condition_name=condition,
                trials=trials,
                redundancy=redundancy,
                target_max_strand_nt=152,
                data_units=int(archive.metadata["data_strands"]),
                encoded_nucleotides=stats.encoded_nucleotides,
                logical_bits_per_nucleotide=stats.logical_bits_per_nucleotide,
                encode_seconds=encode_seconds,
                condition_timeout_seconds=condition_timeout_seconds,
                extra={
                    "redundancy_allocation": {
                        "xor": round(xor_share, 6),
                        "fountain": round(fountain_share, 6),
                    },
                    "profile_name": profile_name,
                    "compact_framing": config.compact_framing,
                    "inline_mask_framing": config.inline_mask_framing,
                    "rs_nsym": config.rs_nsym,
                    "chunk_size": config.chunk_size,
                },
            )
            for condition in conditions
        ]
    elif method == "dna-fountain-cleanroom":
        config = DnaFountainBaselineConfig(redundancy=redundancy)
        archive = encode_dna_fountain_baseline(payload, config)
        encode_seconds = time.perf_counter() - started
        rows = [
            _condition_row(
                method=method,
                payload=payload,
                sequences=archive.sequences,
                decode=lambda reads: decode_dna_fountain_baseline(archive, reads),
                condition_name=condition,
                trials=trials,
                redundancy=redundancy,
                target_max_strand_nt=152,
                data_units=archive.chunk_count,
                encoded_nucleotides=archive.encoded_nucleotides,
                logical_bits_per_nucleotide=archive.logical_bits_per_nucleotide,
                encode_seconds=encode_seconds,
                condition_timeout_seconds=condition_timeout_seconds,
                extra={"accepted_droplet_attempts": archive.attempts},
            )
            for condition in conditions
        ]
    elif method == "goldman-rotating-xor":
        config = RotatingTernaryBaselineConfig(redundancy=redundancy)
        archive = encode_rotating_ternary_baseline(payload, config)
        encode_seconds = time.perf_counter() - started
        rows = [
            _condition_row(
                method=method,
                payload=payload,
                sequences=archive.sequences,
                decode=lambda reads: decode_rotating_ternary_baseline(archive, reads),
                condition_name=condition,
                trials=trials,
                redundancy=redundancy,
                target_max_strand_nt=config.max_strand_nt,
                data_units=archive.data_count,
                encoded_nucleotides=archive.encoded_nucleotides,
                logical_bits_per_nucleotide=archive.logical_bits_per_nucleotide,
                encode_seconds=encode_seconds,
                condition_timeout_seconds=condition_timeout_seconds,
            )
            for condition in conditions
        ]
    else:
        raise ValueError(f"unknown method: {method}")

    for row in rows:
        row["payload_sha256"] = payload_sha256
    return rows


def _isolated(
    method: str,
    size: int,
    redundancy: float,
    trials: int,
    conditions: tuple[str, ...],
    timeout_seconds: int,
    condition_timeout_seconds: int,
) -> list[dict[str, object]]:
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--method",
        method,
        "--size",
        str(size),
        "--redundancy",
        str(redundancy),
        "--trials",
        str(trials),
        "--conditions",
        ",".join(conditions),
        "--condition-timeout-seconds",
        str(condition_timeout_seconds),
    ]
    worker_timeout_seconds = (
        timeout_seconds + condition_timeout_seconds * trials * len(conditions)
    )
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=worker_timeout_seconds,
        )
    except subprocess.TimeoutExpired:
        return [
            {
                "method": method,
                "size_bytes": size,
                "redundancy_budget": redundancy,
                "condition": condition,
                "trials": trials,
                "successes": 0,
                "recovery_rate": 0.0,
                "timed_out": True,
                "timeout_seconds": worker_timeout_seconds,
                "timeout_stage": "encode-or-worker-orchestration",
                "error": (
                    f"method/size worker exceeded derived {worker_timeout_seconds}s budget "
                    f"({timeout_seconds}s encode/orchestration + "
                    f"{condition_timeout_seconds}s x {trials} trials x "
                    f"{len(conditions)} conditions)"
                ),
            }
            for condition in conditions
        ]
    if completed.returncode:
        # An encoding failure prevents every downstream recovery attempt. Preserve all
        # requested trials as explicit negative outcomes instead of dropping them from
        # the evidence artifact. This is NOT a successful codec validation.
        error = completed.stderr.strip() or completed.stdout.strip()
        return [
            {
                "method": method,
                "size_bytes": size,
                "redundancy_budget": redundancy,
                "condition": condition,
                "trials": trials,
                "completed_trials": 0,
                "failed_trials": trials,
                "successes": 0,
                "recovery_rate": 0.0,
                "sha256_verified_success_definition": True,
                "encoding_failed": True,
                "error": error,
                "trial_results": [
                    {
                        "trial": trial,
                        "seed": 20_260_000 + trial,
                        "success": False,
                        "sha256_verified": False,
                        "not_executed": True,
                        "failure_stage": "encode-or-worker",
                        "error": "encoding failed; recovery trial could not run",
                    }
                    for trial in range(trials)
                ],
            }
            for condition in conditions
        ]
    value = json.loads(completed.stdout)
    if not isinstance(value, list):
        raise ValueError("matched codec worker did not return a row list")
    return value


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


def _trial_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    flattened: list[dict[str, object]] = []
    for row in rows:
        trials = row.get("trial_results")
        if not isinstance(trials, list):
            continue
        identity = {
            "method": row.get("method"),
            "size_bytes": row.get("size_bytes"),
            "condition": row.get("condition"),
            "payload_sha256": row.get("payload_sha256"),
            "target_max_strand_nt": row.get("target_max_strand_nt"),
            "redundancy_budget": row.get("redundancy_budget"),
        }
        for trial in trials:
            if isinstance(trial, dict):
                flattened.append({**identity, **trial})
    return flattened


def _summary(
    rows: list[dict[str, object]],
    conditions: tuple[str, ...],
    methods: tuple[str, ...],
) -> dict[str, object]:
    sizes = sorted({int(row["size_bytes"]) for row in rows})
    common_sizes: list[int] = []
    for size in sizes:
        complete = True
        for method in methods:
            selected = [
                row
                for row in rows
                if row.get("method") == method and int(row["size_bytes"]) == size
            ]
            if (
                len(selected) != len(conditions)
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
            and "recovery_rate" in row
            and row.get("recovery_rate") is not None
            and "error" not in row
        ]
        if clean:
            density_winner = max(clean, key=lambda row: float(row["logical_bits_per_nucleotide"]))
            comparison["best_clean_density_method"] = density_winner["method"]
            comparison["best_clean_density_bits_per_nt"] = density_winner[
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
                "timed_out": bool(row.get("timed_out")),
                "error": row.get("error"),
            }
            for row in target_rows
            if "error" in row
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
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=1200,
        help="encode/orchestration watchdog budget before per-condition budgets are added",
    )
    parser.add_argument(
        "--condition-timeout-seconds",
        type=int,
        default=300,
        help="hard per-trial deadline; retained timeouts count as failed trials",
    )
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--method", choices=METHODS)
    parser.add_argument("--size", type=int)
    parser.add_argument("--conditions", type=str)
    args = parser.parse_args()

    if args.profile == "ci":
        default_conditions = ("clean", "dropout-5", "indel-low", "mixed")
        sizes = (1 * KIB,)
        default_trials = 3
    elif args.profile == "full":
        default_conditions = tuple(CONDITIONS)
        sizes = FULL_SIZES
        default_trials = 10
    elif args.profile == "research":
        default_conditions = tuple(CONDITIONS)
        sizes = RESEARCH_SIZES
        default_trials = 10
    elif args.profile == "research-10":
        default_conditions = tuple(CONDITIONS)
        sizes = (10 * MIB,)
        default_trials = 10
    elif args.profile == "research-100":
        default_conditions = tuple(CONDITIONS)
        sizes = (100 * MIB,)
        default_trials = 10
    else:
        default_conditions = tuple(CONDITIONS)
        sizes = SCALE_SIZES
        default_trials = 10

    trials = args.trials if args.trials is not None else default_trials

    if args.worker:
        if args.method is None or args.size is None:
            raise SystemExit("worker requires --method and --size")
        conditions = (
            tuple(args.conditions.split(","))
            if args.conditions
            else default_conditions
        )
        print(
            json.dumps(
                _worker(
                    args.method,
                    args.size,
                    args.redundancy,
                    trials,
                    conditions,
                    args.condition_timeout_seconds,
                )
            )
        )
        return

    if args.size is not None:
        sizes = (args.size,)
    selected_methods = (args.method,) if args.method is not None else METHODS
    conditions = default_conditions
    args.output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for size in sizes:
        for method in selected_methods:
            rows.extend(
                _isolated(
                    method,
                    size,
                    args.redundancy,
                    trials,
                    conditions,
                    args.timeout_seconds,
                    args.condition_timeout_seconds,
                )
            )

    metadata: dict[str, Any] = {
        "oligoark_version": __version__,
        "repository_commit": _repository_commit(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "profile": args.profile,
        "payload_sizes_bytes": list(sizes),
        "redundancy_budget": args.redundancy,
        "target_max_strand_nt": 152,
        "trials_per_condition": trials,
        "conditions": list(conditions),
        "methods": list(selected_methods),
        "trial_seeds": [20_260_000 + trial for trial in range(trials)],
        "argv": [sys.executable, *sys.argv],
        "channel_engine": "sparse-geometric-v1",
        "encode_orchestration_timeout_seconds": args.timeout_seconds,
        "trial_timeout_seconds": args.condition_timeout_seconds,
        "timeout_scope": "per-trial",
        "derived_worker_timeout_seconds": (
            args.timeout_seconds
            + args.condition_timeout_seconds * trials * len(conditions)
        ),
        "fairness": (
            "same deterministic payload bytes, 152-nt ceiling, nominal 25% redundancy budget, "
            "channel rates, trial seeds and SHA-256 exact-recovery definition; clean trials "
            "reuse immutable encoded artifacts, dropout-only trials skip mutation work, and "
            "noisy channels use sparse event sampling; OligoArk splits "
            "the budget between XOR parity and fountain symbols"
        ),
        "claim_scope": (
            "software codec comparison only; no wet-lab performance or historical "
            "bit-compatibility is claimed. 100 MiB is a target size: explicit timeout or "
            "resource failure remains a reportable negative result."
        ),
    }
    summary = _summary(rows, conditions, selected_methods)
    trial_rows = _trial_rows(rows)
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
    (args.output / "trial-results.json").write_text(
        json.dumps(trial_rows, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _write_csv(rows, args.output / "results.csv")
    _write_csv(trial_rows, args.output / "trial-results.csv")
    _write_plots(rows, args.output)
    print(json.dumps({"metadata": metadata, "summary": summary}, indent=2))


if __name__ == "__main__":
    main()

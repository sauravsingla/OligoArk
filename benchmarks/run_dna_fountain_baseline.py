"""Fair software comparison between OligoArk and a clean-room DNA Fountain baseline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import random
import resource
import subprocess
import sys
import time
from pathlib import Path

from oligoark import __version__
from oligoark.archive import archive_bytes, archive_statistics, recover_bytes
from oligoark.baselines import (
    DnaFountainBaselineConfig,
    decode_dna_fountain_baseline,
    encode_dna_fountain_baseline,
)
from oligoark.experiments import wilson_interval
from oligoark.profiles import physical_strand_profile

METHODS = ("oligoark-fountain", "dna-fountain-cleanroom")


def _payload(size: int) -> bytes:
    rng = random.Random(2026)
    text = (
        b"OligoArk/DNA-Fountain controlled comparison\n"
        b"text,json,source,image-like,binary,compressed\n"
    ) * 64
    structured = b"".join(
        json.dumps({"id": i, "value": round(i / 17, 6)}, sort_keys=True).encode() + b"\n"
        for i in range(128)
    )
    source = (
        b"def encode(value: bytes) -> bytes:\n"
        b"    return value\n"
    ) * 128
    binary = bytes(rng.randrange(256) for _ in range(64 * 1024))
    block = text + structured + source + binary
    repeats, remainder = divmod(size, len(block))
    return block * repeats + block[:remainder]


def _peak_rss_mib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return value / (1024 * 1024)
    return value / 1024


def _drop_sequences(
    sequences: tuple[str, ...] | list[str],
    dropout_rate: float,
    seed: int,
) -> list[str]:
    selected = list(sequences)
    count = min(len(selected), round(len(selected) * dropout_rate))
    rng = random.Random(seed)
    dropped = set(rng.sample(range(len(selected)), count)) if count else set()
    return [sequence for index, sequence in enumerate(selected) if index not in dropped]


def _oligoark(
    payload: bytes,
    redundancy: float,
    dropout_rate: float,
    trials: int,
) -> dict[str, object]:
    profile = physical_strand_profile("oligoark-152").with_scheme("fountain")
    config = profile.to_archive_config(fountain_redundancy=redundancy)
    started = time.perf_counter()
    archive = archive_bytes(payload, config)
    encode_seconds = time.perf_counter() - started
    stats = archive_statistics(archive)

    successes = 0
    decode_seconds = 0.0
    for trial in range(trials):
        reads = _drop_sequences(archive.strands, dropout_rate, 20_260_000 + trial)
        trial_started = time.perf_counter()
        try:
            recovered = recover_bytes(archive, reads)
            successes += int(recovered == payload)
        except ValueError:
            pass
        decode_seconds += time.perf_counter() - trial_started

    low, high = wilson_interval(successes, trials)
    return {
        "method": "oligoark-fountain",
        "size_bytes": len(payload),
        "target_max_strand_nt": 152,
        "max_strand_nt": max(map(len, archive.strands)),
        "redundancy_budget": redundancy,
        "dropout_rate": dropout_rate,
        "trials": trials,
        "successes": successes,
        "recovery_rate": round(successes / trials, 6),
        "recovery_ci95_low": round(low, 6),
        "recovery_ci95_high": round(high, 6),
        "sha256_verified_success_definition": True,
        "strand_count": len(archive.strands),
        "encoded_nucleotides": stats.encoded_nucleotides,
        "logical_bits_per_nucleotide": stats.logical_bits_per_nucleotide,
        "encode_seconds": round(encode_seconds, 6),
        "mean_decode_seconds": round(decode_seconds / trials, 6),
        "peak_rss_mib": round(_peak_rss_mib(), 3),
    }


def _dna_fountain(
    payload: bytes,
    redundancy: float,
    dropout_rate: float,
    trials: int,
) -> dict[str, object]:
    config = DnaFountainBaselineConfig(redundancy=redundancy)
    started = time.perf_counter()
    archive = encode_dna_fountain_baseline(payload, config)
    encode_seconds = time.perf_counter() - started

    successes = 0
    decode_seconds = 0.0
    for trial in range(trials):
        reads = _drop_sequences(archive.sequences, dropout_rate, 20_260_000 + trial)
        trial_started = time.perf_counter()
        try:
            recovered = decode_dna_fountain_baseline(archive, reads)
            successes += int(recovered == payload)
        except ValueError:
            pass
        decode_seconds += time.perf_counter() - trial_started

    low, high = wilson_interval(successes, trials)
    return {
        "method": "dna-fountain-cleanroom",
        "size_bytes": len(payload),
        "target_max_strand_nt": 152,
        "max_strand_nt": max(map(len, archive.sequences)),
        "redundancy_budget": redundancy,
        "dropout_rate": dropout_rate,
        "trials": trials,
        "successes": successes,
        "recovery_rate": round(successes / trials, 6),
        "recovery_ci95_low": round(low, 6),
        "recovery_ci95_high": round(high, 6),
        "sha256_verified_success_definition": True,
        "strand_count": len(archive.sequences),
        "encoded_nucleotides": archive.encoded_nucleotides,
        "logical_bits_per_nucleotide": round(archive.logical_bits_per_nucleotide, 6),
        "encode_seconds": round(encode_seconds, 6),
        "mean_decode_seconds": round(decode_seconds / trials, 6),
        "peak_rss_mib": round(_peak_rss_mib(), 3),
        "accepted_droplet_attempts": archive.attempts,
    }


def _worker(
    method: str,
    size: int,
    redundancy: float,
    dropout_rate: float,
    trials: int,
) -> dict[str, object]:
    payload = _payload(size)
    if method == "oligoark-fountain":
        row = _oligoark(payload, redundancy, dropout_rate, trials)
    elif method == "dna-fountain-cleanroom":
        row = _dna_fountain(payload, redundancy, dropout_rate, trials)
    else:
        raise ValueError(f"unknown method: {method}")
    row["payload_sha256"] = hashlib.sha256(payload).hexdigest()
    return row


def _isolated(
    method: str,
    size: int,
    redundancy: float,
    dropout_rate: float,
    trials: int,
) -> dict[str, object]:
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
        "--dropout",
        str(dropout_rate),
        "--trials",
        str(trials),
    ]
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode:
        return {
            "method": method,
            "size_bytes": size,
            "redundancy_budget": redundancy,
            "dropout_rate": dropout_rate,
            "trials": trials,
            "successes": 0,
            "recovery_rate": 0.0,
            "error": completed.stderr.strip() or completed.stdout.strip(),
        }
    return json.loads(completed.stdout)


def _write_csv(rows: list[dict[str, object]], path: Path) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _write_plot(rows: list[dict[str, object]], output: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    for method in METHODS:
        selected = sorted(
            (row for row in rows if row["method"] == method),
            key=lambda row: float(row["dropout_rate"]),
        )
        plt.plot(
            [float(row["dropout_rate"]) for row in selected],
            [float(row["recovery_rate"]) for row in selected],
            marker="o",
            label=method,
        )
    plt.xlabel("Strand dropout rate")
    plt.ylabel("SHA-256 verified recovery rate")
    plt.ylim(-0.02, 1.02)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output / "dna_fountain_recovery_comparison.png", dpi=160)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("ci", "full"), default="ci")
    parser.add_argument("--output", type=Path, default=Path("dna-fountain-results"))
    parser.add_argument("--size", type=int)
    parser.add_argument("--redundancy", type=float, default=0.25)
    parser.add_argument("--trials", type=int)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--method", choices=METHODS)
    parser.add_argument("--dropout", type=float)
    args = parser.parse_args()

    if args.worker:
        if args.method is None or args.size is None or args.dropout is None or args.trials is None:
            raise SystemExit("worker requires method, size, dropout and trials")
        print(
            json.dumps(
                _worker(
                    args.method,
                    args.size,
                    args.redundancy,
                    args.dropout,
                    args.trials,
                )
            )
        )
        return

    size = args.size if args.size is not None else (1024 if args.profile == "ci" else 8192)
    trials = args.trials if args.trials is not None else (3 if args.profile == "ci" else 20)
    dropouts = (0.0, 0.05) if args.profile == "ci" else (0.0, 0.01, 0.05, 0.10)
    args.output.mkdir(parents=True, exist_ok=True)
    rows = [
        _isolated(method, size, args.redundancy, dropout, trials)
        for method in METHODS
        for dropout in dropouts
    ]
    metadata = {
        "oligoark_version": __version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "payload_size_bytes": size,
        "redundancy_budget": args.redundancy,
        "target_max_strand_nt": 152,
        "trials_per_condition": trials,
        "claim_scope": (
            "software codec comparison only; the DNA Fountain baseline is an independent "
            "research implementation and is not claimed bit-compatible with TeamErlich"
        ),
        "fairness": (
            "same payload, 152-nt ceiling, nominal redundancy budget, dropout rates, "
            "trial seeds and SHA-256 exact-recovery gate"
        ),
    }
    (args.output / "results.json").write_text(
        json.dumps(rows, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (args.output / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _write_csv(rows, args.output / "results.csv")
    _write_plot(rows, args.output)
    print(json.dumps({"metadata": metadata, "results": rows}, indent=2))


if __name__ == "__main__":
    main()

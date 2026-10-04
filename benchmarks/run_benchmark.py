"""Reproducible software benchmark. Results are simulation outputs, not wet-lab claims."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import random
import time
from pathlib import Path

from oligoark import __version__
from oligoark.archive import ArchiveConfig, archive_bytes, recover_from_reads
from oligoark.policy import ChannelProfile, recommend_codec_policy
from oligoark.simulator import SimulationConfig, simulate_channel

SEED = 2026
PAYLOAD_BYTES = 8192


def run() -> list[dict[str, object]]:
    rng = random.Random(SEED)
    payload = bytes(rng.randrange(256) for _ in range(PAYLOAD_BYTES))
    regimes = [
        ("clean", ChannelProfile(), 0.0),
        ("substitution-low", ChannelProfile(substitution_rate=0.001), 0.10),
        ("substitution-moderate", ChannelProfile(substitution_rate=0.01), 0.25),
        ("insertion-low", ChannelProfile(insertion_rate=0.001), 0.25),
        ("deletion-low", ChannelProfile(deletion_rate=0.001), 0.25),
        ("dropout-low", ChannelProfile(dropout_rate=0.02), 0.10),
        (
            "mixed",
            ChannelProfile(
                substitution_rate=0.002,
                insertion_rate=0.0005,
                deletion_rate=0.0005,
                dropout_rate=0.02,
            ),
            0.25,
        ),
    ]
    rows: list[dict[str, object]] = []
    for name, channel, duplicate_rate in regimes:
        adaptive = recommend_codec_policy(channel)
        policies = {
            "fixed": ArchiveConfig(96, 8, 8, False),
            "adaptive": ArchiveConfig(
                adaptive.chunk_size,
                adaptive.rs_nsym,
                adaptive.parity_group_size,
                adaptive.adaptive_masks,
            ),
        }
        for policy_name, config in policies.items():
            started = time.perf_counter()
            archive = archive_bytes(payload, config)
            reads = simulate_channel(
                archive.strands,
                SimulationConfig(
                    substitution_rate=channel.substitution_rate,
                    insertion_rate=channel.insertion_rate,
                    deletion_rate=channel.deletion_rate,
                    dropout_rate=channel.dropout_rate,
                    duplicate_rate=duplicate_rate,
                    seed=SEED,
                ),
            )
            recovered = False
            graph_used = False
            try:
                recovered_payload, report = recover_from_reads(archive, reads)
                recovered = recovered_payload == payload
                graph_used = report.graph_reconstruction_used
            except ValueError:
                recovered = False
            measured = archive.metadata["measured"]
            if not isinstance(measured, dict):
                raise ValueError("archive measured metadata must be an object")
            encoded_nucleotides = int(measured["encoded_nucleotides"])
            ideal_nucleotides = max(1, len(payload) * 4)
            rows.append(
                {
                    "regime": name,
                    "policy": policy_name,
                    "recovered": recovered,
                    "graph_reconstruction_used": graph_used,
                    "encoded_nucleotides": encoded_nucleotides,
                    "redundancy_ratio_vs_2bit_ideal": round(
                        encoded_nucleotides / ideal_nucleotides,
                        6,
                    ),
                    "strand_count": len(archive.strands),
                    "read_count": len(reads),
                    "runtime_seconds": round(time.perf_counter() - started, 6),
                }
            )
    return rows


def _write_plot(rows: list[dict[str, object]], output: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    labels = [f"{row['regime']}\n{row['policy']}" for row in rows]
    recovery = [1 if row["recovered"] else 0 for row in rows]
    plt.figure(figsize=(13, 4))
    plt.bar(labels, recovery)
    plt.ylabel("Recovery success (0/1)")
    plt.xticks(rotation=40, ha="right")
    plt.tight_layout()
    plt.savefig(output / "recovery_by_regime.png", dpi=160)
    plt.close()

    overhead = [float(row["redundancy_ratio_vs_2bit_ideal"]) for row in rows]
    plt.figure(figsize=(13, 4))
    plt.bar(labels, overhead)
    plt.ylabel("Encoded nt / ideal 2-bit nt")
    plt.xticks(rotation=40, ha="right")
    plt.tight_layout()
    plt.savefig(output / "overhead_by_regime.png", dpi=160)
    plt.close()


def main() -> None:
    rows = run()
    output = Path("benchmark-results")
    output.mkdir(exist_ok=True)
    (output / "results.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    with (output / "results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    rng = random.Random(SEED)
    payload = bytes(rng.randrange(256) for _ in range(PAYLOAD_BYTES))
    metadata = {
        "oligoark_version": __version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "seed": SEED,
        "payload_bytes": PAYLOAD_BYTES,
        "payload_sha256": hashlib.sha256(payload).hexdigest(),
        "claim_scope": "software simulation only; no wet-lab performance is implied",
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    _write_plot(rows, output)
    print(json.dumps({"metadata": metadata, "results": rows}, indent=2))


if __name__ == "__main__":
    main()

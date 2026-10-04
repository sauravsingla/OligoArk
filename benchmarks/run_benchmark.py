"""Reproducible software benchmark. Results are simulation outputs, not wet-lab claims."""

from __future__ import annotations
import csv
import json
import random
import time
from pathlib import Path
from oligoark.archive import ArchiveConfig, archive_bytes, recover_bytes
from oligoark.policy import ChannelProfile, recommend_codec_policy
from oligoark.simulator import SimulationConfig, simulate_channel


def run() -> list[dict[str, object]]:
    rng = random.Random(2026)
    payload = bytes(rng.randrange(256) for _ in range(8192))
    regimes = [
        ("clean", ChannelProfile()),
        ("substitution-low", ChannelProfile(substitution_rate=0.001)),
        ("substitution-moderate", ChannelProfile(substitution_rate=0.01)),
        ("dropout-low", ChannelProfile(dropout_rate=0.02)),
    ]
    rows: list[dict[str, object]] = []
    for name, channel in regimes:
        adaptive = recommend_codec_policy(channel)
        policies = {
            "fixed": ArchiveConfig(96, 8, 8, False),
            "adaptive": ArchiveConfig(adaptive.chunk_size, adaptive.rs_nsym,
                                      adaptive.parity_group_size, adaptive.adaptive_masks),
        }
        for policy_name, cfg in policies.items():
            started = time.perf_counter()
            archive = archive_bytes(payload, cfg)
            reads = simulate_channel(archive.strands, SimulationConfig(
                substitution_rate=channel.substitution_rate,
                insertion_rate=channel.insertion_rate,
                deletion_rate=channel.deletion_rate,
                dropout_rate=channel.dropout_rate, seed=2026))
            recovered = False
            try:
                recovered = recover_bytes(archive, reads) == payload
            except ValueError:
                pass
            rows.append({
                "regime": name, "policy": policy_name, "recovered": recovered,
                "encoded_nucleotides": archive.metadata["measured"]["encoded_nucleotides"],
                "strand_count": len(archive.strands),
                "runtime_seconds": round(time.perf_counter() - started, 6),
            })
    return rows


def main() -> None:
    rows = run()
    out = Path("benchmark-results")
    out.mkdir(exist_ok=True)
    (out / "results.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    with (out / "results.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    try:
        import matplotlib.pyplot as plt
        labels = [f"{row['regime']}\n{row['policy']}" for row in rows]
        values = [1 if row["recovered"] else 0 for row in rows]
        plt.figure(figsize=(10, 4))
        plt.bar(labels, values)
        plt.ylabel("Recovery success (0/1)")
        plt.xticks(rotation=35, ha="right")
        plt.tight_layout()
        plt.savefig(out / "recovery_by_regime.png", dpi=160)
        plt.close()
    except ImportError:
        pass
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()

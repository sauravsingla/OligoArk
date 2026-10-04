"""Verified adaptive optimisation -> corruption -> graph/alignment recovery demonstration."""

from __future__ import annotations

import json
import random
from pathlib import Path

from oligoark.archive import archive_bytes, recover_from_reads
from oligoark.optimizer import CodecSearchSpace, optimize_codec
from oligoark.policy import ChannelProfile
from oligoark.reconstruct import GraphConsensusReconstructor
from oligoark.simulator import SimulationConfig, simulate_channel
from oligoark.tiering import WorkloadProfile

PAYLOAD_SEED = 404
CALIBRATION_SEED = 2026
EVALUATION_SEED = 2027


def main() -> None:
    rng = random.Random(PAYLOAD_SEED)
    payload = bytes(rng.randrange(256) for _ in range(192))
    workload = WorkloadProfile(
        retention_years=100,
        accesses_per_year=0.1,
        mutability=0.0,
        retrieval_urgency=0.2,
        durability_priority=1.0,
        energy_priority=0.5,
        redundancy_priority=0.9,
        cost_priority=0.4,
        data_size_gb=len(payload) / 1_000_000_000,
    )
    scenarios = {
        "substitution": (ChannelProfile(substitution_rate=0.004), 0.30),
        "indel": (
            ChannelProfile(insertion_rate=0.0008, deletion_rate=0.0008),
            0.60,
        ),
        "dropout": (ChannelProfile(dropout_rate=0.04), 0.20),
        "mixed": (
            ChannelProfile(
                substitution_rate=0.002,
                insertion_rate=0.0004,
                deletion_rate=0.0004,
                dropout_rate=0.02,
            ),
            0.50,
        ),
    }
    search = CodecSearchSpace(
        chunk_sizes=(32, 48, 64),
        rs_nsyms=(8, 16, 24),
        redundancy_schemes=("xor", "fountain", "hybrid"),
        parity_group_sizes=(3, 5),
        fountain_redundancies=(0.35, 0.75),
        reconstruction_modes=("graph",),
        max_candidates=12,
        mask_search_limit=128,
    )
    results: list[dict[str, object]] = []

    for name, (channel, duplicate_rate) in scenarios.items():
        optimization = optimize_codec(
            payload,
            channel,
            workload,
            search_space=search,
            seeds=(CALIBRATION_SEED,),
        )
        archive = archive_bytes(payload, optimization.best_config)
        reads = simulate_channel(
            archive.strands,
            SimulationConfig(
                substitution_rate=channel.substitution_rate,
                insertion_rate=channel.insertion_rate,
                deletion_rate=channel.deletion_rate,
                dropout_rate=channel.dropout_rate,
                duplicate_rate=duplicate_rate,
                seed=EVALUATION_SEED,
            ),
        )
        recovered = False
        graph_used = False
        failure: str | None = None
        try:
            decoded, report = recover_from_reads(
                archive,
                reads,
                reconstructor=GraphConsensusReconstructor(
                    threshold=0.88,
                    consensus_mode="alignment",
                ),
            )
            recovered = decoded == payload
            graph_used = report.graph_reconstruction_used
        except ValueError as exc:
            failure = str(exc)
        results.append(
            {
                "scenario": name,
                "sha256_verified_recovery": recovered,
                "graph_reconstruction_used": graph_used,
                "selected_redundancy_scheme": optimization.best_config.redundancy_scheme,
                "selected_rs_nsym": optimization.best_config.rs_nsym,
                "selected_chunk_size": optimization.best_config.chunk_size,
                "selected_gc_range": [
                    optimization.best_config.min_gc_fraction,
                    optimization.best_config.max_gc_fraction,
                ],
                "selected_max_homopolymer": optimization.best_config.max_homopolymer,
                "calibration_seed": CALIBRATION_SEED,
                "evaluation_seed": EVALUATION_SEED,
                "failure": failure,
            }
        )

    output = Path("verified-demo-results.json")
    output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()

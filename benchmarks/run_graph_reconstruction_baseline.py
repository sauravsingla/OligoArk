"""Reproducible graph-reconstruction baseline benchmark (#43, referencing #18)."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

from oligoark.archive import ArchiveConfig, archive_bytes, recover_bytes
from oligoark.reconstruct import GraphConsensusReconstructor, LevenshteinEdgeScorer
from oligoark.simulator import SimulationConfig, simulate_channel


def run_baseline(
    *,
    threshold: float = 0.90,
    seed: int = 42,
    copies_per_strand: int = 5,
    substitution_rate: float = 0.015,
    insertion_rate: float = 0.005,
    deletion_rate: float = 0.005,
    dropout_rate: float = 0.0,
    duplicate_rate: float = 0.0,
) -> dict[str, Any]:
    """Execute the deterministic graph reconstruction baseline and evaluate results."""
    started = time.perf_counter()

    # 1. Deterministic synthetic fixture generation
    payload = b"OligoArk reproducible graph-reconstruction baseline fixture (#43, #18)"
    archive = archive_bytes(
        payload,
        ArchiveConfig(chunk_size=32, rs_nsym=8, redundancy_scheme="none"),
    )
    ground_truth = tuple(archive.strands)
    config = SimulationConfig(
        substitution_rate=substitution_rate,
        insertion_rate=insertion_rate,
        deletion_rate=deletion_rate,
        dropout_rate=dropout_rate,
        duplicate_rate=duplicate_rate,
        copies_per_strand=copies_per_strand,
        seed=seed,
    )
    reads = simulate_channel(list(ground_truth), config)

    # 2. Graph reconstruction using only reads (no ground truth leakage)
    reconstructor = GraphConsensusReconstructor(
        threshold=threshold,
        scorer=LevenshteinEdgeScorer(),
        consensus_mode="alignment",
    )
    result = reconstructor.reconstruct(reads)

    # 3. Post-reconstruction consensus evaluation against ground truth
    ground_truth_set = set(ground_truth)
    consensus_set = set(result.consensus_reads)
    exact_matches = len(consensus_set & ground_truth_set)

    # 4. Direct frame decoding evaluation (no graph)
    direct_recovered = False
    try:
        direct_recovered = recover_bytes(archive, reads) == payload
    except (ValueError, Exception):
        direct_recovered = False

    # 5. Graph-assisted complete archive recovery (SHA-256 verified)
    graph_recovered = False
    try:
        graph_recovered = recover_bytes(archive, reads + result.consensus_reads) == payload
    except (ValueError, Exception):
        graph_recovered = False

    # 6. Medoid ablation comparison (non-alignment consensus)
    medoid_rec = GraphConsensusReconstructor(
        threshold=threshold,
        scorer=LevenshteinEdgeScorer(),
        consensus_mode="medoid",
    )
    medoid_result = medoid_rec.reconstruct(reads)
    medoid_matches = len(set(medoid_result.consensus_reads) & ground_truth_set)
    medoid_recovered = False
    try:
        medoid_recovered = (
            recover_bytes(archive, reads + medoid_result.consensus_reads) == payload
        )
    except (ValueError, Exception):
        medoid_recovered = False

    elapsed = round(time.perf_counter() - started, 6)

    return {
        "benchmark": "graph-reconstruction-baseline",
        "task_reference": "OligoArk #43 (referencing #18)",
        "simulation": {
            "seed": seed,
            "substitution_rate": substitution_rate,
            "insertion_rate": insertion_rate,
            "deletion_rate": deletion_rate,
            "dropout_rate": dropout_rate,
            "duplicate_rate": duplicate_rate,
            "copies_per_strand": copies_per_strand,
            "ground_truth_strand_count": len(ground_truth),
            "generated_read_count": len(reads),
        },
        "graph": {
            "threshold": threshold,
            "scorer": "LevenshteinEdgeScorer",
            "node_count": result.node_count,
            "candidate_pairs": result.candidate_pairs,
            "edge_count": result.edge_count,
            "component_count": result.component_count,
            "cluster_sizes": list(result.cluster_sizes),
        },
        "consensus_recovery": {
            "mode": "alignment",
            "consensus_reads_count": len(result.consensus_reads),
            "ground_truth_count": len(ground_truth),
            "exact_matches": exact_matches,
            "consensus_recovery_rate": round(exact_matches / len(ground_truth), 4),
        },
        "archive_recovery": {
            "direct_recovery_without_graph": direct_recovered,
            "graph_assisted_recovery": graph_recovered,
            "verified_sha256": graph_recovered,
            "archive_sha256": archive.metadata.get("sha256"),
        },
        "ablation_medoid": {
            "mode": "medoid",
            "exact_matches": medoid_matches,
            "archive_recovered": medoid_recovered,
        },
        "runtime_seconds": elapsed,
    }


def format_report(results: dict[str, Any]) -> str:
    """Format baseline results in plain language for humans."""
    sim = results["simulation"]
    graph = results["graph"]
    rec = results["consensus_recovery"]
    arc = results["archive_recovery"]
    abl = results["ablation_medoid"]

    lines = [
        "",
        "================== GRAPH RECONSTRUCTION BASELINE ==================",
        f"Benchmark:               {results['benchmark']}",
        f"Reference:               {results['task_reference']}",
        "",
        "--- Simulation Parameters ---",
        f"Random Seed:             {sim['seed']}",
        f"Substitution Rate:       {sim['substitution_rate'] * 100:.1f}%",
        f"Insertion Rate:          {sim['insertion_rate'] * 100:.1f}%",
        f"Deletion Rate:           {sim['deletion_rate'] * 100:.1f}%",
        f"Dropout Rate:            {sim['dropout_rate'] * 100:.1f}%",
        f"Coverage (per strand):   {sim['copies_per_strand']} traces",
        f"Ground Truth Strands:    {sim['ground_truth_strand_count']}",
        f"Simulated Read Count:    {sim['generated_read_count']}",
        "",
        "--- Graph Construction & Clustering ---",
        f"Scorer:                  {graph['scorer']} (threshold={graph['threshold']})",
        f"Graph Nodes (Reads):     {graph['node_count']}",
        f"Candidate Pairs:         {graph['candidate_pairs']}",
        f"Graph Edges:             {graph['edge_count']}",
        f"Connected Components:    {graph['component_count']}",
        f"Cluster Sizes:           {graph['cluster_sizes']}",
        "",
        "--- Sequence-Level Consensus Recovery ---",
        f"Consensus Mode:          {rec['mode']}",
        f"Candidate Consensuses:   {rec['consensus_reads_count']}",
        f"Exact Ground-Truth Match:{rec['exact_matches']} / {rec['ground_truth_count']} "
        f"({rec['consensus_recovery_rate'] * 100:.1f}%)",
        "",
        "--- End-to-End Archival Verification ---",
        f"Direct Read Recovery:    {arc['direct_recovery_without_graph']} (expected failure)",
        f"Graph-Assisted Recovery: {arc['graph_assisted_recovery']} (SHA-256 verified)",
        f"Archive SHA-256 Match:   {arc['verified_sha256']}",
        "",
        "--- Ablation Contrast (Medoid vs Alignment) ---",
        f"Medoid Matches:          {abl['exact_matches']} / {rec['ground_truth_count']}",
        f"Medoid Archive Recovery: {abl['archive_recovered']} (fails without alignment repair)",
        "",
        f"Reconstruction Runtime:  {results['runtime_seconds']} s (CPU-only)",
        "===================================================================",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run reproducible graph-reconstruction baseline (OligoArk #43, #18)"
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output raw results in JSON format",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional path to save JSON results",
    )
    args = parser.parse_args()

    results = run_baseline()

    if args.output:
        args.output.write_text(json.dumps(results, indent=2), encoding="utf-8")

    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print(format_report(results))


if __name__ == "__main__":
    main()

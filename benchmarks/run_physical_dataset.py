"""Evaluate OligoArk reconstruction baselines on supplied physical FASTA/FASTQ data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from oligoark.physical import (
    PhysicalDatasetManifest,
    evaluate_physical_reconstruction,
    read_sequences,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="datasets/dna_aeon.json")
    parser.add_argument("--reads", required=True)
    parser.add_argument("--references", required=True)
    parser.add_argument("--max-reads", type=int, default=5000)
    parser.add_argument("--assignment-threshold", type=float, default=0.70)
    parser.add_argument("--output", default="physical-reconstruction-results.json")
    args = parser.parse_args()

    manifest = PhysicalDatasetManifest.load(args.manifest)
    reads = read_sequences(args.reads, max_sequences=args.max_reads)
    references = read_sequences(args.references)
    summary = evaluate_physical_reconstruction(
        manifest,
        reads,
        references,
        assignment_threshold=args.assignment_threshold,
    )
    output = Path(args.output)
    output.write_text(json.dumps(summary.to_dict(), indent=2), encoding="utf-8")
    print(json.dumps(summary.to_dict(), indent=2))


if __name__ == "__main__":
    main()

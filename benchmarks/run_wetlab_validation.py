"""Prepare or verify an OligoArk wet-lab archive experiment.

This script does not synthesize or sequence DNA. The prepare command creates the exact oligo
and metadata package to hand to a lab. The recover command consumes preprocessed FASTA/FASTQ
reads and requires exact SHA-256 archive recovery.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from oligoark.wetlab import prepare_wetlab_bundle, recover_wetlab_reads


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("source", type=Path)
    prepare.add_argument("--output", type=Path, default=Path("wetlab-bundle"))
    prepare.add_argument(
        "--profile",
        choices=(
            "oligoark-152-compact",
            "oligoark-152-efficient-v1",
            "oligoark-200-compact",
            "oligoark-248-compact",
        ),
        default="oligoark-200-compact",
    )
    prepare.add_argument(
        "--scheme",
        choices=("none", "xor", "fountain", "hybrid"),
        default="hybrid",
    )
    prepare.add_argument("--fountain-redundancy", type=float, default=0.125)

    recover = subparsers.add_parser("recover")
    recover.add_argument("bundle", type=Path)
    recover.add_argument("reads", type=Path)
    recover.add_argument("output", type=Path)
    recover.add_argument("--similarity-threshold", type=float, default=0.86)
    recover.add_argument("--trace-rounds", type=int, default=4)

    args = parser.parse_args()
    if args.command == "prepare":
        result = prepare_wetlab_bundle(
            args.source,
            args.output,
            profile_name=args.profile,
            redundancy_scheme=args.scheme,
            fountain_redundancy=args.fountain_redundancy,
        )
    else:
        result = recover_wetlab_reads(
            args.bundle,
            args.reads,
            args.output,
            similarity_threshold=args.similarity_threshold,
            trace_rounds=args.trace_rounds,
        )
    print(json.dumps(result.to_dict(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

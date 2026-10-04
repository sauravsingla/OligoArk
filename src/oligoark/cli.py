"""Command-line interface for OligoArk."""

from __future__ import annotations
import argparse
import json
from pathlib import Path
from .archive import ArchiveConfig, DNAArchive, archive_bytes, recover_bytes
from .policy import ChannelProfile, recommend_codec_policy
from .simulator import SimulationConfig, simulate_channel
from .tiering import WorkloadProfile, recommend_storage_tier


def _archive(args: argparse.Namespace) -> None:
    source = Path(args.input)
    archive = archive_bytes(source.read_bytes(),
        ArchiveConfig(args.chunk_size, args.rs_nsym, args.parity_group_size, True))
    output = Path(args.output or f"{source}.oligoark.json")
    archive.save(output)
    print(json.dumps({"output": str(output), **archive.metadata}, indent=2))


def _recover(args: argparse.Namespace) -> None:
    archive = DNAArchive.load(args.archive)
    recovered = recover_bytes(archive)
    Path(args.output).write_bytes(recovered)
    print(json.dumps({"output": args.output, "sha256": archive.metadata["sha256"],
                      "verified": True}, indent=2))


def _simulate(args: argparse.Namespace) -> None:
    archive = DNAArchive.load(args.archive)
    cfg = SimulationConfig(args.substitution, args.insertion, args.deletion,
                           args.dropout, args.duplicate, args.seed)
    reads = simulate_channel(archive.strands, cfg)
    Path(args.output).write_text("\n".join(reads) + "\n", encoding="utf-8")
    print(json.dumps({"reads": len(reads), "output": args.output}, indent=2))


def _recover_reads(args: argparse.Namespace) -> None:
    archive = DNAArchive.load(args.archive)
    reads = [x.strip() for x in Path(args.reads).read_text(encoding="utf-8").splitlines()
             if x.strip()]
    Path(args.output).write_bytes(recover_bytes(archive, reads))
    print(json.dumps({"output": args.output, "verified": True}, indent=2))


def _recommend(args: argparse.Namespace) -> None:
    p = WorkloadProfile(args.retention_years, args.accesses_per_year, args.mutability,
                        args.retrieval_urgency, args.durability_priority, args.energy_priority)
    print(json.dumps(recommend_storage_tier(p).to_dict(), indent=2))


def _policy(args: argparse.Namespace) -> None:
    p = ChannelProfile(args.substitution, args.insertion, args.deletion, args.dropout)
    print(json.dumps(recommend_codec_policy(p).to_dict(), indent=2))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="oligoark",
        description="AI-native DNA archival storage research framework")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("archive", help="encode a file into an OligoArk DNA archive")
    p.add_argument("input"); p.add_argument("--output")
    p.add_argument("--chunk-size", type=int, default=96)
    p.add_argument("--rs-nsym", type=int, default=8)
    p.add_argument("--parity-group-size", type=int, default=8)
    p.set_defaults(func=_archive)

    p = sub.add_parser("recover", help="recover and verify a file from an archive")
    p.add_argument("archive"); p.add_argument("--output", required=True); p.set_defaults(func=_recover)

    p = sub.add_parser("simulate", help="simulate a noisy DNA channel")
    p.add_argument("archive"); p.add_argument("--output", default="reads.txt")
    p.add_argument("--substitution", type=float, default=0.0)
    p.add_argument("--insertion", type=float, default=0.0)
    p.add_argument("--deletion", type=float, default=0.0)
    p.add_argument("--dropout", type=float, default=0.0)
    p.add_argument("--duplicate", type=float, default=0.0)
    p.add_argument("--seed", type=int, default=7); p.set_defaults(func=_simulate)

    p = sub.add_parser("recover-reads", help="recover from simulated read sequences")
    p.add_argument("archive"); p.add_argument("reads")
    p.add_argument("--output", required=True); p.set_defaults(func=_recover_reads)

    p = sub.add_parser("recommend", help="recommend a storage tier")
    p.add_argument("--retention-years", type=float, required=True)
    p.add_argument("--accesses-per-year", type=float, default=0.0)
    p.add_argument("--mutability", type=float, default=0.0)
    p.add_argument("--retrieval-urgency", type=float, default=0.0)
    p.add_argument("--durability-priority", type=float, default=1.0)
    p.add_argument("--energy-priority", type=float, default=0.5); p.set_defaults(func=_recommend)

    p = sub.add_parser("policy", help="recommend an adaptive codec policy")
    p.add_argument("--substitution", type=float, default=0.0)
    p.add_argument("--insertion", type=float, default=0.0)
    p.add_argument("--deletion", type=float, default=0.0)
    p.add_argument("--dropout", type=float, default=0.0); p.set_defaults(func=_policy)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

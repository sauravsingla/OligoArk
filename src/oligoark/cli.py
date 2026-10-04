"""Command-line interface for OligoArk."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import cast

from .archive import (
    ArchiveConfig,
    DNAArchive,
    archive_bytes,
    archive_statistics,
    recover_bytes,
    recover_from_reads,
)
from .intelligence import plan_archive
from .logging_utils import configure_logging
from .policy import ChannelProfile, PolicyObjective, recommend_codec_policy
from .simulator import SimulationConfig, simulate_channel
from .tiering import EconomicAssumptions, WorkloadProfile, recommend_storage_tier


def _archive(args: argparse.Namespace) -> None:
    source = Path(args.input)
    config = ArchiveConfig(args.chunk_size, args.rs_nsym, args.parity_group_size, True)
    archive = archive_bytes(source.read_bytes(), config)
    output = Path(args.output or f"{source}.oligoark.json")
    archive.save(output)
    print(json.dumps({"output": str(output), **archive.metadata}, indent=2))


def _recover(args: argparse.Namespace) -> None:
    archive = DNAArchive.load(args.archive)
    recovered = recover_bytes(archive)
    output = Path(args.output)
    output.write_bytes(recovered)
    result = {
        "output": str(output),
        "sha256": archive.metadata["sha256"],
        "verified": True,
    }
    print(json.dumps(result, indent=2))


def _simulate(args: argparse.Namespace) -> None:
    archive = DNAArchive.load(args.archive)
    cfg = SimulationConfig(
        args.substitution,
        args.insertion,
        args.deletion,
        args.dropout,
        args.duplicate,
        args.seed,
    )
    reads = simulate_channel(archive.strands, cfg)
    Path(args.output).write_text("\n".join(reads) + "\n", encoding="utf-8")
    print(json.dumps({"reads": len(reads), "output": args.output}, indent=2))


def _recover_reads(args: argparse.Namespace) -> None:
    archive = DNAArchive.load(args.archive)
    text = Path(args.reads).read_text(encoding="utf-8")
    reads = [line.strip() for line in text.splitlines() if line.strip()]
    recovered, report = recover_from_reads(
        archive,
        reads,
        similarity_threshold=args.similarity_threshold,
    )
    Path(args.output).write_bytes(recovered)
    print(json.dumps({"output": args.output, **report.to_dict()}, indent=2))


def _inspect(args: argparse.Namespace) -> None:
    archive = DNAArchive.load(args.archive)
    print(json.dumps(archive_statistics(archive).to_dict(), indent=2))


def _load_economics(path: str | None) -> EconomicAssumptions | None:
    if path is None:
        return None
    parsed = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(parsed, dict):
        raise ValueError("Economic assumptions JSON must contain an object")
    return EconomicAssumptions.from_mapping(cast(dict[str, object], parsed))


def _workload(args: argparse.Namespace) -> WorkloadProfile:
    return WorkloadProfile(
        args.retention_years,
        args.accesses_per_year,
        args.mutability,
        args.retrieval_urgency,
        args.durability_priority,
        args.energy_priority,
        args.redundancy_priority,
        args.cost_priority,
    )


def _channel(args: argparse.Namespace) -> ChannelProfile:
    return ChannelProfile(
        args.substitution,
        args.insertion,
        args.deletion,
        args.dropout,
    )


def _recommend(args: argparse.Namespace) -> None:
    economics = _load_economics(args.economics_json)
    print(json.dumps(recommend_storage_tier(_workload(args), economics).to_dict(), indent=2))


def _policy(args: argparse.Namespace) -> None:
    objective = PolicyObjective(
        args.durability_priority,
        args.storage_overhead_priority,
        args.retrieval_speed_priority,
    )
    print(json.dumps(recommend_codec_policy(_channel(args), objective).to_dict(), indent=2))


def _plan(args: argparse.Namespace) -> None:
    economics = _load_economics(args.economics_json)
    result = plan_archive(_workload(args), _channel(args), economics=economics)
    print(json.dumps(result.to_dict(), indent=2))


def _add_workload_arguments(command: argparse.ArgumentParser) -> None:
    command.add_argument("--retention-years", type=float, required=True)
    command.add_argument("--accesses-per-year", type=float, default=0.0)
    command.add_argument("--mutability", type=float, default=0.0)
    command.add_argument("--retrieval-urgency", type=float, default=0.0)
    command.add_argument("--durability-priority", type=float, default=1.0)
    command.add_argument("--energy-priority", type=float, default=0.5)
    command.add_argument("--redundancy-priority", type=float, default=0.5)
    command.add_argument("--cost-priority", type=float, default=0.5)
    command.add_argument(
        "--economics-json",
        help="JSON file containing normalized per-tier storage/retrieval cost indices",
    )


def _add_channel_arguments(command: argparse.ArgumentParser) -> None:
    command.add_argument("--substitution", type=float, default=0.0)
    command.add_argument("--insertion", type=float, default=0.0)
    command.add_argument("--deletion", type=float, default=0.0)
    command.add_argument("--dropout", type=float, default=0.0)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="oligoark",
        description="AI-native DNA archival storage research framework",
    )
    parser.add_argument("--log-level", default="WARNING")
    sub = parser.add_subparsers(dest="command", required=True)

    command = sub.add_parser("archive", help="encode a file into an OligoArk DNA archive")
    command.add_argument("input")
    command.add_argument("--output")
    command.add_argument("--chunk-size", type=int, default=96)
    command.add_argument("--rs-nsym", type=int, default=8)
    command.add_argument("--parity-group-size", type=int, default=8)
    command.set_defaults(func=_archive)

    command = sub.add_parser("recover", help="recover and verify a file from an archive")
    command.add_argument("archive")
    command.add_argument("--output", required=True)
    command.set_defaults(func=_recover)

    command = sub.add_parser("simulate", help="simulate a noisy DNA channel")
    command.add_argument("archive")
    command.add_argument("--output", default="reads.txt")
    _add_channel_arguments(command)
    command.add_argument("--duplicate", type=float, default=0.0)
    command.add_argument("--seed", type=int, default=7)
    command.set_defaults(func=_simulate)

    command = sub.add_parser("recover-reads", help="reconstruct, recover, and verify reads")
    command.add_argument("archive")
    command.add_argument("reads")
    command.add_argument("--output", required=True)
    command.add_argument("--similarity-threshold", type=float, default=0.90)
    command.set_defaults(func=_recover_reads)

    command = sub.add_parser("inspect", help="show archive encoding and sequence statistics")
    command.add_argument("archive")
    command.set_defaults(func=_inspect)

    command = sub.add_parser("recommend", help="recommend a storage tier")
    _add_workload_arguments(command)
    command.set_defaults(func=_recommend)

    command = sub.add_parser("policy", help="recommend an adaptive codec policy")
    _add_channel_arguments(command)
    command.add_argument("--durability-priority", type=float, default=0.7)
    command.add_argument("--storage-overhead-priority", type=float, default=0.2)
    command.add_argument("--retrieval-speed-priority", type=float, default=0.1)
    command.set_defaults(func=_policy)

    command = sub.add_parser(
        "plan",
        help="create one explainable storage-tier and DNA-codec archival plan",
    )
    _add_workload_arguments(command)
    _add_channel_arguments(command)
    command.set_defaults(func=_plan)

    return parser


def main() -> None:
    args = build_parser().parse_args()
    configure_logging(args.log_level)
    args.func(args)


if __name__ == "__main__":
    main()

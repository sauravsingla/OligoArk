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
from .intelligence import (
    evaluate_optimized_archive_plan,
    optimize_archive_plan,
    plan_archive,
)
from .logging_utils import configure_logging
from .optimizer import CodecSearchSpace, OptimizationWeights
from .policy import ChannelProfile, PolicyObjective, recommend_codec_policy
from .simulator import SimulationConfig, simulate_channel
from .tiering import (
    EconomicAssumptions,
    LifecycleAssumptions,
    WorkloadProfile,
    recommend_storage_tier,
)
from .validation import compare_reconstruction_modes


def _archive(args: argparse.Namespace) -> None:
    source = Path(args.input)
    config = ArchiveConfig(
        chunk_size=args.chunk_size,
        rs_nsym=args.rs_nsym,
        parity_group_size=args.parity_group_size,
        adaptive_masks=not args.disable_adaptive_masks,
        redundancy_scheme=args.redundancy_scheme,
        fountain_redundancy=args.fountain_redundancy,
        min_gc_fraction=args.min_gc_fraction,
        max_gc_fraction=args.max_gc_fraction,
        max_homopolymer=args.max_homopolymer,
        mask_search_limit=args.mask_search_limit,
    )
    archive = archive_bytes(source.read_bytes(), config)
    output = Path(args.output or f"{source}.oligoark.json")
    archive.save(output)
    print(json.dumps({"output": str(output), **archive.metadata}, indent=2))


def _recover(args: argparse.Namespace) -> None:
    archive = DNAArchive.load(args.archive)
    recovered = recover_bytes(archive)
    output = Path(args.output)
    output.write_bytes(recovered)
    print(
        json.dumps(
            {"output": str(output), "sha256": archive.metadata["sha256"], "verified": True},
            indent=2,
        )
    )


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


def _read_reads(path: str) -> list[str]:
    text = Path(path).read_text(encoding="utf-8")
    return [line.strip() for line in text.splitlines() if line.strip()]


def _recover_reads(args: argparse.Namespace) -> None:
    archive = DNAArchive.load(args.archive)
    recovered, report = recover_from_reads(
        archive,
        _read_reads(args.reads),
        similarity_threshold=args.similarity_threshold,
    )
    Path(args.output).write_bytes(recovered)
    print(json.dumps({"output": args.output, **report.to_dict()}, indent=2))


def _diagnose_reconstruction(args: argparse.Namespace) -> None:
    result = compare_reconstruction_modes(
        DNAArchive.load(args.archive),
        _read_reads(args.reads),
        threshold=args.similarity_threshold,
    )
    print(json.dumps(result.to_dict(), indent=2))


def _inspect(args: argparse.Namespace) -> None:
    print(json.dumps(archive_statistics(DNAArchive.load(args.archive)).to_dict(), indent=2))


def _load_mapping(path: str, description: str) -> dict[str, object]:
    parsed = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(parsed, dict):
        raise ValueError(f"{description} JSON must contain an object")
    return cast(dict[str, object], parsed)


def _load_economics(path: str | None) -> EconomicAssumptions | None:
    return (
        None
        if path is None
        else EconomicAssumptions.from_mapping(_load_mapping(path, "Economic assumptions"))
    )


def _load_lifecycle(path: str | None) -> LifecycleAssumptions | None:
    if path is None:
        return None
    mapping = _load_mapping(path, "Lifecycle assumptions")
    tiers = mapping.get("tiers", mapping)
    if not isinstance(tiers, dict):
        raise ValueError("Lifecycle assumptions must contain a tier mapping")
    return LifecycleAssumptions.from_mapping(cast(dict[str, object], tiers))


def _load_weights(path: str | None) -> OptimizationWeights | None:
    if path is None:
        return None
    return OptimizationWeights.from_mapping(_load_mapping(path, "Optimization weights"))


def _parse_seeds(value: str) -> tuple[int, ...]:
    seeds = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    if not seeds:
        raise ValueError("seed list must not be empty")
    return seeds


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
        args.data_size_gb,
        args.expected_access_probability,
    )


def _channel(args: argparse.Namespace) -> ChannelProfile:
    return ChannelProfile(
        args.substitution,
        args.insertion,
        args.deletion,
        args.dropout,
    )


def _recommend(args: argparse.Namespace) -> None:
    result = recommend_storage_tier(
        _workload(args),
        _load_economics(args.economics_json),
        _load_lifecycle(args.lifecycle_json),
    )
    print(json.dumps(result.to_dict(), indent=2))


def _policy(args: argparse.Namespace) -> None:
    objective = PolicyObjective(
        args.durability_priority,
        args.storage_overhead_priority,
        args.retrieval_speed_priority,
    )
    print(json.dumps(recommend_codec_policy(_channel(args), objective).to_dict(), indent=2))


def _plan(args: argparse.Namespace) -> None:
    result = plan_archive(
        _workload(args),
        _channel(args),
        economics=_load_economics(args.economics_json),
        lifecycle=_load_lifecycle(args.lifecycle_json),
    )
    print(json.dumps(result.to_dict(), indent=2))


def _search_space(args: argparse.Namespace) -> CodecSearchSpace:
    return CodecSearchSpace(
        max_candidates=args.max_candidates,
        search_method=args.search_method,
        search_seed=args.search_seed,
    )


def _optimize_plan(args: argparse.Namespace) -> None:
    payload = Path(args.input).read_bytes()
    calibration_text = args.calibration_seeds or args.seeds
    calibration_seeds = _parse_seeds(calibration_text)
    economics = _load_economics(args.economics_json)
    lifecycle = _load_lifecycle(args.lifecycle_json)
    search_space = _search_space(args)
    weights = _load_weights(args.weights_json)
    if args.evaluation_seeds:
        held_out = evaluate_optimized_archive_plan(
            payload,
            _workload(args),
            _channel(args),
            calibration_seeds=calibration_seeds,
            evaluation_seeds=_parse_seeds(args.evaluation_seeds),
            duplicate_rate=args.duplicate,
            economics=economics,
            lifecycle=lifecycle,
            search_space=search_space,
            weights=weights,
        )
        print(json.dumps(held_out.to_dict(), indent=2))
        return

    optimized = optimize_archive_plan(
        payload,
        _workload(args),
        _channel(args),
        seeds=calibration_seeds,
        economics=economics,
        lifecycle=lifecycle,
        search_space=search_space,
        weights=weights,
        duplicate_rate=args.duplicate,
    )
    print(json.dumps(optimized.to_dict(), indent=2))


def _add_workload_arguments(command: argparse.ArgumentParser) -> None:
    command.add_argument("--retention-years", type=float, required=True)
    command.add_argument("--accesses-per-year", type=float, default=0.0)
    command.add_argument("--mutability", type=float, default=0.0)
    command.add_argument("--retrieval-urgency", type=float, default=0.0)
    command.add_argument("--durability-priority", type=float, default=1.0)
    command.add_argument("--energy-priority", type=float, default=0.5)
    command.add_argument("--redundancy-priority", type=float, default=0.5)
    command.add_argument("--cost-priority", type=float, default=0.5)
    command.add_argument("--data-size-gb", type=float, default=1.0)
    command.add_argument("--expected-access-probability", type=float)
    command.add_argument("--economics-json")
    command.add_argument("--lifecycle-json")


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
    command.add_argument(
        "--redundancy-scheme",
        choices=("none", "xor", "fountain", "hybrid"),
        default="xor",
    )
    command.add_argument("--fountain-redundancy", type=float, default=0.25)
    command.add_argument("--min-gc-fraction", type=float, default=0.35)
    command.add_argument("--max-gc-fraction", type=float, default=0.65)
    command.add_argument("--max-homopolymer", type=int, default=4)
    command.add_argument("--mask-search-limit", type=int, default=64)
    command.add_argument("--disable-adaptive-masks", action="store_true")
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

    command = sub.add_parser(
        "diagnose-reconstruction",
        help="compare direct, medoid-graph, and alignment-graph recovery",
    )
    command.add_argument("archive")
    command.add_argument("reads")
    command.add_argument("--similarity-threshold", type=float, default=0.90)
    command.set_defaults(func=_diagnose_reconstruction)

    command = sub.add_parser("inspect", help="show archive encoding and sequence statistics")
    command.add_argument("archive")
    command.set_defaults(func=_inspect)

    command = sub.add_parser("recommend", help="recommend a storage tier")
    _add_workload_arguments(command)
    command.set_defaults(func=_recommend)

    command = sub.add_parser("policy", help="recommend a heuristic codec policy")
    _add_channel_arguments(command)
    command.add_argument("--durability-priority", type=float, default=0.7)
    command.add_argument("--storage-overhead-priority", type=float, default=0.2)
    command.add_argument("--retrieval-speed-priority", type=float, default=0.1)
    command.set_defaults(func=_policy)

    command = sub.add_parser("plan", help="create an explainable heuristic tier + codec plan")
    _add_workload_arguments(command)
    _add_channel_arguments(command)
    command.set_defaults(func=_plan)

    command = sub.add_parser(
        "optimize-plan",
        help="calibrate a measured codec plan and optionally evaluate it on disjoint seeds",
    )
    command.add_argument("input")
    _add_workload_arguments(command)
    _add_channel_arguments(command)
    command.add_argument(
        "--seeds",
        default="2026,2027",
        help="backward-compatible alias for calibration seeds",
    )
    command.add_argument("--calibration-seeds")
    command.add_argument("--evaluation-seeds")
    command.add_argument("--duplicate", type=float, default=0.0)
    command.add_argument("--max-candidates", type=int, default=24)
    command.add_argument("--search-method", choices=("balanced", "full_grid"), default="balanced")
    command.add_argument("--search-seed", type=int, default=5050)
    command.add_argument("--weights-json")
    command.set_defaults(func=_optimize_plan)

    return parser


def main() -> None:
    args = build_parser().parse_args()
    configure_logging(args.log_level)
    args.func(args)


if __name__ == "__main__":
    main()

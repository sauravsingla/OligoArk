"""Reusable deterministic experiment and ablation framework for OligoArk."""

from __future__ import annotations

import hashlib
import math
import random
import statistics
import time
from dataclasses import asdict, dataclass

from .archive import ArchiveConfig, archive_bytes, recover_bytes, recover_from_reads
from .optimizer import CodecSearchSpace, optimize_codec
from .policy import ChannelProfile, recommend_codec_policy
from .simulator import SimulationConfig, simulate_channel
from .tiering import WorkloadProfile


@dataclass(frozen=True)
class ExperimentScenario:
    name: str
    channel: ChannelProfile
    duplicate_rate: float = 0.0


@dataclass(frozen=True)
class ExperimentProfile:
    seeds: tuple[int, ...]
    payload_sizes: tuple[int, ...]
    scenarios: tuple[ExperimentScenario, ...]
    strategies: tuple[str, ...] = (
        "fixed",
        "adaptive",
        "adaptive_fountain",
        "adaptive_graph",
        "combined",
    )
    optimizer_max_candidates: int = 12

    def validate(self) -> None:
        if not self.seeds or not self.payload_sizes or not self.scenarios:
            raise ValueError("experiment seeds, payload sizes, and scenarios must not be empty")
        if any(size <= 0 for size in self.payload_sizes):
            raise ValueError("payload sizes must be positive")
        valid = {
            "fixed",
            "adaptive",
            "adaptive_fountain",
            "adaptive_graph",
            "combined",
        }
        if any(strategy not in valid for strategy in self.strategies):
            raise ValueError(f"experiment strategies must be drawn from {sorted(valid)}")
        if self.optimizer_max_candidates < 1:
            raise ValueError("optimizer_max_candidates must be positive")


@dataclass(frozen=True)
class ExperimentRecord:
    strategy: str
    scenario: str
    seed: int
    payload_size: int
    payload_sha256: str
    recovered: bool
    encoded_nucleotides: int
    overhead_ratio: float
    strand_count: int
    read_count: int
    runtime_seconds: float
    graph_reconstruction_used: bool
    redundancy_scheme: str
    rs_nsym: int
    chunk_size: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ExperimentSummary:
    strategy: str
    scenario: str
    payload_size: int
    trials: int
    successes: int
    recovery_rate: float
    recovery_ci95_low: float
    recovery_ci95_high: float
    mean_overhead_ratio: float
    mean_runtime_seconds: float
    graph_use_rate: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def wilson_interval(successes: int, trials: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial recovery proportion."""
    if trials <= 0:
        raise ValueError("trials must be positive")
    if not 0 <= successes <= trials:
        raise ValueError("successes must be between 0 and trials")
    proportion = successes / trials
    denominator = 1.0 + z * z / trials
    center = (proportion + z * z / (2 * trials)) / denominator
    margin = (
        z
        * math.sqrt(
            proportion * (1.0 - proportion) / trials
            + z * z / (4 * trials * trials)
        )
        / denominator
    )
    return max(0.0, center - margin), min(1.0, center + margin)


def smoke_profile() -> ExperimentProfile:
    return ExperimentProfile(
        seeds=(2026,),
        payload_sizes=(256,),
        scenarios=(
            ExperimentScenario("clean", ChannelProfile()),
            ExperimentScenario(
                "substitution",
                ChannelProfile(substitution_rate=0.01),
                duplicate_rate=0.20,
            ),
            ExperimentScenario(
                "indel",
                ChannelProfile(insertion_rate=0.001, deletion_rate=0.001),
                duplicate_rate=0.35,
            ),
            ExperimentScenario(
                "dropout",
                ChannelProfile(dropout_rate=0.05),
                duplicate_rate=0.10,
            ),
        ),
        optimizer_max_candidates=8,
    )


def publication_profile() -> ExperimentProfile:
    return ExperimentProfile(
        seeds=(2026, 2027, 2028, 2029, 2030),
        payload_sizes=(512, 2048, 8192),
        scenarios=(
            ExperimentScenario("clean", ChannelProfile()),
            ExperimentScenario(
                "substitution-0.1pct",
                ChannelProfile(substitution_rate=0.001),
                duplicate_rate=0.10,
            ),
            ExperimentScenario(
                "substitution-1pct",
                ChannelProfile(substitution_rate=0.01),
                duplicate_rate=0.25,
            ),
            ExperimentScenario(
                "indel-low",
                ChannelProfile(insertion_rate=0.0005, deletion_rate=0.0005),
                duplicate_rate=0.35,
            ),
            ExperimentScenario(
                "indel-moderate",
                ChannelProfile(insertion_rate=0.0015, deletion_rate=0.0015),
                duplicate_rate=0.50,
            ),
            ExperimentScenario(
                "dropout-2pct",
                ChannelProfile(dropout_rate=0.02),
                duplicate_rate=0.10,
            ),
            ExperimentScenario(
                "dropout-10pct",
                ChannelProfile(dropout_rate=0.10),
                duplicate_rate=0.20,
            ),
            ExperimentScenario(
                "mixed",
                ChannelProfile(
                    substitution_rate=0.002,
                    insertion_rate=0.0005,
                    deletion_rate=0.0005,
                    dropout_rate=0.02,
                ),
                duplicate_rate=0.35,
            ),
        ),
        optimizer_max_candidates=24,
    )


def _adaptive_config(
    channel: ChannelProfile,
    *,
    redundancy_scheme: str = "xor",
) -> ArchiveConfig:
    policy = recommend_codec_policy(channel)
    return ArchiveConfig(
        chunk_size=policy.chunk_size,
        rs_nsym=policy.rs_nsym,
        parity_group_size=policy.parity_group_size,
        adaptive_masks=policy.adaptive_masks,
        redundancy_scheme=redundancy_scheme,
        fountain_redundancy=0.35 if redundancy_scheme in {"fountain", "hybrid"} else 0.0,
        mask_search_limit=96,
    )


def _strategy_config(
    strategy: str,
    payload: bytes,
    scenario: ExperimentScenario,
    workload: WorkloadProfile,
    *,
    optimizer_max_candidates: int,
    seed: int,
) -> tuple[ArchiveConfig, bool]:
    if strategy == "fixed":
        return ArchiveConfig(
            96,
            8,
            8,
            False,
            redundancy_scheme="xor",
            min_gc_fraction=0.0,
            max_gc_fraction=1.0,
            max_homopolymer=100,
            mask_search_limit=1,
        ), False
    if strategy == "adaptive":
        return _adaptive_config(scenario.channel, redundancy_scheme="xor"), False
    if strategy == "adaptive_fountain":
        return _adaptive_config(scenario.channel, redundancy_scheme="hybrid"), False
    if strategy == "adaptive_graph":
        return _adaptive_config(scenario.channel, redundancy_scheme="xor"), True
    if strategy == "combined":
        calibration = payload[: min(256, len(payload))]
        search = CodecSearchSpace(
            chunk_sizes=(48, 64, 96),
            rs_nsyms=(8, 16, 24),
            redundancy_schemes=("xor", "fountain", "hybrid"),
            parity_group_sizes=(3, 5),
            fountain_redundancies=(0.25, 0.40),
            reconstruction_modes=("direct", "graph"),
            max_candidates=optimizer_max_candidates,
        )
        result = optimize_codec(
            calibration,
            scenario.channel,
            workload,
            search_space=search,
            seeds=(seed,),
        )
        return result.best_config, result.reconstruction_mode == "graph"
    raise ValueError(f"unknown experiment strategy: {strategy}")


def _run_one(
    payload: bytes,
    scenario: ExperimentScenario,
    strategy: str,
    *,
    seed: int,
    optimizer_max_candidates: int,
) -> ExperimentRecord:
    workload = WorkloadProfile(
        retention_years=100,
        accesses_per_year=0.1,
        mutability=0.0,
        retrieval_urgency=0.2,
        durability_priority=1.0,
        energy_priority=0.5,
        redundancy_priority=0.8,
        cost_priority=0.5,
        data_size_gb=max(1e-9, len(payload) / 1_000_000_000),
    )
    started = time.perf_counter()
    config, use_graph = _strategy_config(
        strategy,
        payload,
        scenario,
        workload,
        optimizer_max_candidates=optimizer_max_candidates,
        seed=seed,
    )
    archive = archive_bytes(payload, config)
    reads = simulate_channel(
        archive.strands,
        SimulationConfig(
            substitution_rate=scenario.channel.substitution_rate,
            insertion_rate=scenario.channel.insertion_rate,
            deletion_rate=scenario.channel.deletion_rate,
            dropout_rate=scenario.channel.dropout_rate,
            duplicate_rate=scenario.duplicate_rate,
            seed=seed,
        ),
    )
    recovered = False
    graph_used = False
    try:
        if use_graph:
            decoded, report = recover_from_reads(archive, reads)
            graph_used = report.graph_reconstruction_used
        else:
            decoded = recover_bytes(archive, reads)
        recovered = decoded == payload
    except ValueError:
        recovered = False
    measured = archive.metadata.get("measured")
    if not isinstance(measured, dict):
        raise ValueError("archive measured metadata must be an object")
    encoded = int(measured["encoded_nucleotides"])
    return ExperimentRecord(
        strategy=strategy,
        scenario=scenario.name,
        seed=seed,
        payload_size=len(payload),
        payload_sha256=hashlib.sha256(payload).hexdigest(),
        recovered=recovered,
        encoded_nucleotides=encoded,
        overhead_ratio=round(encoded / max(1, len(payload) * 4), 6),
        strand_count=len(archive.strands),
        read_count=len(reads),
        runtime_seconds=round(time.perf_counter() - started, 6),
        graph_reconstruction_used=graph_used,
        redundancy_scheme=config.redundancy_scheme,
        rs_nsym=config.rs_nsym,
        chunk_size=config.chunk_size,
    )


def run_experiments(profile: ExperimentProfile) -> list[ExperimentRecord]:
    profile.validate()
    records: list[ExperimentRecord] = []
    for payload_size in profile.payload_sizes:
        for seed in profile.seeds:
            rng = random.Random(seed * 1000003 + payload_size)
            payload = bytes(rng.randrange(256) for _ in range(payload_size))
            for scenario in profile.scenarios:
                for strategy in profile.strategies:
                    records.append(
                        _run_one(
                            payload,
                            scenario,
                            strategy,
                            seed=seed,
                            optimizer_max_candidates=profile.optimizer_max_candidates,
                        )
                    )
    return records


def aggregate_experiments(records: list[ExperimentRecord]) -> list[ExperimentSummary]:
    grouped: dict[tuple[str, str, int], list[ExperimentRecord]] = {}
    for record in records:
        grouped.setdefault(
            (record.strategy, record.scenario, record.payload_size),
            [],
        ).append(record)

    summaries: list[ExperimentSummary] = []
    for (strategy, scenario, payload_size), group in sorted(grouped.items()):
        successes = sum(record.recovered for record in group)
        low, high = wilson_interval(successes, len(group))
        summaries.append(
            ExperimentSummary(
                strategy=strategy,
                scenario=scenario,
                payload_size=payload_size,
                trials=len(group),
                successes=successes,
                recovery_rate=round(successes / len(group), 6),
                recovery_ci95_low=round(low, 6),
                recovery_ci95_high=round(high, 6),
                mean_overhead_ratio=round(
                    statistics.fmean(record.overhead_ratio for record in group),
                    6,
                ),
                mean_runtime_seconds=round(
                    statistics.fmean(record.runtime_seconds for record in group),
                    6,
                ),
                graph_use_rate=round(
                    statistics.fmean(
                        1.0 if record.graph_reconstruction_used else 0.0
                        for record in group
                    ),
                    6,
                ),
            )
        )
    return summaries

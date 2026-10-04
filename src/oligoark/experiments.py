"""Held-out deterministic experiment and ablation framework for OligoArk."""

from __future__ import annotations

import hashlib
import math
import random
import statistics
import time
from dataclasses import asdict, dataclass

from .archive import ArchiveConfig, archive_bytes, recover_bytes, recover_from_reads
from .optimizer import CodecSearchSpace, OptimizationResult, optimize_codec
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
    """Experiment profile with disjoint calibration and held-out evaluation seeds.

    The seeds field is retained for backward compatibility and means evaluation/test
    seeds in v0.5.
    """

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
    calibration_seeds: tuple[int, ...] = (9101, 9102, 9103)
    optimizer_search_method: str = "balanced"
    optimizer_search_seed: int = 5050

    @property
    def evaluation_seeds(self) -> tuple[int, ...]:
        return self.seeds

    def validate(self) -> None:
        if not self.seeds or not self.calibration_seeds:
            raise ValueError("calibration and evaluation seeds must not be empty")
        if set(self.seeds) & set(self.calibration_seeds):
            raise ValueError("calibration and evaluation seeds must be disjoint")
        if not self.payload_sizes or not self.scenarios:
            raise ValueError("payload sizes and scenarios must not be empty")
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
        if self.optimizer_search_method not in {"balanced", "full_grid"}:
            raise ValueError("optimizer_search_method must be balanced or full_grid")


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
    parity_group_size: int
    adaptive_masks: bool
    substitution_rate: float
    insertion_rate: float
    deletion_rate: float
    dropout_rate: float
    duplicate_rate: float
    calibration_seeds: tuple[int, ...] = ()
    selection_score: float | None = None
    selection_search_method: str | None = None
    selection_calibration_recovery_rate: float | None = None

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
    graph_rescue_rate: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class StrategyEffectSummary:
    strategy: str
    baseline_strategy: str
    scenario: str
    payload_size: int
    paired_trials: int
    recovery_rate_difference: float
    overhead_ratio_difference: float
    runtime_seconds_difference: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class CalibrationRecord:
    scenario: str
    payload_size: int
    calibration_payload_size: int
    optimization: OptimizationResult

    def to_dict(self) -> dict[str, object]:
        return {
            "scenario": self.scenario,
            "payload_size": self.payload_size,
            "calibration_payload_size": self.calibration_payload_size,
            "optimization": self.optimization.to_dict(),
        }


@dataclass(frozen=True)
class ExperimentBundle:
    records: tuple[ExperimentRecord, ...]
    calibrations: tuple[CalibrationRecord, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "records": [record.to_dict() for record in self.records],
            "calibrations": [record.to_dict() for record in self.calibrations],
        }


def wilson_interval(successes: int, trials: int, z: float = 1.96) -> tuple[float, float]:
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
        seeds=(2026, 2027),
        calibration_seeds=(9001, 9002),
        payload_sizes=(256,),
        scenarios=(
            ExperimentScenario("clean", ChannelProfile()),
            ExperimentScenario(
                "substitution", ChannelProfile(substitution_rate=0.01), duplicate_rate=0.20
            ),
            ExperimentScenario(
                "indel",
                ChannelProfile(insertion_rate=0.001, deletion_rate=0.001),
                duplicate_rate=0.35,
            ),
            ExperimentScenario("dropout", ChannelProfile(dropout_rate=0.05), duplicate_rate=0.10),
        ),
        optimizer_max_candidates=8,
    )


def publication_profile() -> ExperimentProfile:
    return ExperimentProfile(
        seeds=(2026, 2027, 2028, 2029, 2030, 2031, 2032, 2033),
        calibration_seeds=(9201, 9202, 9203, 9204),
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
                "dropout-2pct", ChannelProfile(dropout_rate=0.02), duplicate_rate=0.10
            ),
            ExperimentScenario(
                "dropout-10pct", ChannelProfile(dropout_rate=0.10), duplicate_rate=0.20
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
        optimizer_search_method="balanced",
        optimizer_search_seed=5050,
    )


def _workload(payload_size: int) -> WorkloadProfile:
    return WorkloadProfile(
        retention_years=100,
        accesses_per_year=0.1,
        mutability=0.0,
        retrieval_urgency=0.2,
        durability_priority=1.0,
        energy_priority=0.5,
        redundancy_priority=0.8,
        cost_priority=0.5,
        data_size_gb=max(1e-9, payload_size / 1_000_000_000),
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


def _optimizer_search(profile: ExperimentProfile) -> CodecSearchSpace:
    return CodecSearchSpace(
        chunk_sizes=(48, 64, 96),
        rs_nsyms=(8, 16, 24),
        redundancy_schemes=("xor", "fountain", "hybrid"),
        parity_group_sizes=(3, 5),
        fountain_redundancies=(0.25, 0.40),
        reconstruction_modes=("direct", "graph"),
        max_candidates=profile.optimizer_max_candidates,
        search_method=profile.optimizer_search_method,
        search_seed=profile.optimizer_search_seed,
    )


def _calibrate_combined(
    payload: bytes,
    scenario: ExperimentScenario,
    profile: ExperimentProfile,
) -> OptimizationResult:
    calibration = payload[: min(256, len(payload))]
    return optimize_codec(
        calibration,
        scenario.channel,
        _workload(len(payload)),
        search_space=_optimizer_search(profile),
        seeds=profile.calibration_seeds,
        duplicate_rate=scenario.duplicate_rate,
    )


def _strategy_config(
    strategy: str,
    scenario: ExperimentScenario,
    combined: OptimizationResult | None,
) -> tuple[ArchiveConfig, bool, float | None, str | None, float | None]:
    if strategy == "fixed":
        return (
            ArchiveConfig(
                96,
                8,
                8,
                False,
                redundancy_scheme="xor",
                min_gc_fraction=0.0,
                max_gc_fraction=1.0,
                max_homopolymer=100,
                mask_search_limit=1,
            ),
            False,
            None,
            None,
            None,
        )
    if strategy == "adaptive":
        return (
            _adaptive_config(scenario.channel, redundancy_scheme="xor"),
            False,
            None,
            None,
            None,
        )
    if strategy == "adaptive_fountain":
        return (
            _adaptive_config(scenario.channel, redundancy_scheme="hybrid"),
            False,
            None,
            None,
            None,
        )
    if strategy == "adaptive_graph":
        return (
            _adaptive_config(scenario.channel, redundancy_scheme="xor"),
            True,
            None,
            None,
            None,
        )
    if strategy == "combined":
        if combined is None:
            raise ValueError("combined strategy requires calibration result")
        winner = next(
            evaluation
            for evaluation in combined.evaluations
            if evaluation.rejected_reason is None
            and evaluation.config == combined.best_config
            and evaluation.reconstruction_mode == combined.reconstruction_mode
            and evaluation.score == combined.best_score
        )
        return (
            combined.best_config,
            combined.reconstruction_mode == "graph",
            combined.best_score,
            combined.search_method,
            winner.recovery_rate,
        )
    raise ValueError(f"unknown experiment strategy: {strategy}")


def _payload(payload_size: int, seed: int) -> bytes:
    rng = random.Random(seed * 1_000_003 + payload_size)
    return bytes(rng.randrange(256) for _ in range(payload_size))


def _run_one(
    payload: bytes,
    scenario: ExperimentScenario,
    strategy: str,
    config: ArchiveConfig,
    use_graph: bool,
    *,
    seed: int,
    calibration_seeds: tuple[int, ...],
    selection_score: float | None,
    selection_search_method: str | None,
    selection_calibration_recovery_rate: float | None,
) -> ExperimentRecord:
    started = time.perf_counter()
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
    graph_rescue = False
    try:
        if use_graph:
            decoded, report = recover_from_reads(archive, reads)
            graph_rescue = report.rescue_changed_result
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
        graph_reconstruction_used=graph_rescue,
        redundancy_scheme=config.redundancy_scheme,
        rs_nsym=config.rs_nsym,
        chunk_size=config.chunk_size,
        parity_group_size=config.parity_group_size,
        adaptive_masks=config.adaptive_masks,
        substitution_rate=scenario.channel.substitution_rate,
        insertion_rate=scenario.channel.insertion_rate,
        deletion_rate=scenario.channel.deletion_rate,
        dropout_rate=scenario.channel.dropout_rate,
        duplicate_rate=scenario.duplicate_rate,
        calibration_seeds=calibration_seeds if strategy == "combined" else (),
        selection_score=selection_score,
        selection_search_method=selection_search_method,
        selection_calibration_recovery_rate=selection_calibration_recovery_rate,
    )


def run_experiment_bundle(profile: ExperimentProfile) -> ExperimentBundle:
    """Run held-out trials and retain full optimizer calibration evidence."""
    profile.validate()
    records: list[ExperimentRecord] = []
    calibrations: list[CalibrationRecord] = []
    for payload_size in profile.payload_sizes:
        calibration_payload = _payload(payload_size, profile.calibration_seeds[0])
        combined_by_scenario: dict[str, OptimizationResult] = {}
        if "combined" in profile.strategies:
            for scenario in profile.scenarios:
                result = _calibrate_combined(calibration_payload, scenario, profile)
                combined_by_scenario[scenario.name] = result
                calibrations.append(
                    CalibrationRecord(
                        scenario=scenario.name,
                        payload_size=payload_size,
                        calibration_payload_size=min(256, len(calibration_payload)),
                        optimization=result,
                    )
                )
        for seed in profile.evaluation_seeds:
            payload = _payload(payload_size, seed)
            for scenario in profile.scenarios:
                combined = combined_by_scenario.get(scenario.name)
                for strategy in profile.strategies:
                    config, use_graph, score, method, calibration_recovery = _strategy_config(
                        strategy, scenario, combined
                    )
                    records.append(
                        _run_one(
                            payload,
                            scenario,
                            strategy,
                            config,
                            use_graph,
                            seed=seed,
                            calibration_seeds=profile.calibration_seeds,
                            selection_score=score,
                            selection_search_method=method,
                            selection_calibration_recovery_rate=calibration_recovery,
                        )
                    )
    return ExperimentBundle(tuple(records), tuple(calibrations))


def run_experiments(profile: ExperimentProfile) -> list[ExperimentRecord]:
    """Backward-compatible convenience wrapper returning only held-out trial records."""
    return list(run_experiment_bundle(profile).records)


def aggregate_experiments(records: list[ExperimentRecord]) -> list[ExperimentSummary]:
    grouped: dict[tuple[str, str, int], list[ExperimentRecord]] = {}
    for record in records:
        grouped.setdefault((record.strategy, record.scenario, record.payload_size), []).append(
            record
        )

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
                    statistics.fmean(record.overhead_ratio for record in group), 6
                ),
                mean_runtime_seconds=round(
                    statistics.fmean(record.runtime_seconds for record in group), 6
                ),
                graph_rescue_rate=round(
                    statistics.fmean(
                        1.0 if record.graph_reconstruction_used else 0.0 for record in group
                    ),
                    6,
                ),
            )
        )
    return summaries


def paired_strategy_effects(
    records: list[ExperimentRecord],
    *,
    baseline_strategy: str = "fixed",
) -> list[StrategyEffectSummary]:
    """Paired mean differences on identical scenario/payload/seed realizations."""
    index = {
        (record.strategy, record.scenario, record.payload_size, record.seed): record
        for record in records
    }
    groups = sorted(
        {(record.strategy, record.scenario, record.payload_size) for record in records}
    )
    effects: list[StrategyEffectSummary] = []
    for strategy, scenario, payload_size in groups:
        if strategy == baseline_strategy:
            continue
        pairs: list[tuple[ExperimentRecord, ExperimentRecord]] = []
        for record in records:
            if (
                record.strategy == strategy
                and record.scenario == scenario
                and record.payload_size == payload_size
            ):
                baseline = index.get((baseline_strategy, scenario, payload_size, record.seed))
                if baseline is not None:
                    pairs.append((record, baseline))
        if not pairs:
            continue
        effects.append(
            StrategyEffectSummary(
                strategy=strategy,
                baseline_strategy=baseline_strategy,
                scenario=scenario,
                payload_size=payload_size,
                paired_trials=len(pairs),
                recovery_rate_difference=round(
                    statistics.fmean(
                        float(left.recovered) - float(right.recovered) for left, right in pairs
                    ),
                    6,
                ),
                overhead_ratio_difference=round(
                    statistics.fmean(
                        left.overhead_ratio - right.overhead_ratio for left, right in pairs
                    ),
                    6,
                ),
                runtime_seconds_difference=round(
                    statistics.fmean(
                        left.runtime_seconds - right.runtime_seconds for left, right in pairs
                    ),
                    6,
                ),
            )
        )
    return effects

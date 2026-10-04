"""AI-Native Archival Intelligence Layer.

This module composes OligoArk's tiering, deterministic policy, empirical learning, and real
search-based codec optimisation into explainable public planning APIs.
"""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass

from .archive import archive_bytes, recover_bytes, recover_from_reads
from .learning import EmpiricalPolicyModel
from .optimizer import (
    CodecSearchSpace,
    LifecycleObjectiveInputs,
    OptimizationResult,
    OptimizationWeights,
    optimize_codec,
)
from .policy import ChannelProfile, CodecPolicy, PolicyObjective, recommend_codec_policy
from .simulator import SimulationConfig, simulate_channel
from .tiering import (
    EconomicAssumptions,
    LifecycleAssumptions,
    TierRecommendation,
    WorkloadProfile,
    recommend_storage_tier,
)


@dataclass(frozen=True)
class ArchivalIntelligencePlan:
    """Combined explainable storage-tier and heuristic codec-policy recommendation."""

    tier: TierRecommendation
    codec_policy: CodecPolicy
    policy_source: str
    learned_confidence: float | None
    rationale: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class OptimizedArchivalPlan:
    """Storage-tier recommendation plus measured search-based codec optimisation."""

    tier: TierRecommendation
    optimization: OptimizationResult
    selected_redundancy_scheme: str
    selected_reconstruction_strategy: str
    selected_sequence_constraints: dict[str, object]
    rationale: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class HeldOutTrial:
    seed: int
    recovered: bool
    runtime_seconds: float
    graph_rescue_used: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class HeldOutOptimizedPlan:
    plan: OptimizedArchivalPlan
    calibration_seeds: tuple[int, ...]
    evaluation_seeds: tuple[int, ...]
    trials: tuple[HeldOutTrial, ...]
    held_out_recovery_rate: float
    mean_runtime_seconds: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def objective_from_workload(workload: WorkloadProfile) -> PolicyObjective:
    """Derive deterministic codec objectives from workload priorities."""
    workload.validate()
    overhead_priority = min(1.0, (workload.cost_priority + workload.energy_priority) / 2.0)
    return PolicyObjective(
        durability_priority=workload.durability_priority,
        storage_overhead_priority=overhead_priority,
        retrieval_speed_priority=workload.retrieval_urgency,
    )


def plan_archive(
    workload: WorkloadProfile,
    channel: ChannelProfile,
    *,
    economics: EconomicAssumptions | None = None,
    lifecycle: LifecycleAssumptions | None = None,
    objective: PolicyObjective | None = None,
    empirical_model: EmpiricalPolicyModel | None = None,
) -> ArchivalIntelligencePlan:
    """Create one explainable tier + heuristic/empirical codec plan.

    This function is intentionally lightweight and does not call itself an optimizer.
    Use optimize_archive_plan when candidates should actually be encoded, simulated,
    recovered, verified, measured, and ranked.
    """
    workload.validate()
    channel.validate()
    tier = recommend_storage_tier(workload, economics, lifecycle)
    resolved_objective = objective or objective_from_workload(workload)

    if empirical_model is None:
        codec_policy = recommend_codec_policy(channel, resolved_objective)
        policy_source = "deterministic-heuristic"
        learned_confidence = None
        source_rationale = "codec policy selected by deterministic channel/workload rules"
    else:
        learned = empirical_model.recommend(channel)
        codec_policy = learned.policy
        policy_source = "empirical-instance"
        learned_confidence = learned.confidence
        source_rationale = (
            "codec policy selected from caller-supplied observations; confidence is relative, "
            "not a calibrated probability"
        )

    rationale = (
        f"storage tier={tier.recommended_tier}",
        f"codec policy source={policy_source}",
        (
            "workload-derived objective: "
            f"durability={resolved_objective.durability_priority:.3f}, "
            f"overhead={resolved_objective.storage_overhead_priority:.3f}, "
            f"retrieval={resolved_objective.retrieval_speed_priority:.3f}"
        ),
        source_rationale,
        "this lightweight plan is heuristic; optimize_archive_plan performs measured search",
    )
    return ArchivalIntelligencePlan(
        tier=tier,
        codec_policy=codec_policy,
        policy_source=policy_source,
        learned_confidence=learned_confidence,
        rationale=rationale,
    )


def optimize_archive_plan(
    payload: bytes,
    workload: WorkloadProfile,
    channel: ChannelProfile,
    *,
    economics: EconomicAssumptions | None = None,
    lifecycle: LifecycleAssumptions | None = None,
    search_space: CodecSearchSpace | None = None,
    weights: OptimizationWeights | None = None,
    seeds: tuple[int, ...] = (2026, 2027),
    duplicate_rate: float | None = None,
) -> OptimizedArchivalPlan:
    """Search real software candidates and combine the winner with storage-tier planning."""
    workload.validate()
    channel.validate()
    tier = recommend_storage_tier(workload, economics, lifecycle)
    lifecycle_objective: LifecycleObjectiveInputs | None = None
    if tier.lifecycle_estimates is not None:
        estimate = tier.lifecycle_estimates[tier.recommended_tier]
        lifecycle_objective = LifecycleObjectiveInputs(
            storage_cost=estimate.total_storage_cost,
            retrieval_cost=estimate.total_retrieval_cost,
            energy_kwh=estimate.total_energy_kwh,
            retrieval_latency_hours=estimate.expected_retrieval_latency_hours,
        )
    optimization = optimize_codec(
        payload,
        channel,
        workload,
        search_space=search_space,
        weights=weights,
        seeds=seeds,
        lifecycle=lifecycle_objective,
        duplicate_rate=duplicate_rate,
    )
    config = optimization.best_config
    constraints: dict[str, object] = {
        "min_gc_fraction": config.min_gc_fraction,
        "max_gc_fraction": config.max_gc_fraction,
        "max_homopolymer": config.max_homopolymer,
        "mask_search_limit": config.mask_search_limit,
    }
    rationale = (
        f"tier={tier.recommended_tier}",
        f"searched codec winner score={optimization.best_score:.8f}",
        f"redundancy={config.redundancy_scheme}",
        f"reconstruction={optimization.reconstruction_mode}",
        (
            "constraints="
            f"GC[{config.min_gc_fraction:.2f},{config.max_gc_fraction:.2f}], "
            f"homopolymer<={config.max_homopolymer}"
        ),
        (
            "codec objective includes caller-supplied lifecycle terms for selected tier"
            if lifecycle_objective is not None
            else "codec objective omits lifecycle physical terms because none were supplied"
        ),
        "all recovery successes are SHA-256-verified software results",
    )
    return OptimizedArchivalPlan(
        tier=tier,
        optimization=optimization,
        selected_redundancy_scheme=config.redundancy_scheme,
        selected_reconstruction_strategy=optimization.reconstruction_mode,
        selected_sequence_constraints=constraints,
        rationale=rationale,
    )


def evaluate_optimized_archive_plan(
    payload: bytes,
    workload: WorkloadProfile,
    channel: ChannelProfile,
    *,
    calibration_seeds: tuple[int, ...],
    evaluation_seeds: tuple[int, ...],
    economics: EconomicAssumptions | None = None,
    lifecycle: LifecycleAssumptions | None = None,
    search_space: CodecSearchSpace | None = None,
    weights: OptimizationWeights | None = None,
    duplicate_rate: float = 0.0,
) -> HeldOutOptimizedPlan:
    """Freeze a policy on calibration seeds, then evaluate it on unseen seeds."""
    if not calibration_seeds or not evaluation_seeds:
        raise ValueError("calibration and evaluation seeds must not be empty")
    if set(calibration_seeds) & set(evaluation_seeds):
        raise ValueError("calibration and evaluation seeds must be disjoint")
    if not 0 <= duplicate_rate <= 1:
        raise ValueError("duplicate_rate must be between 0 and 1")

    plan = optimize_archive_plan(
        payload,
        workload,
        channel,
        economics=economics,
        lifecycle=lifecycle,
        search_space=search_space,
        weights=weights,
        seeds=calibration_seeds,
        duplicate_rate=duplicate_rate,
    )
    config = plan.optimization.best_config
    use_graph = plan.optimization.reconstruction_mode == "graph"
    trials: list[HeldOutTrial] = []
    for seed in evaluation_seeds:
        started = time.perf_counter()
        archive = archive_bytes(payload, config)
        reads = simulate_channel(
            archive.strands,
            SimulationConfig(
                substitution_rate=channel.substitution_rate,
                insertion_rate=channel.insertion_rate,
                deletion_rate=channel.deletion_rate,
                dropout_rate=channel.dropout_rate,
                duplicate_rate=duplicate_rate,
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
        trials.append(
            HeldOutTrial(
                seed=seed,
                recovered=recovered,
                runtime_seconds=round(time.perf_counter() - started, 8),
                graph_rescue_used=graph_rescue,
            )
        )

    return HeldOutOptimizedPlan(
        plan=plan,
        calibration_seeds=calibration_seeds,
        evaluation_seeds=evaluation_seeds,
        trials=tuple(trials),
        held_out_recovery_rate=round(
            sum(trial.recovered for trial in trials) / len(trials),
            6,
        ),
        mean_runtime_seconds=round(
            sum(trial.runtime_seconds for trial in trials) / len(trials),
            8,
        ),
    )

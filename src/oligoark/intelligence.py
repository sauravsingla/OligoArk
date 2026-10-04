"""AI-Native Archival Intelligence Layer.

This module composes OligoArk's tiering, deterministic policy, empirical learning, and real
search-based codec optimisation into explainable public planning APIs.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from .learning import EmpiricalPolicyModel
from .optimizer import (
    CodecSearchSpace,
    OptimizationResult,
    OptimizationWeights,
    optimize_codec,
)
from .policy import ChannelProfile, CodecPolicy, PolicyObjective, recommend_codec_policy
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
) -> OptimizedArchivalPlan:
    """Search real software candidates and combine the winner with storage-tier planning."""
    workload.validate()
    channel.validate()
    tier = recommend_storage_tier(workload, economics, lifecycle)
    optimization = optimize_codec(
        payload,
        channel,
        workload,
        search_space=search_space,
        weights=weights,
        seeds=seeds,
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

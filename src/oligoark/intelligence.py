"""AI-Native Archival Intelligence Layer.

This module composes OligoArk's explainable tiering and codec-policy systems into one public
planning API. It remains deterministic by default and can optionally use an empirical policy
model trained from explicitly supplied observations.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from .learning import EmpiricalPolicyModel
from .policy import ChannelProfile, CodecPolicy, PolicyObjective, recommend_codec_policy
from .tiering import (
    EconomicAssumptions,
    TierRecommendation,
    WorkloadProfile,
    recommend_storage_tier,
)


@dataclass(frozen=True)
class ArchivalIntelligencePlan:
    """Combined explainable storage-tier and codec-policy recommendation."""

    tier: TierRecommendation
    codec_policy: CodecPolicy
    policy_source: str
    learned_confidence: float | None
    rationale: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def objective_from_workload(workload: WorkloadProfile) -> PolicyObjective:
    """Derive codec objectives from workload priorities using a transparent mapping."""
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
    objective: PolicyObjective | None = None,
    empirical_model: EmpiricalPolicyModel | None = None,
) -> ArchivalIntelligencePlan:
    """Create one explainable archival plan from workload and channel requirements."""
    workload.validate()
    channel.validate()
    tier = recommend_storage_tier(workload, economics)
    resolved_objective = objective or objective_from_workload(workload)

    if empirical_model is None:
        codec_policy = recommend_codec_policy(channel, resolved_objective)
        policy_source = "deterministic"
        learned_confidence = None
        source_rationale = "codec policy selected by deterministic channel/workload heuristics"
    else:
        learned = empirical_model.recommend(channel)
        codec_policy = learned.policy
        policy_source = "empirical"
        learned_confidence = learned.confidence
        source_rationale = (
            "codec policy selected by caller-supplied empirical observations; confidence is "
            "relative, not a calibrated probability"
        )

    rationale = (
        f"storage tier={tier.recommended_tier}",
        f"codec policy source={policy_source}",
        (
            "workload-derived codec objective: "
            f"durability={resolved_objective.durability_priority:.3f}, "
            f"overhead={resolved_objective.storage_overhead_priority:.3f}, "
            f"retrieval={resolved_objective.retrieval_speed_priority:.3f}"
        ),
        source_rationale,
        "dna_future remains an experimental scenario tier, not a present-day superiority claim",
    )
    return ArchivalIntelligencePlan(
        tier=tier,
        codec_policy=codec_policy,
        policy_source=policy_source,
        learned_confidence=learned_confidence,
        rationale=rationale,
    )

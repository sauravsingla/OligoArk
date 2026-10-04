from oligoark.intelligence import objective_from_workload, plan_archive
from oligoark.learning import EmpiricalPolicyModel, PolicyObservation
from oligoark.policy import ChannelProfile, CodecPolicy
from oligoark.tiering import WorkloadProfile


def _workload() -> WorkloadProfile:
    return WorkloadProfile(
        retention_years=100,
        accesses_per_year=0.1,
        mutability=0.0,
        retrieval_urgency=0.2,
        durability_priority=0.95,
        energy_priority=0.8,
        redundancy_priority=0.8,
        cost_priority=0.6,
    )


def test_workload_is_mapped_into_codec_objective() -> None:
    objective = objective_from_workload(_workload())
    assert objective.durability_priority == 0.95
    assert objective.storage_overhead_priority == 0.7
    assert objective.retrieval_speed_priority == 0.2


def test_integrated_plan_is_explainable_and_deterministic() -> None:
    channel = ChannelProfile(substitution_rate=0.01, dropout_rate=0.05)
    first = plan_archive(_workload(), channel)
    second = plan_archive(_workload(), channel)
    assert first == second
    assert first.policy_source == "deterministic-heuristic"
    assert first.tier.recommended_tier in first.tier.scores
    assert first.codec_policy.rs_nsym >= 16
    assert "heuristic" in first.rationale[-1]


def test_integrated_plan_can_use_empirical_policy_model() -> None:
    lean = CodecPolicy(96, 8, 8, True, ("lean",))
    robust = CodecPolicy(48, 24, 3, True, ("robust",))
    observations = [
        PolicyObservation(ChannelProfile(), lean, True, 1000, 0.1),
        PolicyObservation(ChannelProfile(0.03, 0, 0, 0.1), lean, False, 1000, 0.1),
        PolicyObservation(ChannelProfile(0.03, 0, 0, 0.1), robust, True, 1800, 0.2),
    ]
    model = EmpiricalPolicyModel().fit(observations)
    plan = plan_archive(
        _workload(),
        ChannelProfile(0.028, 0, 0, 0.09),
        empirical_model=model,
    )
    assert plan.policy_source == "empirical-instance"
    assert plan.codec_policy == robust
    assert plan.learned_confidence is not None

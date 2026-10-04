import pytest

from oligoark.tiering import (
    LifecycleAssumptions,
    WorkloadProfile,
    recommend_storage_tier,
)


def _lifecycle() -> LifecycleAssumptions:
    base = {
        "storage_cost_per_gb_year": 1.0,
        "retrieval_cost_per_gb": 1.0,
        "idle_energy_kwh_per_tb_year": 1.0,
        "retrieval_energy_kwh_per_gb": 1.0,
        "retrieval_latency_hours": 1.0,
    }
    return LifecycleAssumptions.from_mapping(
        {
            "ssd": {**base, "retrieval_latency_hours": 0.1},
            "object_archive": {**base, "storage_cost_per_gb_year": 0.6},
            "tape": {**base, "storage_cost_per_gb_year": 0.4},
            "dna_future": {**base, "storage_cost_per_gb_year": 2.0},
        }
    )


def test_tier_score_breakdown_reproduces_each_score() -> None:
    result = recommend_storage_tier(
        WorkloadProfile(
            retention_years=25,
            accesses_per_year=1.0,
            mutability=0.2,
            retrieval_urgency=0.4,
            durability_priority=0.8,
            energy_priority=0.6,
            redundancy_priority=0.7,
            cost_priority=0.5,
            data_size_gb=10,
            expected_access_probability=0.3,
        ),
        lifecycle=_lifecycle(),
    )
    assert result.score_breakdowns is not None
    for tier, breakdown in result.score_breakdowns.items():
        assert result.scores[tier] == pytest.approx(breakdown.final_score)
    assert result.lifecycle_estimates is not None

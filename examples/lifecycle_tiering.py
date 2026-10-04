"""Evaluate lifecycle outputs from explicitly supplied research assumptions."""

from oligoark import LifecycleAssumptions, WorkloadProfile, recommend_storage_tier

lifecycle = LifecycleAssumptions.from_mapping(
    {
        "ssd": {
            "storage_cost_per_gb_year": 2.0,
            "retrieval_cost_per_gb": 0.0,
            "idle_energy_kwh_per_tb_year": 20.0,
            "retrieval_energy_kwh_per_gb": 0.001,
            "retrieval_latency_hours": 0.001,
        },
        "object_archive": {
            "storage_cost_per_gb_year": 0.5,
            "retrieval_cost_per_gb": 0.05,
            "idle_energy_kwh_per_tb_year": 5.0,
            "retrieval_energy_kwh_per_gb": 0.002,
            "retrieval_latency_hours": 1.0,
        },
        "tape": {
            "storage_cost_per_gb_year": 0.2,
            "retrieval_cost_per_gb": 0.02,
            "idle_energy_kwh_per_tb_year": 0.5,
            "retrieval_energy_kwh_per_gb": 0.005,
            "retrieval_latency_hours": 4.0,
        },
        "dna_future": {
            "storage_cost_per_gb_year": 1.0,
            "retrieval_cost_per_gb": 1.0,
            "idle_energy_kwh_per_tb_year": 0.1,
            "retrieval_energy_kwh_per_gb": 0.1,
            "retrieval_latency_hours": 24.0,
        },
    }
)
profile = WorkloadProfile(
    retention_years=10,
    accesses_per_year=2,
    mutability=0.1,
    retrieval_urgency=0.2,
    durability_priority=0.8,
    energy_priority=0.8,
    data_size_gb=10,
    expected_access_probability=0.5,
)
result = recommend_storage_tier(profile, lifecycle=lifecycle)
assert result.lifecycle_estimates is not None
print(result.to_dict())

"""Supply explicit normalized economics without treating them as vendor price claims."""

from oligoark.tiering import EconomicAssumptions, WorkloadProfile, recommend_storage_tier

profile = WorkloadProfile(100, 0.1, 0.0, 0.1, 1.0, 0.8, 0.8, 0.8)
economics = EconomicAssumptions(
    storage_cost_index={"ssd": 0.9, "object_archive": 0.3, "tape": 0.2, "dna_future": 0.5},
    retrieval_cost_index={"ssd": 0.1, "object_archive": 0.4, "tape": 0.7, "dna_future": 0.9},
)
recommendation = recommend_storage_tier(profile, economics)
print(recommendation.to_dict())

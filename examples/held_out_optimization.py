"""Calibrate an OligoArk codec plan and evaluate it on unseen seeds."""

from oligoark import (
    ChannelProfile,
    CodecSearchSpace,
    WorkloadProfile,
    evaluate_optimized_archive_plan,
)

payload = b"held-out OligoArk example payload"
workload = WorkloadProfile(
    retention_years=50,
    accesses_per_year=0.2,
    mutability=0.0,
    retrieval_urgency=0.2,
    durability_priority=0.9,
    energy_priority=0.5,
)
channel = ChannelProfile(substitution_rate=0.002, dropout_rate=0.01)
search = CodecSearchSpace(max_candidates=6, search_method="balanced", search_seed=5050)

result = evaluate_optimized_archive_plan(
    payload,
    workload,
    channel,
    calibration_seeds=(9001, 9002),
    evaluation_seeds=(2026, 2027),
    search_space=search,
)
print(result.to_dict())

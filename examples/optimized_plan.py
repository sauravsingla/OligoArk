"""Run a small measured archival optimisation with deterministic seeds."""

from oligoark import (
    ChannelProfile,
    CodecSearchSpace,
    WorkloadProfile,
    optimize_archive_plan,
)

payload = b"measured optimisation example" * 3
workload = WorkloadProfile(
    retention_years=100,
    accesses_per_year=0.1,
    mutability=0.0,
    retrieval_urgency=0.2,
    durability_priority=1.0,
    energy_priority=0.5,
    redundancy_priority=0.8,
    cost_priority=0.4,
)
channel = ChannelProfile(substitution_rate=0.003, dropout_rate=0.02)
search = CodecSearchSpace(
    chunk_sizes=(48, 64),
    rs_nsyms=(8, 16),
    redundancy_schemes=("xor", "hybrid"),
    parity_group_sizes=(3,),
    fountain_redundancies=(0.5,),
    reconstruction_modes=("direct", "graph"),
    max_candidates=6,
)
plan = optimize_archive_plan(
    payload,
    workload,
    channel,
    search_space=search,
    seeds=(2026,),
)
assert plan.optimization.evaluations
print(plan.to_dict())

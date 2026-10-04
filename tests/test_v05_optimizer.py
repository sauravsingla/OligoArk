from dataclasses import replace

import pytest

from oligoark.intelligence import evaluate_optimized_archive_plan
from oligoark.optimizer import (
    CodecSearchSpace,
    LifecycleObjectiveInputs,
    OptimizationWeights,
    optimize_codec,
    select_candidate_specs,
)
from oligoark.policy import ChannelProfile
from oligoark.tiering import WorkloadProfile


def _workload() -> WorkloadProfile:
    return WorkloadProfile(
        retention_years=50,
        accesses_per_year=0.2,
        mutability=0.0,
        retrieval_urgency=0.3,
        durability_priority=0.9,
        energy_priority=0.5,
    )


def _space() -> CodecSearchSpace:
    return CodecSearchSpace(
        chunk_sizes=(32, 48),
        rs_nsyms=(8, 12),
        redundancy_schemes=("xor", "fountain", "hybrid"),
        parity_group_sizes=(3, 5),
        fountain_redundancies=(0.25, 0.5),
        reconstruction_modes=("direct", "graph"),
        max_candidates=8,
        search_method="balanced",
        search_seed=77,
    )


def test_balanced_candidate_sampling_is_order_independent() -> None:
    left = _space()
    right = replace(
        left,
        chunk_sizes=tuple(reversed(left.chunk_sizes)),
        rs_nsyms=tuple(reversed(left.rs_nsyms)),
        redundancy_schemes=tuple(reversed(left.redundancy_schemes)),
        parity_group_sizes=tuple(reversed(left.parity_group_sizes)),
        fountain_redundancies=tuple(reversed(left.fountain_redundancies)),
        reconstruction_modes=tuple(reversed(left.reconstruction_modes)),
        constraints=tuple(reversed(left.constraints)),
    )
    assert [spec.key() for spec in select_candidate_specs(left)] == [
        spec.key() for spec in select_candidate_specs(right)
    ]


def test_full_grid_ignores_candidate_budget() -> None:
    space = replace(_space(), max_candidates=1, search_method="full_grid")
    assert len(select_candidate_specs(space)) > 1


def test_objective_breakdown_reproduces_score_and_lifecycle_is_explicit() -> None:
    result = optimize_codec(
        b"v05 optimizer objective" * 3,
        ChannelProfile(substitution_rate=0.002),
        _workload(),
        search_space=_space(),
        weights=OptimizationWeights(),
        seeds=(7001,),
        lifecycle=LifecycleObjectiveInputs(
            storage_cost=12.0,
            retrieval_cost=3.0,
            energy_kwh=1.5,
            retrieval_latency_hours=2.0,
        ),
    )
    winner = max(
        (item for item in result.evaluations if item.rejected_reason is None),
        key=lambda item: item.score,
    )
    assert winner.score == pytest.approx(winner.objective.total)
    assert winner.objective.lifecycle_storage_cost_penalty >= 0
    assert winner.redundancy_ratio > 0
    assert result.total_possible_candidates >= result.evaluated_candidates
    assert result.search_method == "balanced"


def test_held_out_optimizer_rejects_seed_leakage() -> None:
    with pytest.raises(ValueError, match="disjoint"):
        evaluate_optimized_archive_plan(
            b"held-out",
            _workload(),
            ChannelProfile(),
            calibration_seeds=(1, 2),
            evaluation_seeds=(2, 3),
        )


def test_held_out_optimizer_evaluates_unseen_seed() -> None:
    result = evaluate_optimized_archive_plan(
        b"held-out optimizer" * 2,
        _workload(),
        ChannelProfile(substitution_rate=0.001),
        calibration_seeds=(8101,),
        evaluation_seeds=(8201,),
        search_space=replace(_space(), max_candidates=4),
    )
    assert result.calibration_seeds == (8101,)
    assert result.evaluation_seeds == (8201,)
    assert len(result.trials) == 1
    assert 0 <= result.held_out_recovery_rate <= 1

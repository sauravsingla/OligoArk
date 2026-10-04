from oligoark.intelligence import optimize_archive_plan
from oligoark.optimizer import CodecSearchSpace, OptimizationWeights, optimize_codec
from oligoark.policy import ChannelProfile
from oligoark.tiering import WorkloadProfile


def _workload() -> WorkloadProfile:
    return WorkloadProfile(
        100,
        0.1,
        0.0,
        0.2,
        1.0,
        0.5,
        redundancy_priority=0.8,
        cost_priority=0.4,
        data_size_gb=0.001,
    )


def _small_space() -> CodecSearchSpace:
    return CodecSearchSpace(
        chunk_sizes=(32, 48),
        rs_nsyms=(8, 16),
        redundancy_schemes=("xor", "fountain"),
        parity_group_sizes=(3,),
        fountain_redundancies=(0.75,),
        reconstruction_modes=("direct", "graph"),
        max_candidates=8,
        mask_search_limit=96,
    )


def test_optimizer_actually_evaluates_candidates() -> None:
    result = optimize_codec(
        b"optimizer research payload" * 4,
        ChannelProfile(substitution_rate=0.005, dropout_rate=0.02),
        _workload(),
        search_space=_small_space(),
        weights=OptimizationWeights(),
        seeds=(2026,),
    )
    assert len(result.evaluations) > 1
    assert result.best_config.redundancy_scheme in {"xor", "fountain"}
    assert result.reconstruction_mode in {"direct", "graph"}
    assert any(evaluation.encoded_nucleotides > 0 for evaluation in result.evaluations)
    assert "evaluated" in result.rationale[0]


def test_integrated_optimized_plan_exposes_selected_strategy() -> None:
    result = optimize_archive_plan(
        b"integrated plan payload" * 3,
        _workload(),
        ChannelProfile(substitution_rate=0.003),
        search_space=_small_space(),
        seeds=(2026,),
    )
    assert result.selected_redundancy_scheme in {"xor", "fountain"}
    assert result.selected_reconstruction_strategy in {"direct", "graph"}
    assert "min_gc_fraction" in result.selected_sequence_constraints
    assert result.optimization.evaluations

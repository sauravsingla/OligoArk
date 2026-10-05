from dataclasses import replace

import pytest

import oligoark.optimizer as optimizer_module
from oligoark.optimizer import CodecSearchSpace, OptimizationWeights, optimize_codec
from oligoark.policy import ChannelProfile
from oligoark.tiering import WorkloadProfile


def _workload() -> WorkloadProfile:
    return WorkloadProfile(
        retention_years=80,
        accesses_per_year=0.1,
        mutability=0.0,
        retrieval_urgency=0.2,
        durability_priority=1.0,
        energy_priority=0.5,
    )


def test_balanced_robust_optimizer_tracks_cross_seed_instability() -> None:
    payload = b"robust optimizer calibration" * 2
    search = CodecSearchSpace(
        chunk_sizes=(32, 48),
        rs_nsyms=(8, 16),
        redundancy_schemes=("xor", "hybrid"),
        parity_group_sizes=(3,),
        fountain_redundancies=(0.35,),
        reconstruction_modes=("direct", "trace"),
        max_candidates=6,
        search_method="balanced_robust",
        search_seed=6060,
    )
    result = optimize_codec(
        payload,
        ChannelProfile(insertion_rate=0.0005, deletion_rate=0.0005),
        _workload(),
        search_space=search,
        weights=OptimizationWeights(instability=0.2),
        seeds=(9401, 9402, 9403, 9404),
        duplicate_rate=0.5,
        calibration_payloads=(payload, bytes(reversed(payload))),
    )
    assert result.search_method == "balanced_robust"
    assert result.evaluated_candidates == 6
    assert all(
        0 <= evaluation.recovery_instability <= 1
        for evaluation in result.evaluations
        if evaluation.rejected_reason is None
    )
    assert any(
        evaluation.reconstruction_mode == "trace"
        for evaluation in result.evaluations
    )


def test_calibration_payload_variants_must_match_size() -> None:
    search = replace(
        CodecSearchSpace(),
        reconstruction_modes=("direct",),
        max_candidates=1,
        search_method="balanced_robust",
    )
    try:
        optimize_codec(
            b"same-size",
            ChannelProfile(),
            _workload(),
            search_space=search,
            seeds=(1,),
            calibration_payloads=(b"different-size",),
        )
    except ValueError as exc:
        assert "match the primary payload length" in str(exc)
    else:
        raise AssertionError("mismatched calibration payloads should fail")


def test_zero_runtime_weight_makes_tie_break_wall_clock_independent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    search = CodecSearchSpace(
        chunk_sizes=(48,),
        rs_nsyms=(8,),
        redundancy_schemes=("xor",),
        parity_group_sizes=(3,),
        fountain_redundancies=(0.25,),
        reconstruction_modes=("direct", "graph"),
        max_candidates=8,
        search_method="full_grid",
    )
    weights = OptimizationWeights(runtime=0.0, retrieval=0.0)

    def run_with_times(direct: float, graph: float):
        def fake_simulate(
            payload: bytes,
            spec: optimizer_module.CandidateSpec,
            channel: ChannelProfile,
            *,
            seed: int,
            duplicate_rate: float | None,
            copies_per_strand: int,
        ) -> tuple[bool, int, float, bool]:
            del payload, channel, seed, duplicate_rate, copies_per_strand
            elapsed = direct if spec.reconstruction_mode == "direct" else graph
            return True, 100, elapsed, False

        monkeypatch.setattr(optimizer_module, "_simulate_once", fake_simulate)
        return optimize_codec(
            b"runtime-independent-tie",
            ChannelProfile(),
            _workload(),
            search_space=search,
            weights=weights,
            seeds=(9401, 9402),
        )

    slow_direct = run_with_times(10.0, 1.0)
    fast_direct = run_with_times(1.0, 10.0)
    assert slow_direct.best_config == fast_direct.best_config
    assert slow_direct.reconstruction_mode == fast_direct.reconstruction_mode
    assert slow_direct.best_score == pytest.approx(fast_direct.best_score)

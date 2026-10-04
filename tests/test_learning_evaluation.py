import pytest

from oligoark.experiments import ExperimentRecord
from oligoark.learning import LinearUtilityPolicyModel, PolicyObservation
from oligoark.learning_eval import (
    evaluate_learning_from_records,
    policy_observations_from_records,
)
from oligoark.policy import ChannelProfile, CodecPolicy


def _record(
    strategy: str,
    seed: int,
    recovered: bool,
    chunk: int,
    rs: int,
    parity: int,
    overhead: float,
    runtime: float,
    scenario: str = "sub",
) -> ExperimentRecord:
    return ExperimentRecord(
        strategy=strategy,
        scenario=scenario,
        seed=seed,
        payload_size=128,
        payload_sha256="0" * 64,
        recovered=recovered,
        encoded_nucleotides=int(overhead * 512),
        overhead_ratio=overhead,
        strand_count=2,
        read_count=2,
        runtime_seconds=runtime,
        graph_reconstruction_used=False,
        redundancy_scheme="xor",
        rs_nsym=rs,
        chunk_size=chunk,
        parity_group_size=parity,
        adaptive_masks=True,
        substitution_rate=0.01,
        insertion_rate=0.0,
        deletion_rate=0.0,
        dropout_rate=0.0,
        duplicate_rate=0.0,
        calibration_seeds=(),
    )


def test_experiment_records_convert_to_learning_dataset() -> None:
    records = [
        _record("fixed", 1, False, 96, 8, 8, 1.5, 0.1),
        _record("adaptive", 1, True, 64, 16, 5, 1.8, 0.2),
    ]
    observations = policy_observations_from_records(records)
    assert len(observations) == 2
    assert observations[1].recovered is True


def test_linear_model_serialization_round_trip() -> None:
    lean = CodecPolicy(96, 8, 8, True, ("lean",))
    robust = CodecPolicy(64, 16, 5, True, ("robust",))
    observations = [
        PolicyObservation(ChannelProfile(), lean, True, 100, 0.1),
        PolicyObservation(ChannelProfile(0.01, 0, 0, 0), robust, True, 140, 0.2),
        PolicyObservation(ChannelProfile(0.02, 0, 0, 0), lean, False, 100, 0.1),
    ]
    model = LinearUtilityPolicyModel(0.01).fit(observations)
    restored = LinearUtilityPolicyModel.from_dict(model.to_dict())
    target = ChannelProfile(0.015, 0, 0, 0)
    assert restored.recommend(target).policy == model.recommend(target).policy
    assert restored.coefficients == pytest.approx(model.coefficients)


def test_held_out_learning_pipeline_reports_regret_and_baselines() -> None:
    records = []
    for seed in (1, 2, 3, 4):
        records.extend(
            [
                _record("fixed", seed, False, 96, 8, 8, 1.5, 0.1),
                _record("adaptive", seed, True, 64, 16, 8, 1.8, 0.2),
                _record("combined", seed, True, 64, 16, 8, 1.7, 0.18),
            ]
        )
    result = evaluate_learning_from_records(
        records,
        training_seeds=(1, 2),
        test_seeds=(3, 4),
    )
    methods = {summary.method for summary in result.summaries}
    assert methods == {"heuristic", "empirical", "linear", "measured_search"}
    assert result.training_seeds == (1, 2)
    assert result.test_seeds == (3, 4)
    assert result.test_groups == 2


def test_learning_split_rejects_leakage() -> None:
    records = [
        _record("fixed", 1, True, 96, 8, 8, 1.5, 0.1),
        _record("adaptive", 1, True, 64, 16, 5, 1.8, 0.2),
    ]
    with pytest.raises(ValueError, match="disjoint"):
        evaluate_learning_from_records(records, training_seeds=(1,), test_seeds=(1,))


def test_learning_rejects_insufficient_and_malformed_model_state() -> None:
    policy = CodecPolicy(96, 8, 8, True, ("single",))
    observation = PolicyObservation(ChannelProfile(), policy, True, 100, 0.1)
    with pytest.raises(ValueError, match="at least two"):
        LinearUtilityPolicyModel().fit([observation])
    with pytest.raises(ValueError, match="Unsupported serialized"):
        LinearUtilityPolicyModel.from_dict({"model": "unknown"})


def test_learning_can_hold_out_unseen_channel_regimes() -> None:
    records = [
        _record("fixed", 1, True, 96, 8, 8, 1.5, 0.1, scenario="train"),
        _record("adaptive", 1, True, 96, 8, 8, 1.6, 0.11, scenario="train"),
        _record("fixed", 2, False, 96, 8, 8, 1.5, 0.1, scenario="unseen"),
        _record("adaptive", 2, True, 64, 16, 8, 1.9, 0.2, scenario="unseen"),
        _record("combined", 2, True, 64, 16, 8, 1.8, 0.18, scenario="unseen"),
    ]
    result = evaluate_learning_from_records(
        records,
        training_seeds=(1,),
        test_seeds=(2,),
        training_scenarios=("train",),
        test_scenarios=("unseen",),
    )
    assert result.training_scenarios == ("train",)
    assert result.test_scenarios == ("unseen",)
    assert result.test_groups == 1

    with pytest.raises(ValueError, match="scenarios must be disjoint"):
        evaluate_learning_from_records(
            records,
            training_seeds=(1,),
            test_seeds=(2,),
            training_scenarios=("train",),
            test_scenarios=("train",),
        )

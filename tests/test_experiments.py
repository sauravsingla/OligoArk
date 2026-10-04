from oligoark.experiments import (
    ExperimentProfile,
    ExperimentScenario,
    aggregate_experiments,
    run_experiments,
    wilson_interval,
)
from oligoark.policy import ChannelProfile


def test_wilson_interval_is_bounded() -> None:
    low, high = wilson_interval(7, 10)
    assert 0 <= low <= 0.7 <= high <= 1


def test_small_ablation_experiment_runs_all_strategies() -> None:
    profile = ExperimentProfile(
        seeds=(2026,),
        payload_sizes=(96,),
        scenarios=(
            ExperimentScenario(
                "substitution",
                ChannelProfile(substitution_rate=0.002),
                duplicate_rate=0.2,
            ),
        ),
        optimizer_max_candidates=4,
    )
    records = run_experiments(profile)
    assert {record.strategy for record in records} == set(profile.strategies)
    assert all(record.encoded_nucleotides > 0 for record in records)
    summaries = aggregate_experiments(records)
    assert len(summaries) == len(profile.strategies)
    assert all(summary.trials == 1 for summary in summaries)

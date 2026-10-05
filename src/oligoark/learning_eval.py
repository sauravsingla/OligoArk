"""Held-out evaluation pipeline for transparent policy-learning baselines."""

from __future__ import annotations

import statistics
from collections.abc import Callable
from dataclasses import asdict, dataclass

from .experiments import ExperimentRecord
from .learning import (
    EmpiricalPolicyModel,
    KernelUtilityPolicyModel,
    LinearUtilityPolicyModel,
    PolicyObservation,
)
from .policy import ChannelProfile, CodecPolicy, recommend_codec_policy


@dataclass(frozen=True)
class LearningMethodSummary:
    method: str
    trials: int
    recovery_rate: float
    mean_overhead_ratio: float
    mean_runtime_seconds: float
    selection_accuracy: float
    mean_regret: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class LearningEvaluationResult:
    training_seeds: tuple[int, ...]
    validation_seeds: tuple[int, ...]
    test_seeds: tuple[int, ...]
    training_scenarios: tuple[str, ...]
    test_scenarios: tuple[str, ...]
    training_observations: int
    validation_observations: int
    test_groups: int
    selected_ridge: float
    selected_kernel_bandwidth: float
    validation_regret: dict[str, float]
    summaries: tuple[LearningMethodSummary, ...]
    linear_model_state: dict[str, object]
    kernel_model_state: dict[str, object]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _policy(record: ExperimentRecord) -> CodecPolicy:
    return CodecPolicy(
        chunk_size=record.chunk_size,
        rs_nsym=record.rs_nsym,
        parity_group_size=record.parity_group_size,
        adaptive_masks=record.adaptive_masks,
        rationale=("derived from reproducible experiment record",),
    )


def _policy_key(policy: CodecPolicy) -> tuple[int, int, int, bool]:
    return (
        policy.chunk_size,
        policy.rs_nsym,
        policy.parity_group_size,
        policy.adaptive_masks,
    )


def _channel(record: ExperimentRecord) -> ChannelProfile:
    return ChannelProfile(
        record.substitution_rate,
        record.insertion_rate,
        record.deletion_rate,
        record.dropout_rate,
    )


def policy_observations_from_records(
    records: list[ExperimentRecord],
    *,
    allowed_strategies: tuple[str, ...] = ("fixed", "adaptive"),
) -> list[PolicyObservation]:
    """Convert experiment output into learning observations without hidden data."""
    observations: list[PolicyObservation] = []
    for record in records:
        if record.strategy not in allowed_strategies:
            continue
        observations.append(
            PolicyObservation(
                channel=_channel(record),
                policy=_policy(record),
                recovered=record.recovered,
                encoded_nucleotides=record.encoded_nucleotides,
                runtime_seconds=record.runtime_seconds,
            )
        )
    if not observations:
        raise ValueError("No experiment records matched the allowed learning strategies")
    return observations


def _utility(record: ExperimentRecord, max_overhead: float, max_runtime: float) -> float:
    recovery = 1.0 if record.recovered else 0.0
    overhead = record.overhead_ratio / max(1e-12, max_overhead)
    runtime = record.runtime_seconds / max(1e-12, max_runtime)
    return recovery - 0.20 * overhead - 0.05 * runtime


def _match_record(
    records: list[ExperimentRecord],
    policy: CodecPolicy,
) -> ExperimentRecord:
    key = _policy_key(policy)
    matches = [record for record in records if _policy_key(_policy(record)) == key]
    if not matches:
        raise ValueError("Recommended policy is unavailable in held-out candidate records")
    return max(
        matches,
        key=lambda record: (
            int(record.recovered),
            -record.overhead_ratio,
            -record.runtime_seconds,
        ),
    )


def _selector_validation_regret(
    records: list[ExperimentRecord],
    selector: Callable[[ChannelProfile, list[CodecPolicy]], CodecPolicy],
) -> float:
    grouped: dict[tuple[str, int, int], list[ExperimentRecord]] = {}
    for record in records:
        grouped.setdefault((record.scenario, record.payload_size, record.seed), []).append(record)

    regrets: list[float] = []
    for group in grouped.values():
        codec_candidates = [
            record for record in group if record.strategy in {"fixed", "adaptive"}
        ]
        if len(codec_candidates) < 2:
            continue
        unique_policies: dict[tuple[int, int, int, bool], CodecPolicy] = {}
        for record in codec_candidates:
            policy = _policy(record)
            unique_policies[_policy_key(policy)] = policy
        candidates = list(unique_policies.values())
        channel = _channel(codec_candidates[0])
        max_overhead = max(record.overhead_ratio for record in codec_candidates)
        max_runtime = max(record.runtime_seconds for record in codec_candidates)
        oracle_utility = max(
            _utility(record, max_overhead, max_runtime)
            for record in codec_candidates
        )
        selected = _match_record(codec_candidates, selector(channel, candidates))
        selected_utility = _utility(selected, max_overhead, max_runtime)
        regrets.append(max(0.0, oracle_utility - selected_utility))
    if not regrets:
        raise ValueError("No complete validation groups were available for model selection")
    return statistics.fmean(regrets)


def evaluate_learning_from_records(
    records: list[ExperimentRecord],
    *,
    training_seeds: tuple[int, ...],
    test_seeds: tuple[int, ...],
    validation_seeds: tuple[int, ...] = (),
    training_scenarios: tuple[str, ...] | None = None,
    test_scenarios: tuple[str, ...] | None = None,
) -> LearningEvaluationResult:
    """Evaluate policy learning on disjoint seeds and, optionally, disjoint channels."""
    if not training_seeds or not test_seeds:
        raise ValueError("training and test seeds must not be empty")
    if set(training_seeds) & set(test_seeds):
        raise ValueError("training and test seeds must be disjoint")
    if validation_seeds and (
        set(validation_seeds) & set(training_seeds)
        or set(validation_seeds) & set(test_seeds)
    ):
        raise ValueError("training, validation, and test seeds must be disjoint")

    all_scenarios = tuple(sorted({record.scenario for record in records}))
    resolved_training_scenarios = training_scenarios or all_scenarios
    resolved_test_scenarios = test_scenarios or all_scenarios
    if (
        training_scenarios is not None
        and test_scenarios is not None
        and set(resolved_training_scenarios) & set(resolved_test_scenarios)
    ):
        raise ValueError("training and test scenarios must be disjoint")

    training_records = [
        record
        for record in records
        if record.seed in training_seeds
        and record.scenario in resolved_training_scenarios
    ]
    validation_records = [
        record
        for record in records
        if record.seed in validation_seeds
        and record.scenario in resolved_training_scenarios
    ]
    test_records = [
        record
        for record in records
        if record.seed in test_seeds
        and record.scenario in resolved_test_scenarios
    ]
    if not training_records or not test_records:
        raise ValueError(
            "records do not cover both requested training and test seed/scenario sets"
        )
    if validation_seeds and not validation_records:
        raise ValueError("records do not cover the requested validation seeds")

    training_observations = policy_observations_from_records(training_records)
    selected_ridge = 0.01
    selected_bandwidth = 0.01
    validation_regret: dict[str, float] = {}
    if validation_records:
        ridge_candidates = (0.001, 0.01, 0.1)
        ridge_scores: list[tuple[float, float]] = []
        for ridge in ridge_candidates:
            linear_candidate = LinearUtilityPolicyModel(ridge=ridge).fit(
                training_observations
            )

            def linear_selector(
                channel: ChannelProfile,
                candidates: list[CodecPolicy],
                *,
                model: LinearUtilityPolicyModel = linear_candidate,
            ) -> CodecPolicy:
                return model.recommend(channel, candidates=candidates).policy

            regret = _selector_validation_regret(
                validation_records,
                linear_selector,
            )
            ridge_scores.append((regret, ridge))
        selected_ridge = min(ridge_scores)[1]
        validation_regret["linear"] = round(min(ridge_scores)[0], 6)

        bandwidth_candidates = (0.005, 0.01, 0.02, 0.05)
        bandwidth_scores: list[tuple[float, float]] = []
        for bandwidth in bandwidth_candidates:
            kernel_candidate = KernelUtilityPolicyModel(bandwidth=bandwidth).fit(
                training_observations
            )

            def kernel_selector(
                channel: ChannelProfile,
                candidates: list[CodecPolicy],
                *,
                model: KernelUtilityPolicyModel = kernel_candidate,
            ) -> CodecPolicy:
                return model.recommend(channel, candidates=candidates).policy

            regret = _selector_validation_regret(
                validation_records,
                kernel_selector,
            )
            bandwidth_scores.append((regret, bandwidth))
        selected_bandwidth = min(bandwidth_scores)[1]
        validation_regret["kernel"] = round(min(bandwidth_scores)[0], 6)

    validation_observations = (
        policy_observations_from_records(validation_records)
        if validation_records
        else []
    )
    fitted_observations = training_observations + validation_observations
    empirical = EmpiricalPolicyModel().fit(fitted_observations)
    linear = LinearUtilityPolicyModel(ridge=selected_ridge).fit(fitted_observations)
    kernel = KernelUtilityPolicyModel(bandwidth=selected_bandwidth).fit(
        fitted_observations
    )

    grouped: dict[tuple[str, int, int], list[ExperimentRecord]] = {}
    for record in test_records:
        grouped.setdefault((record.scenario, record.payload_size, record.seed), []).append(record)

    metrics: dict[str, list[tuple[ExperimentRecord, float, bool]]] = {
        "heuristic": [],
        "empirical": [],
        "linear": [],
        "kernel": [],
        "adaptive_fountain": [],
        "measured_search": [],
    }

    for group in grouped.values():
        codec_candidates = [
            record for record in group if record.strategy in {"fixed", "adaptive"}
        ]
        fountain_candidates = [
            record for record in group if record.strategy == "adaptive_fountain"
        ]
        search_candidates = [record for record in group if record.strategy == "combined"]
        if len(codec_candidates) < 2 or not fountain_candidates or not search_candidates:
            continue

        unique_policies: dict[tuple[int, int, int, bool], CodecPolicy] = {}
        for record in codec_candidates:
            policy = _policy(record)
            unique_policies[_policy_key(policy)] = policy
        candidates = list(unique_policies.values())
        channel = _channel(codec_candidates[0])

        max_overhead = max(record.overhead_ratio for record in codec_candidates)
        max_runtime = max(record.runtime_seconds for record in codec_candidates)
        oracle_record = max(
            codec_candidates,
            key=lambda record: _utility(record, max_overhead, max_runtime),
        )
        oracle_utility = _utility(oracle_record, max_overhead, max_runtime)
        oracle_key = _policy_key(_policy(oracle_record))

        heuristic_policy = recommend_codec_policy(channel)
        heuristic_record = _match_record(codec_candidates, heuristic_policy)
        empirical_record = _match_record(
            codec_candidates,
            empirical.recommend(channel, candidates=candidates).policy,
        )
        linear_record = _match_record(
            codec_candidates,
            linear.recommend(channel, candidates=candidates).policy,
        )
        kernel_record = _match_record(
            codec_candidates,
            kernel.recommend(channel, candidates=candidates).policy,
        )
        fountain_record = max(
            fountain_candidates,
            key=lambda record: (
                int(record.recovered),
                -record.overhead_ratio,
                -record.runtime_seconds,
            ),
        )
        search_record = max(
            search_candidates,
            key=lambda record: (
                int(record.recovered),
                -record.overhead_ratio,
                -record.runtime_seconds,
            ),
        )

        for method, selected in (
            ("heuristic", heuristic_record),
            ("empirical", empirical_record),
            ("linear", linear_record),
            ("kernel", kernel_record),
        ):
            selected_utility = _utility(selected, max_overhead, max_runtime)
            metrics[method].append(
                (
                    selected,
                    max(0.0, oracle_utility - selected_utility),
                    _policy_key(_policy(selected)) == oracle_key,
                )
            )

        fountain_utility = (
            (1.0 if fountain_record.recovered else 0.0)
            - 0.20 * fountain_record.overhead_ratio / max(1e-12, max_overhead)
            - 0.05 * fountain_record.runtime_seconds / max(1e-12, max_runtime)
        )
        metrics["adaptive_fountain"].append(
            (
                fountain_record,
                max(0.0, oracle_utility - fountain_utility),
                False,
            )
        )

        search_regret = max(
            0.0,
            oracle_utility
            - (
                (1.0 if search_record.recovered else 0.0)
                - 0.20 * search_record.overhead_ratio / max(1e-12, max_overhead)
                - 0.05 * search_record.runtime_seconds / max(1e-12, max_runtime)
            ),
        )
        metrics["measured_search"].append((search_record, search_regret, False))

    if not metrics["linear"] or not metrics["kernel"]:
        raise ValueError("No complete held-out groups were available for learning evaluation")

    summaries: list[LearningMethodSummary] = []
    for method in (
        "heuristic",
        "empirical",
        "linear",
        "kernel",
        "adaptive_fountain",
        "measured_search",
    ):
        values = metrics[method]
        summaries.append(
            LearningMethodSummary(
                method=method,
                trials=len(values),
                recovery_rate=round(
                    statistics.fmean(float(record.recovered) for record, _, _ in values), 6
                ),
                mean_overhead_ratio=round(
                    statistics.fmean(record.overhead_ratio for record, _, _ in values), 6
                ),
                mean_runtime_seconds=round(
                    statistics.fmean(record.runtime_seconds for record, _, _ in values), 6
                ),
                selection_accuracy=round(
                    statistics.fmean(float(correct) for _, _, correct in values), 6
                )
                if method not in {"adaptive_fountain", "measured_search"}
                else 0.0,
                mean_regret=round(statistics.fmean(regret for _, regret, _ in values), 6),
            )
        )

    return LearningEvaluationResult(
        training_seeds=tuple(training_seeds),
        validation_seeds=tuple(validation_seeds),
        test_seeds=tuple(test_seeds),
        training_scenarios=tuple(resolved_training_scenarios),
        test_scenarios=tuple(resolved_test_scenarios),
        training_observations=len(training_observations),
        validation_observations=len(validation_observations),
        test_groups=len(metrics["linear"]),
        selected_ridge=selected_ridge,
        selected_kernel_bandwidth=selected_bandwidth,
        validation_regret=validation_regret,
        summaries=tuple(summaries),
        linear_model_state=linear.to_dict(),
        kernel_model_state=kernel.to_dict(),
    )

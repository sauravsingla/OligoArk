"""Transparent deterministic policy-learning baselines for OligoArk."""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

from .policy import ChannelProfile, CodecPolicy


@dataclass(frozen=True)
class PolicyObservation:
    channel: ChannelProfile
    policy: CodecPolicy
    recovered: bool
    encoded_nucleotides: int
    runtime_seconds: float


@dataclass(frozen=True)
class LearnedPolicyRecommendation:
    policy: CodecPolicy
    confidence: float
    supporting_observations: int
    rationale: tuple[str, ...]


class EmpiricalPolicyModel:
    """Instance-based policy ranker learned from prior reproducible observations."""

    def __init__(self) -> None:
        self._observations: list[PolicyObservation] = []

    def fit(self, observations: list[PolicyObservation]) -> EmpiricalPolicyModel:
        if not observations:
            raise ValueError("At least one policy observation is required")
        self._observations = list(observations)
        return self

    @staticmethod
    def _distance(left: ChannelProfile, right: ChannelProfile) -> float:
        return (
            abs(left.substitution_rate - right.substitution_rate)
            + abs(left.insertion_rate - right.insertion_rate)
            + abs(left.deletion_rate - right.deletion_rate)
            + abs(left.dropout_rate - right.dropout_rate)
        )

    def recommend(
        self,
        channel: ChannelProfile,
        candidates: list[CodecPolicy] | None = None,
    ) -> LearnedPolicyRecommendation:
        if not self._observations:
            raise ValueError("Model must be fit before recommendation")

        allowed = {_policy_key(policy) for policy in candidates} if candidates else None
        max_nucleotides = max(obs.encoded_nucleotides for obs in self._observations)
        max_runtime = max(obs.runtime_seconds for obs in self._observations) or 1.0
        scores: dict[tuple[int, int, int, bool], float] = {}
        counts: dict[tuple[int, int, int, bool], int] = {}
        policies: dict[tuple[int, int, int, bool], CodecPolicy] = {}

        for obs in self._observations:
            key = _policy_key(obs.policy)
            if allowed is not None and key not in allowed:
                continue
            distance = self._distance(channel, obs.channel)
            weight = 1.0 / (0.001 + distance)
            recovery_score = 1.0 if obs.recovered else -1.0
            overhead_penalty = obs.encoded_nucleotides / max(1, max_nucleotides)
            runtime_penalty = obs.runtime_seconds / max_runtime
            value = recovery_score - 0.15 * overhead_penalty - 0.05 * runtime_penalty
            scores[key] = scores.get(key, 0.0) + weight * value
            counts[key] = counts.get(key, 0) + 1
            policies[key] = obs.policy

        if not scores:
            raise ValueError("No fitted observations match the supplied candidate policies")
        winner = max(scores, key=lambda policy_key: scores[policy_key])
        total_magnitude = sum(abs(value) for value in scores.values()) or 1.0
        confidence = min(1.0, abs(scores[winner]) / total_magnitude)
        rationale = (
            "inverse-distance weighting over prior channel/policy observations",
            "recovery success dominates; encoded overhead and runtime act as penalties",
            "confidence is a relative score, not a calibrated probability",
        )
        return LearnedPolicyRecommendation(
            policy=policies[winner],
            confidence=round(confidence, 4),
            supporting_observations=counts[winner],
            rationale=rationale,
        )


@dataclass(frozen=True)
class LinearPolicyRecommendation:
    policy: CodecPolicy
    predicted_utility: float
    margin: float
    candidate_count: int
    rationale: tuple[str, ...]


def _policy_key(policy: CodecPolicy) -> tuple[int, int, int, bool]:
    return (
        policy.chunk_size,
        policy.rs_nsym,
        policy.parity_group_size,
        policy.adaptive_masks,
    )


def _features(channel: ChannelProfile, policy: CodecPolicy) -> list[float]:
    substitution = channel.substitution_rate
    insertion = channel.insertion_rate
    deletion = channel.deletion_rate
    dropout = channel.dropout_rate
    chunk = policy.chunk_size / 128.0
    rs = policy.rs_nsym / 64.0
    parity = policy.parity_group_size / 16.0
    masks = 1.0 if policy.adaptive_masks else 0.0
    return [
        1.0,
        substitution,
        insertion,
        deletion,
        dropout,
        chunk,
        rs,
        parity,
        masks,
        substitution * rs,
        (insertion + deletion) * chunk,
        dropout / max(0.0625, parity),
        (substitution + insertion + deletion) * masks,
    ]


def _solve_linear_system(matrix: list[list[float]], vector: list[float]) -> list[float]:
    size = len(vector)
    augmented = [row[:] + [vector[index]] for index, row in enumerate(matrix)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-12:
            raise ValueError("policy-learning linear system is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == 0:
                continue
            augmented[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(
                    augmented[row],
                    augmented[column],
                    strict=True,
                )
            ]
    return [augmented[row][-1] for row in range(size)]


class KernelUtilityPolicyModel:
    """Deterministic RBF-kernel utility baseline over reproducible observations.

    The model stores no proprietary state and performs a transparent weighted average
    over observed channel/policy utilities. It is stronger than nearest-neighbor
    ranking because all observations contribute smoothly to each candidate estimate.
    """

    def __init__(self, bandwidth: float = 0.02) -> None:
        if bandwidth <= 0:
            raise ValueError("bandwidth must be positive")
        self.bandwidth = bandwidth
        self._observations: list[PolicyObservation] = []
        self._max_nucleotides = 1.0
        self._max_runtime = 1.0

    def fit(self, observations: list[PolicyObservation]) -> KernelUtilityPolicyModel:
        if len(observations) < 2:
            raise ValueError("KernelUtilityPolicyModel requires at least two observations")
        self._observations = list(observations)
        self._max_nucleotides = float(
            max(max(1, observation.encoded_nucleotides) for observation in observations)
        )
        self._max_runtime = max(
            1e-9,
            max(observation.runtime_seconds for observation in observations),
        )
        return self

    @staticmethod
    def _channel_vector(channel: ChannelProfile) -> tuple[float, float, float, float]:
        return (
            channel.substitution_rate,
            channel.insertion_rate,
            channel.deletion_rate,
            channel.dropout_rate,
        )

    def _utility(self, observation: PolicyObservation) -> float:
        recovery = 1.0 if observation.recovered else 0.0
        overhead = observation.encoded_nucleotides / self._max_nucleotides
        runtime = observation.runtime_seconds / self._max_runtime
        return recovery - 0.20 * overhead - 0.05 * runtime

    def predict_utility(self, channel: ChannelProfile, policy: CodecPolicy) -> float:
        if not self._observations:
            raise ValueError("Model must be fit before prediction")
        channel.validate()
        target = self._channel_vector(channel)
        key = _policy_key(policy)
        numerator = 0.0
        denominator = 0.0
        for observation in self._observations:
            if _policy_key(observation.policy) != key:
                continue
            source = self._channel_vector(observation.channel)
            distance_sq = sum(
                (left - right) ** 2 for left, right in zip(target, source, strict=True)
            )
            weight = math.exp(-distance_sq / (2.0 * self.bandwidth**2))
            numerator += weight * self._utility(observation)
            denominator += weight
        if denominator <= 1e-15:
            return -1_000_000.0
        return numerator / denominator

    def recommend(
        self,
        channel: ChannelProfile,
        candidates: list[CodecPolicy] | None = None,
    ) -> LinearPolicyRecommendation:
        if not self._observations:
            raise ValueError("Model must be fit before recommendation")
        resolved = candidates or list(
            {
                _policy_key(observation.policy): observation.policy
                for observation in self._observations
            }.values()
        )
        if not resolved:
            raise ValueError("At least one candidate policy is required")
        scored = sorted(
            (
                (self.predict_utility(channel, policy), policy)
                for policy in resolved
            ),
            key=lambda item: (
                item[0],
                item[1].rs_nsym,
                -item[1].chunk_size,
            ),
            reverse=True,
        )
        best_score, best_policy = scored[0]
        second_score = scored[1][0] if len(scored) > 1 else best_score
        margin = max(0.0, best_score - second_score)
        return LinearPolicyRecommendation(
            policy=best_policy,
            predicted_utility=round(best_score, 6),
            margin=round(1.0 - math.exp(-margin), 6),
            candidate_count=len(resolved),
            rationale=(
                "deterministic Gaussian-kernel utility regression",
                "all matching-policy observations contribute by channel similarity",
                "margin is relative separation, not a calibrated probability",
            ),
        )

    def to_dict(self) -> dict[str, object]:
        if not self._observations:
            raise ValueError("Model must be fit before serialization")
        return {
            "model": "rbf-kernel-utility",
            "bandwidth": self.bandwidth,
            "max_nucleotides": self._max_nucleotides,
            "max_runtime": self._max_runtime,
            "observations": [
                {
                    "channel": asdict(observation.channel),
                    "policy": asdict(observation.policy),
                    "recovered": observation.recovered,
                    "encoded_nucleotides": observation.encoded_nucleotides,
                    "runtime_seconds": observation.runtime_seconds,
                }
                for observation in self._observations
            ],
        }


class LinearUtilityPolicyModel:
    """Ridge-regression utility model trained only from caller-supplied observations.

    This is a stronger parametric baseline than nearest-neighbor ranking while remaining
    deterministic, dependency-free, inspectable, and intentionally modest in scope.
    """

    def __init__(self, ridge: float = 1e-3) -> None:
        if ridge < 0:
            raise ValueError("ridge must be non-negative")
        self.ridge = ridge
        self._coefficients: list[float] = []
        self._policies: dict[tuple[int, int, int, bool], CodecPolicy] = {}
        self._max_nucleotides = 1.0
        self._max_runtime = 1.0

    def fit(self, observations: list[PolicyObservation]) -> LinearUtilityPolicyModel:
        if len(observations) < 2:
            raise ValueError("LinearUtilityPolicyModel requires at least two observations")
        self._max_nucleotides = float(
            max(max(1, observation.encoded_nucleotides) for observation in observations)
        )
        self._max_runtime = max(
            1e-9,
            max(observation.runtime_seconds for observation in observations),
        )
        rows: list[list[float]] = []
        targets: list[float] = []
        for observation in observations:
            observation.channel.validate()
            rows.append(_features(observation.channel, observation.policy))
            recovery = 1.0 if observation.recovered else 0.0
            overhead = observation.encoded_nucleotides / self._max_nucleotides
            runtime = observation.runtime_seconds / self._max_runtime
            targets.append(recovery - 0.20 * overhead - 0.05 * runtime)
            self._policies[_policy_key(observation.policy)] = observation.policy

        dimension = len(rows[0])
        xtx = [[0.0] * dimension for _ in range(dimension)]
        xty = [0.0] * dimension
        for row, target in zip(rows, targets, strict=True):
            for i in range(dimension):
                xty[i] += row[i] * target
                for j in range(dimension):
                    xtx[i][j] += row[i] * row[j]
        for index in range(dimension):
            xtx[index][index] += self.ridge
        self._coefficients = _solve_linear_system(xtx, xty)
        return self

    @property
    def coefficients(self) -> tuple[float, ...]:
        return tuple(self._coefficients)

    def to_dict(self) -> dict[str, object]:
        if not self._coefficients:
            raise ValueError("Model must be fit before serialization")
        return {
            "model": "linear-utility-ridge",
            "ridge": self.ridge,
            "coefficients": list(self._coefficients),
            "max_nucleotides": self._max_nucleotides,
            "max_runtime": self._max_runtime,
            "policies": [
                asdict(policy)
                for _, policy in sorted(self._policies.items(), key=lambda item: item[0])
            ],
        }

    @classmethod
    def from_dict(cls, value: dict[str, object]) -> LinearUtilityPolicyModel:
        if value.get("model") != "linear-utility-ridge":
            raise ValueError("Unsupported serialized policy model")
        ridge = value.get("ridge")
        coefficients = value.get("coefficients")
        policies = value.get("policies")
        max_nucleotides = value.get("max_nucleotides")
        max_runtime = value.get("max_runtime")
        if isinstance(ridge, bool) or not isinstance(ridge, (int, float)):
            raise ValueError("Serialized ridge must be numeric")
        if not isinstance(coefficients, list) or not coefficients:
            raise ValueError("Serialized coefficients must be a non-empty list")
        if not isinstance(policies, list) or not policies:
            raise ValueError("Serialized policies must be a non-empty list")
        if (
            isinstance(max_nucleotides, bool)
            or not isinstance(max_nucleotides, (int, float))
            or isinstance(max_runtime, bool)
            or not isinstance(max_runtime, (int, float))
        ):
            raise ValueError("Serialized normalization values must be numeric")

        model = cls(float(ridge))
        model._coefficients = []
        for coefficient in coefficients:
            if isinstance(coefficient, bool) or not isinstance(coefficient, (int, float)):
                raise ValueError("Serialized coefficient values must be numeric")
            model._coefficients.append(float(coefficient))
        for raw in policies:
            if not isinstance(raw, dict):
                raise ValueError("Serialized policy entries must be objects")
            rationale = raw.get("rationale", ())
            if not isinstance(rationale, (list, tuple)):
                raise ValueError("Serialized policy rationale must be a list or tuple")
            policy = CodecPolicy(
                chunk_size=int(raw["chunk_size"]),
                rs_nsym=int(raw["rs_nsym"]),
                parity_group_size=int(raw["parity_group_size"]),
                adaptive_masks=bool(raw["adaptive_masks"]),
                rationale=tuple(str(item) for item in rationale),
            )
            model._policies[_policy_key(policy)] = policy
        model._max_nucleotides = float(max_nucleotides)
        model._max_runtime = float(max_runtime)
        return model

    def predict_utility(self, channel: ChannelProfile, policy: CodecPolicy) -> float:
        if not self._coefficients:
            raise ValueError("Model must be fit before prediction")
        channel.validate()
        features = _features(channel, policy)
        return sum(
            coefficient * feature
            for coefficient, feature in zip(
                self._coefficients,
                features,
                strict=True,
            )
        )

    def recommend(
        self,
        channel: ChannelProfile,
        candidates: list[CodecPolicy] | None = None,
    ) -> LinearPolicyRecommendation:
        if not self._coefficients:
            raise ValueError("Model must be fit before recommendation")
        resolved = candidates or list(self._policies.values())
        if not resolved:
            raise ValueError("At least one candidate policy is required")
        scored = sorted(
            (
                (self.predict_utility(channel, policy), policy)
                for policy in resolved
            ),
            key=lambda item: (
                item[0],
                item[1].rs_nsym,
                -item[1].chunk_size,
            ),
            reverse=True,
        )
        best_score, best_policy = scored[0]
        second_score = scored[1][0] if len(scored) > 1 else best_score
        margin = max(0.0, best_score - second_score)
        bounded_margin = 1.0 - math.exp(-margin)
        return LinearPolicyRecommendation(
            policy=best_policy,
            predicted_utility=round(best_score, 6),
            margin=round(bounded_margin, 6),
            candidate_count=len(resolved),
            rationale=(
                "deterministic ridge regression trained only from supplied observations",
                "target utility rewards verified recovery and penalizes overhead/runtime",
                "margin is a relative separation score, not a calibrated probability",
            ),
        )

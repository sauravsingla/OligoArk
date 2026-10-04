"""Deterministic empirical policy-learning baseline.

The model is intentionally dependency-free and transparent. It learns from observed
(simulated or measured) policy outcomes using inverse-distance weighting over channel profiles.
It is a research baseline, not a claim of optimality.
"""

from __future__ import annotations

from dataclasses import dataclass

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
    """Instance-based policy ranker learned from prior observations."""

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

    def recommend(self, channel: ChannelProfile) -> LearnedPolicyRecommendation:
        if not self._observations:
            raise ValueError("Model must be fit before recommendation")

        max_nucleotides = max(obs.encoded_nucleotides for obs in self._observations)
        max_runtime = max(obs.runtime_seconds for obs in self._observations) or 1.0
        scores: dict[tuple[int, int, int, bool], float] = {}
        counts: dict[tuple[int, int, int, bool], int] = {}
        policies: dict[tuple[int, int, int, bool], CodecPolicy] = {}

        for obs in self._observations:
            key = (
                obs.policy.chunk_size,
                obs.policy.rs_nsym,
                obs.policy.parity_group_size,
                obs.policy.adaptive_masks,
            )
            distance = self._distance(channel, obs.channel)
            weight = 1.0 / (0.001 + distance)
            recovery_score = 1.0 if obs.recovered else -1.0
            overhead_penalty = obs.encoded_nucleotides / max(1, max_nucleotides)
            runtime_penalty = obs.runtime_seconds / max_runtime
            value = recovery_score - 0.15 * overhead_penalty - 0.05 * runtime_penalty
            scores[key] = scores.get(key, 0.0) + weight * value
            counts[key] = counts.get(key, 0) + 1
            policies[key] = obs.policy

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

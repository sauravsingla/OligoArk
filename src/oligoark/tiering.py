"""Explainable heterogeneous storage-tier and lifecycle recommendation engine."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from typing import cast


@dataclass(frozen=True)
class WorkloadProfile:
    retention_years: float
    accesses_per_year: float
    mutability: float
    retrieval_urgency: float
    durability_priority: float
    energy_priority: float
    redundancy_priority: float = 0.5
    cost_priority: float = 0.5
    data_size_gb: float = 1.0
    expected_access_probability: float | None = None

    def validate(self) -> None:
        if self.retention_years <= 0 or self.accesses_per_year < 0:
            raise ValueError("retention_years must be > 0 and accesses_per_year >= 0")
        if self.data_size_gb <= 0:
            raise ValueError("data_size_gb must be positive")
        if (
            self.expected_access_probability is not None
            and not 0 <= self.expected_access_probability <= 1
        ):
            raise ValueError("expected_access_probability must be between 0 and 1")
        for name in (
            "mutability",
            "retrieval_urgency",
            "durability_priority",
            "energy_priority",
            "redundancy_priority",
            "cost_priority",
        ):
            value = getattr(self, name)
            if not 0 <= value <= 1:
                raise ValueError(f"{name} must be between 0 and 1")


@dataclass(frozen=True)
class EconomicAssumptions:
    """User-configurable normalized economic indices; lower means lower assumed cost."""

    storage_cost_index: dict[str, float] = field(
        default_factory=lambda: {
            "ssd": 0.5,
            "object_archive": 0.5,
            "tape": 0.5,
            "dna_future": 0.5,
        }
    )
    retrieval_cost_index: dict[str, float] = field(
        default_factory=lambda: {
            "ssd": 0.5,
            "object_archive": 0.5,
            "tape": 0.5,
            "dna_future": 0.5,
        }
    )

    def validate(self) -> None:
        expected = set(_TRAITS)
        for name, values in (
            ("storage_cost_index", self.storage_cost_index),
            ("retrieval_cost_index", self.retrieval_cost_index),
        ):
            if set(values) != expected:
                raise ValueError(f"{name} must define exactly these tiers: {sorted(expected)}")
            if any(not 0 <= value <= 1 for value in values.values()):
                raise ValueError(f"{name} values must be between 0 and 1")

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> EconomicAssumptions:
        storage = values.get("storage_cost_index")
        retrieval = values.get("retrieval_cost_index")
        if not isinstance(storage, dict) or not isinstance(retrieval, dict):
            raise ValueError(
                "Economic assumptions require storage_cost_index and retrieval_cost_index objects"
            )

        def numeric_mapping(raw: dict[object, object]) -> dict[str, float]:
            converted: dict[str, float] = {}
            for key, value in raw.items():
                if isinstance(value, bool) or not isinstance(value, (int, float)):
                    raise ValueError("Economic assumption values must be numeric")
                converted[str(key)] = float(value)
            return converted

        assumptions = cls(
            storage_cost_index=numeric_mapping(cast(dict[object, object], storage)),
            retrieval_cost_index=numeric_mapping(cast(dict[object, object], retrieval)),
        )
        assumptions.validate()
        return assumptions


@dataclass(frozen=True)
class TierLifecycleAssumption:
    """Explicit user-supplied lifecycle inputs for one storage tier."""

    storage_cost_per_gb_year: float
    retrieval_cost_per_gb: float
    idle_energy_kwh_per_tb_year: float
    retrieval_energy_kwh_per_gb: float
    retrieval_latency_hours: float

    def validate(self) -> None:
        for name, value in vars(self).items():
            if value < 0:
                raise ValueError(f"{name} must be non-negative")


@dataclass(frozen=True)
class LifecycleAssumptions:
    """Per-tier lifecycle inputs. OligoArk intentionally provides no price/energy defaults."""

    tiers: dict[str, TierLifecycleAssumption]

    def validate(self) -> None:
        expected = set(_TRAITS)
        if set(self.tiers) != expected:
            raise ValueError(f"lifecycle tiers must be exactly {sorted(expected)}")
        for assumption in self.tiers.values():
            assumption.validate()

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> LifecycleAssumptions:
        tiers: dict[str, TierLifecycleAssumption] = {}
        for tier_name, raw in values.items():
            if not isinstance(raw, dict):
                raise ValueError(f"lifecycle tier {tier_name!r} must be an object")
            numeric: dict[str, float] = {}
            for field_name in (
                "storage_cost_per_gb_year",
                "retrieval_cost_per_gb",
                "idle_energy_kwh_per_tb_year",
                "retrieval_energy_kwh_per_gb",
                "retrieval_latency_hours",
            ):
                value = raw.get(field_name)
                if isinstance(value, bool) or not isinstance(value, (int, float)):
                    raise ValueError(f"lifecycle {tier_name}.{field_name} must be numeric")
                numeric[field_name] = float(value)
            tiers[str(tier_name)] = TierLifecycleAssumption(
                storage_cost_per_gb_year=numeric["storage_cost_per_gb_year"],
                retrieval_cost_per_gb=numeric["retrieval_cost_per_gb"],
                idle_energy_kwh_per_tb_year=numeric["idle_energy_kwh_per_tb_year"],
                retrieval_energy_kwh_per_gb=numeric["retrieval_energy_kwh_per_gb"],
                retrieval_latency_hours=numeric["retrieval_latency_hours"],
            )
        assumptions = cls(tiers)
        assumptions.validate()
        return assumptions


@dataclass(frozen=True)
class LifecycleEstimate:
    total_storage_cost: float
    total_retrieval_cost: float
    total_cost: float
    idle_energy_kwh: float
    retrieval_energy_kwh: float
    total_energy_kwh: float
    expected_retrievals: float
    expected_retrieval_latency_hours: float

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class TierScoreBreakdown:
    retention: float
    access_frequency: float
    access_probability: float
    mutability: float
    retrieval_urgency: float
    durability: float
    redundancy: float
    energy: float
    normalized_storage_cost: float
    normalized_retrieval_cost: float
    lifecycle_storage_cost: float = 0.0
    lifecycle_retrieval_cost: float = 0.0
    lifecycle_energy: float = 0.0
    lifecycle_latency: float = 0.0
    dna_future_penalty: float = 0.0

    @property
    def final_score(self) -> float:
        positive = (
            self.retention
            + self.access_frequency
            + self.access_probability
            + self.mutability
            + self.retrieval_urgency
            + self.durability
            + self.redundancy
            + self.energy
            + self.normalized_storage_cost
            + self.normalized_retrieval_cost
            + self.lifecycle_storage_cost
            + self.lifecycle_retrieval_cost
            + self.lifecycle_energy
            + self.lifecycle_latency
        )
        return round(max(0.0, positive - self.dna_future_penalty), 6)


@dataclass(frozen=True)
class TierRecommendation:
    recommended_tier: str
    scores: dict[str, float]
    rationale: tuple[str, ...]
    assumptions: str
    lifecycle_estimates: dict[str, LifecycleEstimate] | None = None
    score_breakdowns: dict[str, TierScoreBreakdown] | None = None

    def to_dict(self) -> dict[str, object]:
        return cast(dict[str, object], asdict(self))


_TRAITS = {
    "ssd": {
        "latency": 1.00,
        "mutable": 1.00,
        "durability": 0.45,
        "idle_energy": 0.30,
        "long_term": 0.20,
        "redundancy": 0.70,
    },
    "object_archive": {
        "latency": 0.55,
        "mutable": 0.65,
        "durability": 0.80,
        "idle_energy": 0.70,
        "long_term": 0.65,
        "redundancy": 0.85,
    },
    "tape": {
        "latency": 0.25,
        "mutable": 0.30,
        "durability": 0.78,
        "idle_energy": 0.90,
        "long_term": 0.80,
        "redundancy": 0.70,
    },
    "dna_future": {
        "latency": 0.05,
        "mutable": 0.05,
        "durability": 0.95,
        "idle_energy": 1.00,
        "long_term": 1.00,
        "redundancy": 0.90,
    },
}


def _estimate_lifecycle(
    profile: WorkloadProfile,
    assumptions: LifecycleAssumptions,
) -> dict[str, LifecycleEstimate]:
    estimates: dict[str, LifecycleEstimate] = {}
    access_probability = (
        profile.expected_access_probability
        if profile.expected_access_probability is not None
        else min(1.0, profile.accesses_per_year)
    )
    expected_retrievals = profile.accesses_per_year * profile.retention_years * access_probability
    for tier, value in assumptions.tiers.items():
        storage_cost = (
            value.storage_cost_per_gb_year * profile.data_size_gb * profile.retention_years
        )
        retrieval_cost = (
            value.retrieval_cost_per_gb * profile.data_size_gb * expected_retrievals
        )
        idle_energy = (
            value.idle_energy_kwh_per_tb_year
            * (profile.data_size_gb / 1000.0)
            * profile.retention_years
        )
        retrieval_energy = (
            value.retrieval_energy_kwh_per_gb * profile.data_size_gb * expected_retrievals
        )
        estimates[tier] = LifecycleEstimate(
            total_storage_cost=round(storage_cost, 8),
            total_retrieval_cost=round(retrieval_cost, 8),
            total_cost=round(storage_cost + retrieval_cost, 8),
            idle_energy_kwh=round(idle_energy, 8),
            retrieval_energy_kwh=round(retrieval_energy, 8),
            total_energy_kwh=round(idle_energy + retrieval_energy, 8),
            expected_retrievals=round(expected_retrievals, 8),
            expected_retrieval_latency_hours=round(
                expected_retrievals * value.retrieval_latency_hours,
                8,
            ),
        )
    return estimates


def _efficiency(value: float, maximum: float) -> float:
    if maximum <= 0:
        return 1.0
    return max(0.0, 1.0 - value / maximum)


def recommend_storage_tier(
    profile: WorkloadProfile,
    economics: EconomicAssumptions | None = None,
    lifecycle: LifecycleAssumptions | None = None,
) -> TierRecommendation:
    """Return fully decomposable tier scores plus optional lifecycle estimates."""
    profile.validate()
    economics = economics or EconomicAssumptions()
    economics.validate()
    if lifecycle is not None:
        lifecycle.validate()

    retention_need = min(1.0, profile.retention_years / 100.0)
    access_frequency_need = min(1.0, profile.accesses_per_year / 365.0)
    access_probability_need = profile.expected_access_probability or 0.0
    estimates = _estimate_lifecycle(profile, lifecycle) if lifecycle else None

    max_storage_cost = max(
        (value.total_storage_cost for value in estimates.values()), default=1.0
    ) if estimates else 1.0
    max_retrieval_cost = max(
        (value.total_retrieval_cost for value in estimates.values()), default=1.0
    ) if estimates else 1.0
    max_energy = max(
        (value.total_energy_kwh for value in estimates.values()), default=1.0
    ) if estimates else 1.0
    max_latency = max(
        (value.expected_retrieval_latency_hours for value in estimates.values()), default=1.0
    ) if estimates else 1.0

    breakdowns: dict[str, TierScoreBreakdown] = {}
    scores: dict[str, float] = {}
    for name, traits in _TRAITS.items():
        lifecycle_storage = lifecycle_retrieval = lifecycle_energy = lifecycle_latency = 0.0
        if estimates:
            estimate = estimates[name]
            lifecycle_storage = (
                0.04
                * _efficiency(estimate.total_storage_cost, max_storage_cost)
                * profile.cost_priority
            )
            lifecycle_retrieval = (
                0.03
                * _efficiency(estimate.total_retrieval_cost, max_retrieval_cost)
                * profile.cost_priority
            )
            lifecycle_energy = (
                0.04
                * _efficiency(estimate.total_energy_kwh, max_energy)
                * profile.energy_priority
            )
            lifecycle_latency = (
                0.04
                * _efficiency(estimate.expected_retrieval_latency_hours, max_latency)
                * profile.retrieval_urgency
            )

        penalty = (
            0.25 * profile.retrieval_urgency + 0.20 * profile.mutability
            if name == "dna_future"
            else 0.0
        )
        breakdown = TierScoreBreakdown(
            retention=0.14 * traits["long_term"] * retention_need,
            access_frequency=0.05 * traits["latency"] * access_frequency_need,
            access_probability=0.03 * traits["latency"] * access_probability_need,
            mutability=0.10 * traits["mutable"] * profile.mutability,
            retrieval_urgency=0.12 * traits["latency"] * profile.retrieval_urgency,
            durability=0.15 * traits["durability"] * profile.durability_priority,
            redundancy=0.07 * traits["redundancy"] * profile.redundancy_priority,
            energy=0.08 * traits["idle_energy"] * profile.energy_priority,
            normalized_storage_cost=(
                0.06
                * (1.0 - economics.storage_cost_index[name])
                * profile.cost_priority
            ),
            normalized_retrieval_cost=(
                0.04
                * (1.0 - economics.retrieval_cost_index[name])
                * profile.cost_priority
            ),
            lifecycle_storage_cost=lifecycle_storage,
            lifecycle_retrieval_cost=lifecycle_retrieval,
            lifecycle_energy=lifecycle_energy,
            lifecycle_latency=lifecycle_latency,
            dna_future_penalty=penalty,
        )
        breakdowns[name] = breakdown
        scores[name] = breakdown.final_score

    winner = max(scores, key=lambda tier_name: scores[tier_name])
    rationale = (
        f"retention horizon normalized to {retention_need:.3f}",
        f"access frequency normalized to {access_frequency_need:.3f}",
        f"access probability normalized to {access_probability_need:.3f}",
        f"data size={profile.data_size_gb:.6g} GB",
        f"highest decomposed score: {winner}={scores[winner]:.6f}",
        (
            "lifecycle score terms use only caller-supplied physical inputs"
            if estimates
            else "lifecycle physical terms omitted because none were supplied"
        ),
        "dna_future remains an experimental scenario tier, not a present-day superiority claim",
    )
    assumptions_text = (
        "Technical traits are normalized research heuristics. Normalized economic indices are "
        "neutral unless supplied. Physical lifecycle cost, energy and latency terms are included "
        "only when explicitly supplied by the caller."
    )
    return TierRecommendation(
        recommended_tier=winner,
        scores=scores,
        rationale=rationale,
        assumptions=assumptions_text,
        lifecycle_estimates=estimates,
        score_breakdowns=breakdowns,
    )

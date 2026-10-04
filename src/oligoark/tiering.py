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
                    raise ValueError(
                        f"lifecycle {tier_name}.{field_name} must be numeric"
                    )
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
class TierRecommendation:
    recommended_tier: str
    scores: dict[str, float]
    rationale: tuple[str, ...]
    assumptions: str
    lifecycle_estimates: dict[str, LifecycleEstimate] | None = None

    def to_dict(self) -> dict[str, object]:
        raw = asdict(self)
        return cast(dict[str, object], raw)


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
    expected_retrievals = (
        profile.accesses_per_year * profile.retention_years * access_probability
    )
    for tier, value in assumptions.tiers.items():
        storage_cost = (
            value.storage_cost_per_gb_year
            * profile.data_size_gb
            * profile.retention_years
        )
        retrieval_cost = (
            value.retrieval_cost_per_gb
            * profile.data_size_gb
            * expected_retrievals
        )
        idle_energy = (
            value.idle_energy_kwh_per_tb_year
            * (profile.data_size_gb / 1000.0)
            * profile.retention_years
        )
        retrieval_energy = (
            value.retrieval_energy_kwh_per_gb
            * profile.data_size_gb
            * expected_retrievals
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


def recommend_storage_tier(
    profile: WorkloadProfile,
    economics: EconomicAssumptions | None = None,
    lifecycle: LifecycleAssumptions | None = None,
) -> TierRecommendation:
    """Return explainable tier scores plus optional explicit lifecycle estimates."""
    profile.validate()
    economics = economics or EconomicAssumptions()
    economics.validate()
    if lifecycle is not None:
        lifecycle.validate()

    retention_need = min(1.0, profile.retention_years / 100.0)
    access_need = min(1.0, profile.accesses_per_year / 365.0)
    lifecycle_estimates = _estimate_lifecycle(profile, lifecycle) if lifecycle else None

    lifecycle_cost_max = 1.0
    lifecycle_energy_max = 1.0
    lifecycle_latency_max = 1.0
    if lifecycle_estimates:
        lifecycle_cost_max = max(value.total_cost for value in lifecycle_estimates.values()) or 1.0
        lifecycle_energy_max = (
            max(value.total_energy_kwh for value in lifecycle_estimates.values()) or 1.0
        )
        lifecycle_latency_max = (
            max(
                value.expected_retrieval_latency_hours
                for value in lifecycle_estimates.values()
            )
            or 1.0
        )

    scores: dict[str, float] = {}
    for name, traits in _TRAITS.items():
        storage_affordability = 1.0 - economics.storage_cost_index[name]
        retrieval_affordability = 1.0 - economics.retrieval_cost_index[name]
        economic_score = 0.7 * storage_affordability + 0.3 * retrieval_affordability
        lifecycle_adjustment = 0.0
        if lifecycle_estimates:
            estimate = lifecycle_estimates[name]
            cost_efficiency = 1.0 - estimate.total_cost / lifecycle_cost_max
            energy_efficiency = 1.0 - estimate.total_energy_kwh / lifecycle_energy_max
            latency_efficiency = (
                1.0 - estimate.expected_retrieval_latency_hours / lifecycle_latency_max
            )
            lifecycle_adjustment = (
                0.10 * cost_efficiency * profile.cost_priority
                + 0.08 * energy_efficiency * profile.energy_priority
                + 0.06 * latency_efficiency * profile.retrieval_urgency
            )

        score = (
            0.14 * traits["latency"] * profile.retrieval_urgency
            + 0.11 * traits["mutable"] * profile.mutability
            + 0.16 * traits["durability"] * profile.durability_priority
            + 0.10 * traits["idle_energy"] * profile.energy_priority
            + 0.15 * traits["long_term"] * retention_need
            + 0.07 * traits["latency"] * access_need
            + 0.07 * traits["redundancy"] * profile.redundancy_priority
            + 0.10 * economic_score * profile.cost_priority
            + lifecycle_adjustment
        )
        if name == "dna_future":
            score -= 0.25 * profile.retrieval_urgency + 0.20 * profile.mutability
        scores[name] = round(max(0.0, score), 6)

    winner = max(scores, key=lambda tier_name: scores[tier_name])
    rationale = (
        f"retention horizon normalized to {retention_need:.3f}",
        f"access frequency normalized to {access_need:.3f}",
        f"data size={profile.data_size_gb:.6g} GB",
        (
            "expected access probability="
            f"{profile.expected_access_probability:.3f}"
            if profile.expected_access_probability is not None
            else "expected access probability not supplied; frequency remains explicit"
        ),
        f"highest transparent score: {winner}={scores[winner]:.6f}",
        (
            "explicit lifecycle estimates included from caller-supplied inputs"
            if lifecycle_estimates
            else "no explicit lifecycle price/energy inputs supplied"
        ),
        "dna_future is an experimental planning tier, not a claim of present-day superiority",
    )
    assumptions_text = (
        "Technical traits are normalized research heuristics. Economic indices are neutral "
        "unless explicitly supplied. Lifecycle price/energy/latency values are included only "
        "when supplied by the caller; OligoArk does not invent vendor or physical-DNA values."
    )
    return TierRecommendation(
        winner,
        scores,
        rationale,
        assumptions_text,
        lifecycle_estimates,
    )

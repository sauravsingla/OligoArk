"""Explainable heterogeneous storage-tier recommendation engine."""

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

    def validate(self) -> None:
        if self.retention_years <= 0 or self.accesses_per_year < 0:
            raise ValueError("retention_years must be > 0 and accesses_per_year >= 0")
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
    """User-configurable normalized economic inputs; lower index means lower cost.

    Defaults are deliberately neutral and equal for every tier so OligoArk does not fabricate
    present-day vendor pricing or future DNA economics.
    """

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
class TierRecommendation:
    recommended_tier: str
    scores: dict[str, float]
    rationale: tuple[str, ...]
    assumptions: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


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


def recommend_storage_tier(
    profile: WorkloadProfile,
    economics: EconomicAssumptions | None = None,
) -> TierRecommendation:
    profile.validate()
    economics = economics or EconomicAssumptions()
    economics.validate()
    retention_need = min(1.0, profile.retention_years / 100.0)
    access_need = min(1.0, profile.accesses_per_year / 365.0)
    scores: dict[str, float] = {}

    for name, traits in _TRAITS.items():
        storage_affordability = 1.0 - economics.storage_cost_index[name]
        retrieval_affordability = 1.0 - economics.retrieval_cost_index[name]
        economic_score = 0.7 * storage_affordability + 0.3 * retrieval_affordability
        score = (
            0.16 * traits["latency"] * profile.retrieval_urgency
            + 0.12 * traits["mutable"] * profile.mutability
            + 0.17 * traits["durability"] * profile.durability_priority
            + 0.12 * traits["idle_energy"] * profile.energy_priority
            + 0.17 * traits["long_term"] * retention_need
            + 0.08 * traits["latency"] * access_need
            + 0.08 * traits["redundancy"] * profile.redundancy_priority
            + 0.10 * economic_score * profile.cost_priority
        )
        if name == "dna_future":
            score -= 0.25 * profile.retrieval_urgency + 0.20 * profile.mutability
        scores[name] = round(max(0.0, score), 4)

    winner = max(scores, key=lambda tier_name: scores[tier_name])
    rationale = (
        f"retention horizon normalized to {retention_need:.2f}",
        f"access intensity normalized to {access_need:.2f}",
        f"highest transparent heuristic score: {winner}={scores[winner]:.4f}",
        "economic inputs are normalized user assumptions; defaults are neutral across tiers",
        "dna_future is an experimental planning tier, not a claim of present-day superiority",
    )
    assumptions = (
        "Technical traits are documented normalized research heuristics. Economic indices are "
        "neutral unless explicitly supplied by the user; use real operational data for decisions."
    )
    return TierRecommendation(winner, scores, rationale, assumptions)

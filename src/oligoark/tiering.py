"""Explainable heterogeneous storage-tier recommendation engine."""

from __future__ import annotations
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class WorkloadProfile:
    retention_years: float
    accesses_per_year: float
    mutability: float
    retrieval_urgency: float
    durability_priority: float
    energy_priority: float

    def validate(self) -> None:
        if self.retention_years <= 0 or self.accesses_per_year < 0:
            raise ValueError("retention_years must be > 0 and accesses_per_year >= 0")
        for field in ("mutability", "retrieval_urgency", "durability_priority", "energy_priority"):
            if not 0 <= getattr(self, field) <= 1:
                raise ValueError(f"{field} must be between 0 and 1")


@dataclass(frozen=True)
class TierRecommendation:
    recommended_tier: str
    scores: dict[str, float]
    rationale: tuple[str, ...]
    assumptions: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


_TRAITS = {
    "ssd": {"latency": 1.00, "mutable": 1.00, "durability": 0.45, "idle_energy": 0.30, "long_term": 0.20},
    "object_archive": {"latency": 0.55, "mutable": 0.65, "durability": 0.80, "idle_energy": 0.70, "long_term": 0.65},
    "tape": {"latency": 0.25, "mutable": 0.30, "durability": 0.78, "idle_energy": 0.90, "long_term": 0.80},
    "dna_future": {"latency": 0.05, "mutable": 0.05, "durability": 0.95, "idle_energy": 1.00, "long_term": 1.00},
}


def recommend_storage_tier(profile: WorkloadProfile) -> TierRecommendation:
    profile.validate()
    retention_need = min(1.0, profile.retention_years / 100.0)
    access_need = min(1.0, profile.accesses_per_year / 365.0)
    scores: dict[str, float] = {}
    for name, t in _TRAITS.items():
        score = (0.20 * t["latency"] * profile.retrieval_urgency
                 + 0.15 * t["mutable"] * profile.mutability
                 + 0.20 * t["durability"] * profile.durability_priority
                 + 0.15 * t["idle_energy"] * profile.energy_priority
                 + 0.20 * t["long_term"] * retention_need
                 + 0.10 * t["latency"] * access_need)
        if name == "dna_future":
            score -= 0.25 * profile.retrieval_urgency + 0.20 * profile.mutability
        scores[name] = round(max(0.0, score), 4)
    winner = max(scores, key=scores.get)
    rationale = (
        f"retention horizon normalized to {retention_need:.2f}",
        f"access intensity normalized to {access_need:.2f}",
        f"highest transparent heuristic score: {winner}={scores[winner]:.4f}",
        "dna_future is an experimental planning tier, not a claim of present-day superiority",
    )
    return TierRecommendation(winner, scores, rationale,
        "Scores use documented normalized heuristics only; configure real economic and operational data before production decisions.")

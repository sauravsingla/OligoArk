"""Adaptive, explainable codec policy selection."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class ChannelProfile:
    substitution_rate: float = 0.0
    insertion_rate: float = 0.0
    deletion_rate: float = 0.0
    dropout_rate: float = 0.0

    def validate(self) -> None:
        for name, value in vars(self).items():
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1")

    @property
    def total_symbol_error_rate(self) -> float:
        return self.substitution_rate + self.insertion_rate + self.deletion_rate


@dataclass(frozen=True)
class PolicyObjective:
    durability_priority: float = 0.7
    storage_overhead_priority: float = 0.2
    retrieval_speed_priority: float = 0.1

    def validate(self) -> None:
        for name, value in vars(self).items():
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1")
        if sum(vars(self).values()) <= 0:
            raise ValueError("At least one policy objective priority must be positive")


@dataclass(frozen=True)
class CodecPolicy:
    chunk_size: int
    rs_nsym: int
    parity_group_size: int
    adaptive_masks: bool
    rationale: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def recommend_codec_policy(
    channel: ChannelProfile,
    objective: PolicyObjective | None = None,
) -> CodecPolicy:
    """Choose a deterministic baseline policy from channel and workload-style priorities."""
    channel.validate()
    objective = objective or PolicyObjective()
    objective.validate()
    symbol_error = channel.total_symbol_error_rate

    if symbol_error >= 0.02:
        chunk_size, rs_nsym = 48, 24
    elif symbol_error >= 0.005:
        chunk_size, rs_nsym = 64, 16
    else:
        chunk_size, rs_nsym = 96, 8

    if objective.durability_priority >= 0.85 and rs_nsym < 24:
        rs_nsym += 4
        chunk_size = min(chunk_size, 80)
    if objective.storage_overhead_priority >= 0.70 and symbol_error < 0.005:
        rs_nsym = max(6, rs_nsym - 2)

    if channel.dropout_rate >= 0.10:
        parity_group = 3
    elif channel.dropout_rate >= 0.03:
        parity_group = 5
    else:
        parity_group = 8
    if objective.durability_priority >= 0.85:
        parity_group = max(3, parity_group - 1)

    rationale = (
        f"symbol error rate={symbol_error:.4f} -> chunk_size={chunk_size}, rs_nsym={rs_nsym}",
        f"dropout rate={channel.dropout_rate:.4f} -> XOR parity group={parity_group}",
        "policy objective adjusts durability versus storage overhead deterministically",
        "adaptive payload masks reduce a simple GC/homopolymer quality penalty",
    )
    return CodecPolicy(chunk_size, rs_nsym, parity_group, True, rationale)

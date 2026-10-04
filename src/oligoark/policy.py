"""Adaptive, explainable codec policy selection."""

from __future__ import annotations
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class ChannelProfile:
    substitution_rate: float = 0.0
    insertion_rate: float = 0.0
    deletion_rate: float = 0.0
    dropout_rate: float = 0.0

    @property
    def total_symbol_error_rate(self) -> float:
        return self.substitution_rate + self.insertion_rate + self.deletion_rate


@dataclass(frozen=True)
class CodecPolicy:
    chunk_size: int
    rs_nsym: int
    parity_group_size: int
    adaptive_masks: bool
    rationale: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def recommend_codec_policy(channel: ChannelProfile) -> CodecPolicy:
    symbol_error = channel.total_symbol_error_rate
    if symbol_error >= 0.02:
        chunk_size, rs_nsym = 48, 24
    elif symbol_error >= 0.005:
        chunk_size, rs_nsym = 64, 16
    else:
        chunk_size, rs_nsym = 96, 8
    if channel.dropout_rate >= 0.10:
        parity_group = 3
    elif channel.dropout_rate >= 0.03:
        parity_group = 5
    else:
        parity_group = 8
    rationale = (
        f"symbol error rate={symbol_error:.4f} -> chunk_size={chunk_size}, rs_nsym={rs_nsym}",
        f"dropout rate={channel.dropout_rate:.4f} -> XOR parity group={parity_group}",
        "adaptive payload masks enabled to reduce simple GC/homopolymer quality penalty",
    )
    return CodecPolicy(chunk_size, rs_nsym, parity_group, True, rationale)

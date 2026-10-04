"""Deterministic DNA channel simulator for software experiments."""

from __future__ import annotations
import random
from dataclasses import dataclass

DNA = "ACGT"


@dataclass(frozen=True)
class SimulationConfig:
    substitution_rate: float = 0.0
    insertion_rate: float = 0.0
    deletion_rate: float = 0.0
    dropout_rate: float = 0.0
    duplicate_rate: float = 0.0
    seed: int = 7

    def validate(self) -> None:
        for name, value in vars(self).items():
            if name != "seed" and not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1")


def _mutate(sequence: str, cfg: SimulationConfig, rng: random.Random) -> str:
    out: list[str] = []
    for base in sequence:
        if rng.random() < cfg.deletion_rate:
            continue
        if rng.random() < cfg.insertion_rate:
            out.append(rng.choice(DNA))
        if rng.random() < cfg.substitution_rate:
            out.append(rng.choice(DNA.replace(base, "")))
        else:
            out.append(base)
    if rng.random() < cfg.insertion_rate:
        out.append(rng.choice(DNA))
    return "".join(out)


def simulate_channel(strands: list[str], cfg: SimulationConfig) -> list[str]:
    cfg.validate()
    rng = random.Random(cfg.seed)
    reads: list[str] = []
    for strand in strands:
        if rng.random() < cfg.dropout_rate:
            continue
        reads.append(_mutate(strand, cfg, rng))
        if rng.random() < cfg.duplicate_rate:
            reads.append(_mutate(strand, cfg, rng))
    rng.shuffle(reads)
    return reads

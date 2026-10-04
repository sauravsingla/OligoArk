"""Deterministic fountain-style XOR redundancy baseline.

This is not the published DNA Fountain algorithm. It is an independently implemented
LT-style research baseline used to compare fixed parity groups with seeded, overlapping
XOR symbols inside OligoArk's real archive/recovery pipeline.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

from .ecc import xor_bytes


@dataclass(frozen=True)
class FountainSymbol:
    seed: int
    indexes: tuple[int, ...]
    payload: bytes


def indexes_for_seed(total: int, seed: int, max_degree: int = 4) -> tuple[int, ...]:
    """Return the deterministic source-chunk set represented by one fountain seed."""
    if total <= 0:
        raise ValueError("total must be positive")
    if max_degree < 1:
        raise ValueError("max_degree must be positive")
    rng = random.Random(seed)
    degree = min(total, 1 + rng.randrange(max(1, min(max_degree, total))))
    return tuple(sorted(rng.sample(range(total), degree)))


def make_symbols(
    chunks: list[bytes],
    count: int,
    width: int,
    seed: int = 1,
    max_degree: int = 4,
) -> list[FountainSymbol]:
    if not chunks:
        raise ValueError("chunks must not be empty")
    if count < 0:
        raise ValueError("count must be non-negative")
    if width <= 0:
        raise ValueError("width must be positive")
    symbols: list[FountainSymbol] = []
    for offset in range(count):
        symbol_seed = seed + offset
        indexes = indexes_for_seed(len(chunks), symbol_seed, max_degree=max_degree)
        symbols.append(
            FountainSymbol(
                symbol_seed,
                indexes,
                xor_bytes([chunks[index] for index in indexes], width),
            )
        )
    return symbols


def peel_decode(
    known: dict[int, bytes],
    symbols: list[FountainSymbol],
    total: int,
    width: int,
) -> dict[int, bytes]:
    """Iteratively solve degree-one XOR equations against already-known chunks."""
    if total <= 0 or width <= 0:
        raise ValueError("total and width must be positive")
    recovered = dict(known)
    pending = [(set(symbol.indexes), symbol.payload) for symbol in symbols]
    changed = True
    while changed:
        changed = False
        next_pending: list[tuple[set[int], bytes]] = []
        for indexes, payload in pending:
            unknown = set(indexes)
            residual_parts = [payload]
            for index in list(unknown):
                if index in recovered:
                    residual_parts.append(recovered[index])
                    unknown.remove(index)
            residual = xor_bytes(residual_parts, width)
            if len(unknown) == 1:
                index = next(iter(unknown))
                if index not in recovered:
                    recovered[index] = residual
                    changed = True
            elif unknown:
                next_pending.append((unknown, residual))
        pending = next_pending
    return {index: recovered[index] for index in range(total) if index in recovered}

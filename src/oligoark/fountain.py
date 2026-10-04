"""Small deterministic fountain-style XOR redundancy baseline.

This is not an implementation of DNA Fountain. It is an independently implemented LT-style
research baseline used to compare fixed parity groups with seeded, overlapping XOR symbols.
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


def _choose_indexes(total: int, seed: int, max_degree: int = 4) -> tuple[int, ...]:
    if total <= 0:
        raise ValueError("total must be positive")
    rng = random.Random(seed)
    degree = min(total, 1 + rng.randrange(max(1, min(max_degree, total))))
    return tuple(sorted(rng.sample(range(total), degree)))


def make_symbols(
    chunks: list[bytes],
    count: int,
    width: int,
    seed: int = 1,
) -> list[FountainSymbol]:
    if count < 0:
        raise ValueError("count must be non-negative")
    symbols: list[FountainSymbol] = []
    for offset in range(count):
        symbol_seed = seed + offset
        indexes = _choose_indexes(len(chunks), symbol_seed)
        symbols.append(
            FountainSymbol(symbol_seed, indexes, xor_bytes([chunks[i] for i in indexes], width))
        )
    return symbols


def peel_decode(
    known: dict[int, bytes],
    symbols: list[FountainSymbol],
    total: int,
    width: int,
) -> dict[int, bytes]:
    """Iteratively solve degree-one XOR equations against already-known chunks."""
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
    return {i: recovered[i] for i in range(total) if i in recovered}

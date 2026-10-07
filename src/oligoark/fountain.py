"""Deterministic fountain-style XOR redundancy baseline.

This is not the published DNA Fountain algorithm. It is an independently implemented
LT-style research baseline used to compare fixed parity groups with seeded, overlapping
XOR symbols inside OligoArk's real archive/recovery pipeline.
"""

from __future__ import annotations

import random
from collections import deque
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
    unknown_sets: list[set[int]] = []
    residuals: list[bytearray] = []
    incident: list[list[int]] = [[] for _ in range(total)]
    ready: deque[int] = deque()

    for symbol in symbols:
        unknown = {index for index in symbol.indexes if index not in recovered}
        residual = bytearray(symbol.payload)
        for index in symbol.indexes:
            if index in recovered:
                payload = recovered[index]
                for offset, value in enumerate(payload):
                    residual[offset] ^= value
        equation_id = len(unknown_sets)
        unknown_sets.append(unknown)
        residuals.append(residual)
        for index in unknown:
            incident[index].append(equation_id)
        if len(unknown) == 1:
            ready.append(equation_id)

    while ready:
        equation_id = ready.popleft()
        unknown = unknown_sets[equation_id]
        if len(unknown) != 1:
            continue
        index = next(iter(unknown))
        if index in recovered:
            continue
        payload = bytes(residuals[equation_id])
        recovered[index] = payload

        for dependent_id in incident[index]:
            dependent = unknown_sets[dependent_id]
            if index not in dependent:
                continue
            dependent.remove(index)
            residual = residuals[dependent_id]
            for offset in range(width):
                residual[offset] ^= payload[offset]
            if len(dependent) == 1:
                ready.append(dependent_id)

    return {index: recovered[index] for index in range(total) if index in recovered}

"""Deterministic sparse XOR redundancy helpers.

Legacy random LT-style symbols remain available for historical profiles. Research profiles can
use balanced interleaved symbols: deterministic permutation layers that give every data chunk
predictable sparse-check coverage at the same nominal redundancy budget.
"""

from __future__ import annotations

import math
import random
from collections import deque
from dataclasses import dataclass

from .ecc import xor_bytes

FOUNTAIN_LAYOUTS = {"random", "interleaved"}


@dataclass(frozen=True)
class FountainSymbol:
    seed: int
    indexes: tuple[int, ...]
    payload: bytes


def _coprime_multiplier(total: int, layer: int) -> int:
    if total <= 1:
        return 1
    candidate = max(3, 2 * layer + 1)
    while math.gcd(candidate, total) != 1:
        candidate += 2
    return candidate


def interleaved_indexes_for_seed(
    total: int,
    seed: int,
    *,
    group_size: int,
    first_seed: int = 1,
) -> tuple[int, ...]:
    """Map one seed to a balanced sparse-check group.

    Layer zero is supplied by OligoArk's ordinary contiguous XOR parity. Fountain seeds start
    at permutation layer one, so a hybrid profile can combine one contiguous layer with two
    interleaved layers without duplicating the same check groups.
    """
    if total <= 0:
        raise ValueError("total must be positive")
    if group_size < 2:
        raise ValueError("group_size must be at least 2")
    offset = seed - first_seed
    if offset < 0:
        raise ValueError("seed is below first_seed")

    groups_per_layer = math.ceil(total / group_size)
    layer = 1 + offset // groups_per_layer
    group = offset % groups_per_layer
    multiplier = _coprime_multiplier(total, layer)
    shift = (layer * 0x9E3779B1) % total

    start = group * group_size
    stop = min(start + group_size, total)
    return tuple(
        sorted((multiplier * position + shift) % total for position in range(start, stop))
    )


def symbol_count(
    total: int,
    redundancy: float,
    *,
    max_degree: int,
    layout: str,
) -> int:
    """Return the deterministic symbol count for a requested redundancy budget."""
    if total <= 0:
        raise ValueError("total must be positive")
    if redundancy <= 0:
        return 0
    if layout == "interleaved":
        layers = max(1, round(redundancy * max_degree))
        return layers * math.ceil(total / max_degree)
    if layout != "random":
        raise ValueError(f"unknown fountain layout: {layout}")
    return max(1, math.ceil(total * redundancy))


def indexes_for_seed(
    total: int,
    seed: int,
    max_degree: int = 4,
    *,
    layout: str = "random",
    first_seed: int = 1,
) -> tuple[int, ...]:
    """Return the deterministic source-chunk set represented by one fountain seed."""
    if total <= 0:
        raise ValueError("total must be positive")
    if max_degree < 1:
        raise ValueError("max_degree must be positive")
    if layout == "interleaved":
        return interleaved_indexes_for_seed(
            total,
            seed,
            group_size=max_degree,
            first_seed=first_seed,
        )
    if layout != "random":
        raise ValueError(f"unknown fountain layout: {layout}")
    rng = random.Random(seed)
    degree = min(total, 1 + rng.randrange(max(1, min(max_degree, total))))
    return tuple(sorted(rng.sample(range(total), degree)))


def make_symbols(
    chunks: list[bytes],
    count: int,
    width: int,
    seed: int = 1,
    max_degree: int = 4,
    *,
    layout: str = "random",
) -> list[FountainSymbol]:
    if not chunks:
        raise ValueError("chunks must not be empty")
    if count < 0:
        raise ValueError("count must be non-negative")
    if width <= 0:
        raise ValueError("width must be positive")
    if layout not in FOUNTAIN_LAYOUTS:
        raise ValueError(f"unknown fountain layout: {layout}")

    symbols: list[FountainSymbol] = []
    for offset in range(count):
        symbol_seed = seed + offset
        indexes = indexes_for_seed(
            len(chunks),
            symbol_seed,
            max_degree=max_degree,
            layout=layout,
            first_seed=seed,
        )
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

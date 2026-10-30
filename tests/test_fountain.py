from __future__ import annotations

import random

from oligoark.archive import archive_bytes, recover_bytes
from oligoark.fountain import (
    indexes_for_seed,
    make_symbols,
    peel_decode,
    symbol_count,
)
from oligoark.profiles import physical_strand_profile


def test_fountain_style_symbol_can_recover_missing_chunk() -> None:
    chunks = [bytes([i]) * 8 for i in range(5)]
    symbols = make_symbols(chunks, count=30, width=8, seed=10)
    known = {i: chunk for i, chunk in enumerate(chunks) if i != 2}
    recovered = peel_decode(known, symbols, total=5, width=8)
    assert recovered[2] == chunks[2]


def test_interleaved_layout_covers_every_chunk_once_per_layer() -> None:
    total = 125
    degree = 12
    first_seed = 7
    count = symbol_count(
        total,
        1.0 / 6.0,
        max_degree=degree,
        layout="interleaved",
    )
    coverage = [0] * total
    for offset in range(count):
        indexes = indexes_for_seed(
            total,
            first_seed + offset,
            max_degree=degree,
            layout="interleaved",
            first_seed=first_seed,
        )
        assert 1 <= len(indexes) <= degree
        assert len(indexes) == len(set(indexes))
        for index in indexes:
            coverage[index] += 1
    assert set(coverage) == {2}


def test_v3_hybrid_recovers_deterministic_five_percent_strand_dropout() -> None:
    profile = physical_strand_profile("oligoark-152-compact-v3").with_scheme("hybrid")
    config = profile.to_archive_config()
    payload = bytes(range(256)) * 20
    archive = archive_bytes(payload, config)

    for seed in range(10):
        rng = random.Random(20_260_000 + seed)
        drop_count = round(len(archive.strands) * 0.05)
        dropped = set(rng.sample(range(len(archive.strands)), drop_count))
        reads = [
            sequence
            for index, sequence in enumerate(archive.strands)
            if index not in dropped
        ]
        assert recover_bytes(archive, reads) == payload


def _reference_peel_with_satisfied_equations(
    known: dict[int, bytes], symbols: list, total: int, width: int
) -> dict[int, bytes]:
    """Independent legacy peel including equations with no unknown indexes."""
    from collections import deque

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
                for offset, value in enumerate(recovered[index]):
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


def test_fountain_peel_randomized_parity_with_legacy_equations() -> None:
    for seed in range(25):
        rng = random.Random(2026 + seed)
        total = rng.randrange(3, 50)
        width = rng.randrange(1, 48)
        chunks = [rng.randbytes(width) for _ in range(total)]
        layout = "interleaved" if seed % 2 else "random"
        symbols = make_symbols(
            chunks, count=total * 2, width=width, seed=seed + 1,
            max_degree=4, layout=layout,
        )
        missing = set(rng.sample(range(total), rng.randrange(total + 1)))
        known = {index: chunk for index, chunk in enumerate(chunks) if index not in missing}
        before = dict(known)
        expected = _reference_peel_with_satisfied_equations(known, symbols, total, width)
        actual = peel_decode(known, symbols, total, width)
        assert actual == expected
        assert known == before
        assert all(actual[index] == chunks[index] for index in actual)

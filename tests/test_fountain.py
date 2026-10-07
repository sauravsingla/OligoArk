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

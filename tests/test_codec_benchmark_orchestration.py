from __future__ import annotations

import importlib.util
import random
from pathlib import Path
from types import ModuleType


def _benchmark_module() -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "benchmarks" / "run_dna_fountain_baseline.py"
    spec = importlib.util.spec_from_file_location("run_dna_fountain_baseline", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_BENCHMARK = _benchmark_module()
_bernoulli_positions = _BENCHMARK._bernoulli_positions
_channel_sequences = _BENCHMARK._channel_sequences
_mutate_sequence = _BENCHMARK._mutate_sequence


def test_clean_channel_reuses_immutable_sequences() -> None:
    sequences = ("ACGT" * 10, "TGCA" * 10)

    reads = _channel_sequences(sequences, "clean", 20260000)

    assert reads is sequences


def test_dropout_only_channel_is_deterministic() -> None:
    sequences = tuple("ACGT" * 10 for _ in range(100))

    first = _channel_sequences(sequences, "dropout-5", 20260003)
    second = _channel_sequences(sequences, "dropout-5", 20260003)

    assert first == second
    assert len(first) == 95


def test_sparse_bernoulli_sampler_is_deterministic_and_bounded() -> None:
    first = _bernoulli_positions(10_000, 0.001, random.Random(1234))
    second = _bernoulli_positions(10_000, 0.001, random.Random(1234))

    assert first == second
    assert all(0 <= position < 10_000 for position in first)
    assert 1 <= len(first) <= 40


def test_sparse_mutation_preserves_dna_alphabet_and_is_deterministic() -> None:
    sequence = "ACGT" * 100
    kwargs = {
        "substitution_rate": 0.01,
        "insertion_rate": 0.005,
        "deletion_rate": 0.005,
    }

    first = _mutate_sequence(sequence, rng=random.Random(99), **kwargs)
    second = _mutate_sequence(sequence, rng=random.Random(99), **kwargs)

    assert first == second
    assert set(first) <= {"A", "C", "G", "T"}
    assert first != sequence

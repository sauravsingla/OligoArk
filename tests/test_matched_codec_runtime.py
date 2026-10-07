from __future__ import annotations

from benchmarks.run_dna_fountain_baseline import CONDITIONS, _channel_sequences


def test_clean_channel_reuses_immutable_encoded_artifact() -> None:
    sequences = ("ACGT" * 38, "TGCA" * 38)

    selected = _channel_sequences(sequences, "clean", 20260000)

    assert selected is sequences


def test_dropout_channel_is_deterministic_and_does_not_mutate_survivors() -> None:
    sequences = tuple(("ACGT" * 38) for _ in range(1000))

    first = _channel_sequences(sequences, "dropout-5", 20260001)
    second = _channel_sequences(sequences, "dropout-5", 20260001)

    assert first == second
    assert len(first) == 950
    assert set(first) == {"ACGT" * 38}


def test_sparse_noise_channel_is_deterministic_and_preserves_dna_alphabet() -> None:
    sequences = tuple(("ACGT" * 38) for _ in range(256))

    first = _channel_sequences(sequences, "mixed", 20260002)
    second = _channel_sequences(sequences, "mixed", 20260002)

    assert first == second
    assert len(first) == len(sequences) - round(
        len(sequences) * CONDITIONS["mixed"]["dropout_rate"]
    )
    assert all(set(sequence) <= {"A", "C", "G", "T"} for sequence in first)
    assert any(sequence != "ACGT" * 38 for sequence in first)

"""Low-level binary/DNA conversion and sequence quality helpers."""

from __future__ import annotations

from dataclasses import dataclass

BITS_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}
BASE_TO_BITS = {v: k for k, v in BITS_TO_BASE.items()}


def bytes_to_dna(data: bytes) -> str:
    """Encode bytes into a reversible 2-bit DNA alphabet representation."""
    out: list[str] = []
    for value in data:
        out.extend(
            (
                BITS_TO_BASE[(value >> 6) & 0b11],
                BITS_TO_BASE[(value >> 4) & 0b11],
                BITS_TO_BASE[(value >> 2) & 0b11],
                BITS_TO_BASE[value & 0b11],
            )
        )
    return "".join(out)


def dna_to_bytes(sequence: str) -> bytes:
    """Decode a 2-bit DNA string back to bytes."""
    sequence = sequence.strip().upper()
    if len(sequence) % 4:
        raise ValueError("DNA sequence length must be divisible by 4")
    try:
        vals = [BASE_TO_BITS[b] for b in sequence]
    except KeyError as exc:
        raise ValueError(f"Invalid DNA base: {exc.args[0]!r}") from exc
    return bytes(
        (vals[i] << 6) | (vals[i + 1] << 4) | (vals[i + 2] << 2) | vals[i + 3]
        for i in range(0, len(vals), 4)
    )


@dataclass(frozen=True)
class SequenceMetrics:
    gc_fraction: float
    max_homopolymer: int


def sequence_metrics(sequence: str) -> SequenceMetrics:
    if not sequence:
        return SequenceMetrics(0.0, 0)
    sequence = sequence.upper()
    gc = sum(base in {"G", "C"} for base in sequence) / len(sequence)
    max_run = run = 1
    for prev, cur in zip(sequence, sequence[1:], strict=False):
        run = run + 1 if prev == cur else 1
        max_run = max(max_run, run)
    return SequenceMetrics(gc_fraction=gc, max_homopolymer=max_run)


def quality_penalty(sequence: str, gc_target: float = 0.5, max_homopolymer: int = 4) -> float:
    """A transparent software heuristic; not a wet-lab fitness claim."""
    m = sequence_metrics(sequence)
    gc_penalty = abs(m.gc_fraction - gc_target)
    homopolymer_penalty = max(0, m.max_homopolymer - max_homopolymer) * 0.05
    return gc_penalty + homopolymer_penalty

"""Low-level binary/DNA conversion and sequence quality/constraint helpers."""

from __future__ import annotations

from dataclasses import asdict, dataclass

BITS_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}
BASE_TO_BITS = {v: k for k, v in BITS_TO_BASE.items()}


class SequenceConstraintError(ValueError):
    """Raised when no deterministic encoding candidate satisfies hard sequence constraints."""


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
        vals = [BASE_TO_BITS[base] for base in sequence]
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

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class SequenceConstraints:
    """Hard DNA-like sequence constraints used by the software encoder.

    These defaults are research settings, not biochemical claims about a specific synthesis
    or sequencing platform.
    """

    min_gc_fraction: float = 0.35
    max_gc_fraction: float = 0.65
    max_homopolymer: int = 4

    def validate(self) -> None:
        if not 0.0 <= self.min_gc_fraction <= 1.0:
            raise ValueError("min_gc_fraction must be between 0 and 1")
        if not 0.0 <= self.max_gc_fraction <= 1.0:
            raise ValueError("max_gc_fraction must be between 0 and 1")
        if self.min_gc_fraction > self.max_gc_fraction:
            raise ValueError("min_gc_fraction must be <= max_gc_fraction")
        if self.max_homopolymer < 1:
            raise ValueError("max_homopolymer must be positive")

    def accepts(self, sequence: str) -> bool:
        self.validate()
        metrics = sequence_metrics(sequence)
        return (
            self.min_gc_fraction <= metrics.gc_fraction <= self.max_gc_fraction
            and metrics.max_homopolymer <= self.max_homopolymer
        )

    def violations(self, sequence: str) -> tuple[str, ...]:
        self.validate()
        metrics = sequence_metrics(sequence)
        failures: list[str] = []
        if metrics.gc_fraction < self.min_gc_fraction:
            failures.append(
                f"gc_fraction={metrics.gc_fraction:.4f} < min={self.min_gc_fraction:.4f}"
            )
        if metrics.gc_fraction > self.max_gc_fraction:
            failures.append(
                f"gc_fraction={metrics.gc_fraction:.4f} > max={self.max_gc_fraction:.4f}"
            )
        if metrics.max_homopolymer > self.max_homopolymer:
            failures.append(
                f"max_homopolymer={metrics.max_homopolymer} > limit={self.max_homopolymer}"
            )
        return tuple(failures)


def sequence_metrics(sequence: str) -> SequenceMetrics:
    """Measure GC fraction and longest same-base run."""
    if not sequence:
        return SequenceMetrics(0.0, 0)
    sequence = sequence.upper()
    invalid = set(sequence) - {"A", "C", "G", "T"}
    if invalid:
        raise ValueError(f"Invalid DNA base(s): {sorted(invalid)}")
    gc = sum(base in {"G", "C"} for base in sequence) / len(sequence)
    max_run = run = 1
    for previous, current in zip(sequence, sequence[1:], strict=False):
        run = run + 1 if previous == current else 1
        max_run = max(max_run, run)
    return SequenceMetrics(gc_fraction=gc, max_homopolymer=max_run)


def quality_penalty(
    sequence: str,
    gc_target: float = 0.5,
    max_homopolymer: int = 4,
) -> float:
    """Return a transparent soft score used only to break ties among valid candidates."""
    if not 0.0 <= gc_target <= 1.0:
        raise ValueError("gc_target must be between 0 and 1")
    if max_homopolymer < 1:
        raise ValueError("max_homopolymer must be positive")
    metrics = sequence_metrics(sequence)
    gc_penalty = abs(metrics.gc_fraction - gc_target)
    homopolymer_penalty = max(0, metrics.max_homopolymer - max_homopolymer) * 0.05
    return gc_penalty + homopolymer_penalty

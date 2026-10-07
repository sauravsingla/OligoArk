"""Independent reference baselines for DNA-storage experiments.

The DNA Fountain implementation here is a clean-room research baseline derived from
Erlich & Zielinski (Science, 2017): robust-soliton LT droplets, 32-bit seeds, XOR payloads,
Reed-Solomon protection, 2-bit DNA mapping, and GC/homopolymer screening. It does not copy
or claim bit compatibility with the historical GPL-licensed implementation.
"""

from __future__ import annotations

import bisect
import hashlib
import math
import random
from dataclasses import asdict, dataclass
from functools import lru_cache

from .dna import SequenceConstraints, bytes_to_dna, dna_to_bytes
from .ecc import rs_decode, rs_encode, xor_bytes


@dataclass(frozen=True)
class DnaFountainBaselineConfig:
    chunk_size: int = 32
    redundancy: float = 0.25
    rs_nsym: int = 2
    c: float = 0.025
    delta: float = 0.001
    gc_tolerance: float = 0.05
    max_homopolymer: int = 3
    first_seed: int = 1
    max_attempt_factor: int = 100

    def validate(self) -> None:
        if self.chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
        if not 0.0 <= self.redundancy <= 5.0:
            raise ValueError("redundancy must be between 0 and 5")
        if not 0 <= self.rs_nsym <= 64:
            raise ValueError("rs_nsym must be between 0 and 64")
        if self.chunk_size + 4 + self.rs_nsym > 255:
            raise ValueError("seed + payload + rs_nsym must fit one RS codeword")
        if self.c <= 0 or not 0 < self.delta < 1:
            raise ValueError("c must be positive and delta must be between 0 and 1")
        if not 0 <= self.gc_tolerance <= 0.5:
            raise ValueError("gc_tolerance must be between 0 and 0.5")
        if self.max_homopolymer < 1:
            raise ValueError("max_homopolymer must be positive")
        if self.max_attempt_factor < 1:
            raise ValueError("max_attempt_factor must be positive")


@dataclass(frozen=True)
class DnaFountainBaselineArchive:
    original_size: int
    sha256: str
    chunk_count: int
    sequences: tuple[str, ...]
    attempts: int
    config: DnaFountainBaselineConfig

    @property
    def encoded_nucleotides(self) -> int:
        return sum(map(len, self.sequences))

    @property
    def logical_bits_per_nucleotide(self) -> float:
        if not self.encoded_nucleotides:
            return 0.0
        return self.original_size * 8 / self.encoded_nucleotides

    def to_dict(self) -> dict[str, object]:
        return {
            "original_size": self.original_size,
            "sha256": self.sha256,
            "chunk_count": self.chunk_count,
            "sequence_count": len(self.sequences),
            "attempts": self.attempts,
            "encoded_nucleotides": self.encoded_nucleotides,
            "logical_bits_per_nucleotide": round(self.logical_bits_per_nucleotide, 6),
            "config": asdict(self.config),
        }


@lru_cache(maxsize=32)
def _robust_soliton_cdf(k: int, c: float, delta: float) -> tuple[float, ...]:
    if k <= 0:
        raise ValueError("k must be positive")
    if k == 1:
        return (1.0,)
    s = max(c * math.log(k / delta) * math.sqrt(k), 1e-12)
    pivot = max(1, min(k, int(math.floor(k / s))))
    weights: list[float] = []
    for degree in range(1, k + 1):
        rho = 1.0 / k if degree == 1 else 1.0 / (degree * (degree - 1))
        if degree < pivot:
            tau = s / (k * degree)
        elif degree == pivot:
            tau = (s / k) * math.log(s / delta)
        else:
            tau = 0.0
        weights.append(max(0.0, rho + tau))
    total = sum(weights)
    running = 0.0
    cdf: list[float] = []
    for weight in weights:
        running += weight / total
        cdf.append(min(1.0, running))
    cdf[-1] = 1.0
    return tuple(cdf)


def _droplet_indexes(k: int, seed: int, c: float, delta: float) -> tuple[int, ...]:
    rng = random.Random(seed)
    degree = bisect.bisect_left(_robust_soliton_cdf(k, c, delta), rng.random()) + 1
    return tuple(sorted(rng.sample(range(k), min(k, max(1, degree)))))


def encode_dna_fountain_baseline(
    data: bytes,
    config: DnaFountainBaselineConfig | None = None,
) -> DnaFountainBaselineArchive:
    resolved = config or DnaFountainBaselineConfig()
    resolved.validate()
    chunks = [data[i : i + resolved.chunk_size] for i in range(0, len(data), resolved.chunk_size)]
    chunks = chunks or [b""]
    target = max(1, math.ceil(len(chunks) * (1.0 + resolved.redundancy)))
    constraints = SequenceConstraints(
        min_gc_fraction=0.5 - resolved.gc_tolerance,
        max_gc_fraction=0.5 + resolved.gc_tolerance,
        max_homopolymer=resolved.max_homopolymer,
    )
    accepted: list[str] = []
    seed = resolved.first_seed & 0xFFFFFFFF
    attempts = 0
    max_attempts = target * resolved.max_attempt_factor
    while len(accepted) < target and attempts < max_attempts:
        indexes = _droplet_indexes(len(chunks), seed, resolved.c, resolved.delta)
        payload = xor_bytes([chunks[index] for index in indexes], resolved.chunk_size)
        sequence = bytes_to_dna(
            rs_encode(seed.to_bytes(4, "big") + payload, resolved.rs_nsym)
        )
        attempts += 1
        # A full-width Weyl step avoids the long zero-byte prefixes produced by small
        # sequential seeds while keeping seed generation deterministic and reversible.
        seed = (seed + 0x9E3779B9) & 0xFFFFFFFF
        if constraints.accepts(sequence):
            accepted.append(sequence)
    if len(accepted) < target:
        raise ValueError(
            f"DNA Fountain baseline screening budget exhausted ({len(accepted)}/{target})"
        )
    return DnaFountainBaselineArchive(
        original_size=len(data),
        sha256=hashlib.sha256(data).hexdigest(),
        chunk_count=len(chunks),
        sequences=tuple(accepted),
        attempts=attempts,
        config=resolved,
    )


def _mask_indexes(mask: int) -> tuple[int, ...]:
    indexes: list[int] = []
    while mask:
        bit = mask & -mask
        indexes.append(bit.bit_length() - 1)
        mask ^= bit
    return tuple(indexes)


def _gaussian_recover(
    equations: list[tuple[set[int], bytes]],
    known: dict[int, bytes],
    *,
    total: int,
    width: int,
) -> dict[int, bytes]:
    """Solve residual XOR equations exactly when LT peeling stalls."""
    pivots: dict[int, tuple[int, bytes]] = {}
    for indexes, payload in equations:
        mask = 0
        residual_parts = [payload]
        for index in indexes:
            if index in known:
                residual_parts.append(known[index])
            else:
                mask |= 1 << index
        residual = xor_bytes(residual_parts, width)
        while mask:
            pivot = (mask & -mask).bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (mask, residual)
                break
            pivot_mask, pivot_payload = pivots[pivot]
            mask ^= pivot_mask
            residual = xor_bytes([residual, pivot_payload], width)

    recovered = dict(known)
    for pivot in sorted(pivots, reverse=True):
        mask, payload = pivots[pivot]
        residual_parts = [payload]
        unresolved = False
        for index in _mask_indexes(mask & ~(1 << pivot)):
            if index not in recovered:
                unresolved = True
                break
            residual_parts.append(recovered[index])
        if not unresolved:
            recovered[pivot] = xor_bytes(residual_parts, width)
    return {index: recovered[index] for index in range(total) if index in recovered}


def decode_dna_fountain_baseline(
    archive: DnaFountainBaselineArchive,
    sequences: tuple[str, ...] | list[str] | None = None,
) -> bytes:
    selected = archive.sequences if sequences is None else tuple(sequences)
    equations: list[tuple[set[int], bytes]] = []
    for sequence in selected:
        try:
            packet = rs_decode(dna_to_bytes(sequence), archive.config.rs_nsym)
        except ValueError:
            continue
        if len(packet) != 4 + archive.config.chunk_size:
            continue
        seed = int.from_bytes(packet[:4], "big")
        droplet_indexes = _droplet_indexes(
            archive.chunk_count,
            seed,
            archive.config.c,
            archive.config.delta,
        )
        equations.append((set(droplet_indexes), packet[4:]))

    known: dict[int, bytes] = {}
    pending = equations
    changed = True
    while changed:
        changed = False
        next_pending: list[tuple[set[int], bytes]] = []
        for unknown_indexes, payload in pending:
            unknown = set(unknown_indexes)
            parts = [payload]
            for index in tuple(unknown):
                if index in known:
                    parts.append(known[index])
                    unknown.remove(index)
            residual = xor_bytes(parts, archive.config.chunk_size)
            if len(unknown) == 1:
                index = next(iter(unknown))
                if index not in known:
                    known[index] = residual
                    changed = True
            elif unknown:
                next_pending.append((unknown, residual))
        pending = next_pending

    missing = [index for index in range(archive.chunk_count) if index not in known]
    if missing:
        known = _gaussian_recover(
            equations,
            known,
            total=archive.chunk_count,
            width=archive.config.chunk_size,
        )
        missing = [index for index in range(archive.chunk_count) if index not in known]
    if missing:
        raise ValueError(
            "DNA Fountain baseline is not recoverable; "
            f"missing chunks: {missing[:20]}"
        )
    recovered = b"".join(
        known[index] for index in range(archive.chunk_count)
    )[: archive.original_size]
    if hashlib.sha256(recovered).hexdigest() != archive.sha256:
        raise ValueError("DNA Fountain baseline failed SHA-256 verification")
    return recovered

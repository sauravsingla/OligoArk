"""Independent reference baselines for DNA-storage experiments.

The DNA Fountain implementation here is a clean-room research baseline derived from
Erlich & Zielinski (Science, 2017): robust-soliton LT droplets, 32-bit seeds, XOR payloads,
Reed-Solomon protection, 2-bit DNA mapping, and GC/homopolymer screening. A second
Goldman-style rotating-ternary reference maps trits onto non-repeating DNA bases and uses a
simple XOR erasure layer so dropout can be measured with the same exact-recovery gate.
Neither implementation copies or claims bit compatibility with historical source code.
"""

from __future__ import annotations

import bisect
import hashlib
import math
import random
import struct
import zlib
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


_ROTATING_HEADER = struct.Struct(">BIIH")
_ROTATING_FLAG_PARITY = 0x01
_ROTATING_BASES = "ACGT"


@dataclass(frozen=True)
class RotatingTernaryBaselineConfig:
    """Goldman-style rotating ternary reference with optional XOR erasure parity."""

    chunk_size: int = 12
    redundancy: float = 0.25
    max_strand_nt: int = 152

    def validate(self) -> None:
        if self.chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
        if not 0.0 <= self.redundancy <= 1.0:
            raise ValueError("redundancy must be between 0 and 1")
        if self.max_strand_nt < 32:
            raise ValueError("max_strand_nt is too short")
        probe = _rotating_frame(
            b"x" * self.chunk_size,
            index=0,
            total=1,
            parity=False,
        )
        if len(probe) > self.max_strand_nt:
            raise ValueError(
                "rotating ternary frame exceeds configured strand limit "
                f"({len(probe)} > {self.max_strand_nt})"
            )

    @property
    def parity_group_size(self) -> int | None:
        if self.redundancy <= 0:
            return None
        return max(2, round(1.0 / self.redundancy))


@dataclass(frozen=True)
class RotatingTernaryBaselineArchive:
    original_size: int
    sha256: str
    data_count: int
    sequences: tuple[str, ...]
    config: RotatingTernaryBaselineConfig

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
            "data_count": self.data_count,
            "sequence_count": len(self.sequences),
            "encoded_nucleotides": self.encoded_nucleotides,
            "logical_bits_per_nucleotide": round(self.logical_bits_per_nucleotide, 6),
            "config": asdict(self.config),
        }


def _fixed_trits(value: int, width: int) -> list[int]:
    if value < 0:
        raise ValueError("value must be non-negative")
    out = [0] * width
    for index in range(width - 1, -1, -1):
        out[index] = value % 3
        value //= 3
    if value:
        raise ValueError("value does not fit requested ternary width")
    return out


def _bytes_to_trits(data: bytes) -> list[int]:
    out: list[int] = []
    pair_bytes = len(data) - (len(data) % 2)
    for index in range(0, pair_bytes, 2):
        out.extend(_fixed_trits(int.from_bytes(data[index : index + 2], "big"), 11))
    if len(data) % 2:
        out.extend(_fixed_trits(data[-1], 6))
    return out


def _trits_to_bytes(trits: list[int]) -> bytes:
    remainder = len(trits) % 11
    if remainder not in {0, 6}:
        raise ValueError("invalid rotating ternary length")
    out = bytearray()
    pair_limit = len(trits) - remainder
    for start in range(0, pair_limit, 11):
        value = 0
        for trit in trits[start : start + 11]:
            value = value * 3 + trit
        if value > 0xFFFF:
            raise ValueError("invalid rotating ternary pair")
        out.extend(value.to_bytes(2, "big"))
    if remainder:
        value = 0
        for trit in trits[pair_limit:]:
            value = value * 3 + trit
        if value > 0xFF:
            raise ValueError("invalid rotating ternary byte")
        out.append(value)
    return bytes(out)


def _rotating_encode(data: bytes) -> str:
    previous = "A"
    sequence: list[str] = []
    for trit in _bytes_to_trits(data):
        choices = [base for base in _ROTATING_BASES if base != previous]
        base = choices[trit]
        sequence.append(base)
        previous = base
    return "".join(sequence)


def _rotating_decode(sequence: str) -> bytes:
    previous = "A"
    trits: list[int] = []
    for raw_base in sequence:
        base = raw_base.upper()
        choices = [candidate for candidate in _ROTATING_BASES if candidate != previous]
        try:
            trits.append(choices.index(base))
        except ValueError as exc:
            raise ValueError("invalid rotating ternary sequence") from exc
        previous = base
    return _trits_to_bytes(trits)


def _rotating_frame(
    payload: bytes,
    *,
    index: int,
    total: int,
    parity: bool,
) -> str:
    flags = _ROTATING_FLAG_PARITY if parity else 0
    header = _ROTATING_HEADER.pack(flags, index, total, len(payload))
    crc = (zlib.crc32(payload) & 0xFFFFFFFF).to_bytes(4, "big")
    return _rotating_encode(header + payload + crc)


def _decode_rotating_frame(sequence: str) -> tuple[bool, int, int, bytes]:
    raw = _rotating_decode(sequence)
    if len(raw) < _ROTATING_HEADER.size + 4:
        raise ValueError("rotating ternary frame is too short")
    flags, index, total, payload_len = _ROTATING_HEADER.unpack(raw[: _ROTATING_HEADER.size])
    end = _ROTATING_HEADER.size + payload_len
    if len(raw) != end + 4:
        raise ValueError("rotating ternary frame length mismatch")
    payload = raw[_ROTATING_HEADER.size : end]
    crc = int.from_bytes(raw[end:], "big")
    if (zlib.crc32(payload) & 0xFFFFFFFF) != crc:
        raise ValueError("rotating ternary CRC check failed")
    return bool(flags & _ROTATING_FLAG_PARITY), index, total, payload


def encode_rotating_ternary_baseline(
    data: bytes,
    config: RotatingTernaryBaselineConfig | None = None,
) -> RotatingTernaryBaselineArchive:
    resolved = config or RotatingTernaryBaselineConfig()
    resolved.validate()
    chunks = [
        data[index : index + resolved.chunk_size]
        for index in range(0, len(data), resolved.chunk_size)
    ] or [b""]
    total = len(chunks)
    sequences = [
        _rotating_frame(chunk, index=index, total=total, parity=False)
        for index, chunk in enumerate(chunks)
    ]
    group_size = resolved.parity_group_size
    if group_size is not None:
        for group_index, start in enumerate(range(0, total, group_size)):
            parity_payload = xor_bytes(
                chunks[start : start + group_size],
                resolved.chunk_size,
            )
            sequences.append(
                _rotating_frame(
                    parity_payload,
                    index=group_index,
                    total=total,
                    parity=True,
                )
            )
    if max(map(len, sequences)) > resolved.max_strand_nt:
        raise ValueError("rotating ternary archive exceeded configured strand limit")
    return RotatingTernaryBaselineArchive(
        original_size=len(data),
        sha256=hashlib.sha256(data).hexdigest(),
        data_count=total,
        sequences=tuple(sequences),
        config=resolved,
    )


def decode_rotating_ternary_baseline(
    archive: RotatingTernaryBaselineArchive,
    sequences: tuple[str, ...] | list[str] | None = None,
) -> bytes:
    selected = archive.sequences if sequences is None else tuple(sequences)
    known: dict[int, bytes] = {}
    parity: dict[int, bytes] = {}
    for sequence in selected:
        try:
            is_parity, index, total, payload = _decode_rotating_frame(sequence)
        except ValueError:
            continue
        if total != archive.data_count:
            continue
        if is_parity:
            parity[index] = payload
        elif 0 <= index < archive.data_count:
            known[index] = payload

    group_size = archive.config.parity_group_size
    if group_size is not None:
        for group_index, start in enumerate(range(0, archive.data_count, group_size)):
            end = min(start + group_size, archive.data_count)
            missing = [index for index in range(start, end) if index not in known]
            if len(missing) != 1 or group_index not in parity:
                continue
            parts = [parity[group_index]] + [
                known[index] for index in range(start, end) if index in known
            ]
            known[missing[0]] = xor_bytes(parts, archive.config.chunk_size)

    missing = [index for index in range(archive.data_count) if index not in known]
    if missing:
        raise ValueError(
            "rotating ternary baseline is not recoverable; "
            f"missing chunks: {missing[:20]}"
        )
    recovered = b"".join(
        known[index] for index in range(archive.data_count)
    )[: archive.original_size]
    if hashlib.sha256(recovered).hexdigest() != archive.sha256:
        raise ValueError("rotating ternary baseline failed SHA-256 verification")
    return recovered

"""Self-describing, error-corrected strand framing for OligoArk."""

from __future__ import annotations

import struct
import zlib
from dataclasses import dataclass

from .dna import (
    SequenceConstraintError,
    SequenceConstraints,
    bytes_to_dna,
    dna_to_bytes,
    quality_penalty,
)
from .ecc import ECCDecodeError, rs_decode, rs_encode

MAGIC = b"OA"
VERSION = 1
FLAG_PARITY = 0x01
FLAG_FOUNTAIN = 0x02
_HEADER = struct.Struct(">2sBBIIHI")
_LEGACY_MASKS = (0x00, 0x55, 0xAA, 0xFF)


@dataclass(frozen=True)
class DecodedFrame:
    index: int
    total_data: int
    payload: bytes
    is_parity: bool
    is_fountain: bool = False

    @property
    def kind(self) -> str:
        if self.is_fountain:
            return "fountain"
        if self.is_parity:
            return "parity"
        return "data"


def frame_overhead_bytes(rs_nsym: int) -> int:
    """Return fixed per-strand byte overhead before payload bytes."""
    return 1 + _HEADER.size + rs_nsym


def _mask(data: bytes, mask_id: int) -> bytes:
    """Apply a reversible deterministic whitening mask.

    IDs 0..3 preserve the v0.1-v0.3 constant-mask format. IDs 4..255 use a deterministic
    byte stream so the encoder can search a larger sequence-constraint space without changing
    the one-byte mask identifier stored in the strand.
    """
    if not 0 <= mask_id <= 255:
        raise ValueError("mask_id must fit in one byte")
    if mask_id < len(_LEGACY_MASKS):
        value = _LEGACY_MASKS[mask_id]
        return bytes(byte ^ value for byte in data)

    state = (0x9E3779B9 ^ (mask_id * 0x45D9F3B)) & 0xFFFFFFFF
    out = bytearray(len(data))
    for index, byte in enumerate(data):
        state ^= (state << 13) & 0xFFFFFFFF
        state ^= state >> 17
        state ^= (state << 5) & 0xFFFFFFFF
        out[index] = byte ^ (state & 0xFF)
    return bytes(out)


def encode_frame(
    payload: bytes,
    *,
    index: int,
    total_data: int,
    is_parity: bool,
    rs_nsym: int,
    adaptive_masks: bool,
    is_fountain: bool = False,
    sequence_constraints: SequenceConstraints | None = None,
    mask_search_limit: int = 64,
) -> str:
    """Encode one protected strand while enforcing optional hard sequence constraints."""
    if not 1 <= mask_search_limit <= 256:
        raise ValueError("mask_search_limit must be between 1 and 256")
    constraints = sequence_constraints or SequenceConstraints()
    constraints.validate()

    flags = 0
    if is_parity:
        flags |= FLAG_PARITY
    if is_fountain:
        flags |= FLAG_FOUNTAIN

    header = _HEADER.pack(
        MAGIC,
        VERSION,
        flags,
        index,
        total_data,
        len(payload),
        zlib.crc32(payload) & 0xFFFFFFFF,
    )
    protected = rs_encode(header + payload, rs_nsym)
    candidate_ids = range(mask_search_limit) if adaptive_masks else range(1)
    valid: list[tuple[float, str]] = []

    for mask_id in candidate_ids:
        dna = bytes_to_dna(bytes([mask_id]) + _mask(protected, mask_id))
        if constraints.accepts(dna):
            target = (constraints.min_gc_fraction + constraints.max_gc_fraction) / 2.0
            score = quality_penalty(
                dna,
                gc_target=target,
                max_homopolymer=constraints.max_homopolymer,
            )
            valid.append((score, dna))

    if not valid:
        raise SequenceConstraintError(
            "No deterministic mask candidate satisfied configured GC/homopolymer constraints"
        )
    return min(valid, key=lambda item: item[0])[1]


def _decode_with_mask(raw: bytes, mask_id: int, rs_nsym: int) -> DecodedFrame:
    protected = _mask(raw[1:], mask_id)
    try:
        inner = rs_decode(protected, rs_nsym)
    except ECCDecodeError as exc:
        raise ValueError("Reed-Solomon recovery failed") from exc
    if len(inner) < _HEADER.size:
        raise ValueError("Decoded frame is shorter than the OligoArk header")
    magic, version, flags, index, total, raw_len, crc = _HEADER.unpack(inner[: _HEADER.size])
    if magic != MAGIC or version != VERSION:
        raise ValueError("Not a supported OligoArk strand")
    payload = inner[_HEADER.size : _HEADER.size + raw_len]
    if len(payload) != raw_len:
        raise ValueError("Decoded payload length is inconsistent with frame header")
    if (zlib.crc32(payload) & 0xFFFFFFFF) != crc:
        raise ValueError("Payload CRC check failed")
    return DecodedFrame(
        index=index,
        total_data=total,
        payload=payload,
        is_parity=bool(flags & FLAG_PARITY),
        is_fountain=bool(flags & FLAG_FOUNTAIN),
    )


def decode_frame(
    sequence: str,
    *,
    rs_nsym: int,
    mask_search_limit: int = 256,
) -> DecodedFrame:
    """Decode a strand, with a bounded fallback search for a damaged mask byte."""
    if not 1 <= mask_search_limit <= 256:
        raise ValueError("mask_search_limit must be between 1 and 256")
    raw = dna_to_bytes(sequence)
    if len(raw) < 1 + _HEADER.size + rs_nsym:
        raise ValueError("Strand is shorter than the OligoArk frame")

    indicated = raw[0]
    candidates = [indicated]
    candidates.extend(
        mask_id
        for mask_id in range(mask_search_limit)
        if mask_id != indicated
    )
    errors: list[Exception] = []
    for mask_id in candidates:
        try:
            return _decode_with_mask(raw, mask_id, rs_nsym)
        except ValueError as exc:
            errors.append(exc)
    raise ValueError("Unable to decode OligoArk strand with any deterministic mask") from errors[-1]

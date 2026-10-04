"""Self-describing, error-corrected strand framing for OligoArk."""

from __future__ import annotations
import struct
import zlib
from dataclasses import dataclass
from .dna import bytes_to_dna, dna_to_bytes, quality_penalty
from .ecc import ECCDecodeError, rs_decode, rs_encode

MAGIC = b"OA"
VERSION = 1
FLAG_PARITY = 0x01
_HEADER = struct.Struct(">2sBBIIHI")
_MASKS = (0x00, 0x55, 0xAA, 0xFF)


@dataclass(frozen=True)
class DecodedFrame:
    index: int
    total_data: int
    payload: bytes
    is_parity: bool


def frame_overhead_bytes(rs_nsym: int) -> int:
    return 1 + _HEADER.size + rs_nsym


def _mask(data: bytes, mask_id: int) -> bytes:
    value = _MASKS[mask_id]
    return bytes(byte ^ value for byte in data)


def encode_frame(payload: bytes, *, index: int, total_data: int, is_parity: bool,
                 rs_nsym: int, adaptive_masks: bool) -> str:
    flags = FLAG_PARITY if is_parity else 0
    header = _HEADER.pack(MAGIC, VERSION, flags, index, total_data, len(payload),
                          zlib.crc32(payload) & 0xFFFFFFFF)
    protected = rs_encode(header + payload, rs_nsym)
    candidates = range(len(_MASKS)) if adaptive_masks else range(1)
    best: tuple[float, str] | None = None
    for mask_id in candidates:
        dna = bytes_to_dna(bytes([mask_id]) + _mask(protected, mask_id))
        score = quality_penalty(dna)
        if best is None or score < best[0]:
            best = (score, dna)
    assert best is not None
    return best[1]


def _decode_with_mask(raw: bytes, mask_id: int, rs_nsym: int) -> DecodedFrame:
    try:
        inner = rs_decode(_mask(raw[1:], mask_id), rs_nsym)
    except ECCDecodeError as exc:
        raise ValueError("Reed-Solomon recovery failed") from exc
    if len(inner) < _HEADER.size:
        raise ValueError("Decoded frame is shorter than the OligoArk header")
    magic, version, flags, index, total, raw_len, crc = _HEADER.unpack(inner[:_HEADER.size])
    if magic != MAGIC or version != VERSION:
        raise ValueError("Not a supported OligoArk strand")
    payload = inner[_HEADER.size:_HEADER.size + raw_len]
    if len(payload) != raw_len:
        raise ValueError("Decoded payload length is inconsistent with frame header")
    if (zlib.crc32(payload) & 0xFFFFFFFF) != crc:
        raise ValueError("Payload CRC check failed")
    return DecodedFrame(index, total, payload, bool(flags & FLAG_PARITY))


def decode_frame(sequence: str, *, rs_nsym: int) -> DecodedFrame:
    raw = dna_to_bytes(sequence)
    if len(raw) < 1 + _HEADER.size + rs_nsym:
        raise ValueError("Strand is shorter than the OligoArk frame")
    indicated = raw[0]
    candidates = [indicated] if indicated < len(_MASKS) else []
    candidates.extend(i for i in range(len(_MASKS)) if i not in candidates)
    errors: list[Exception] = []
    for mask_id in candidates:
        try:
            return _decode_with_mask(raw, mask_id, rs_nsym)
        except ValueError as exc:
            errors.append(exc)
    raise ValueError("Unable to decode OligoArk strand with any payload mask") from errors[-1]

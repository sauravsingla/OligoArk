"""Self-describing, error-corrected strand framing for OligoArk."""

from __future__ import annotations

import binascii
import struct
import zlib
from dataclasses import dataclass

from .dna import (
    SequenceConstraintError,
    SequenceConstraints,
    bytes_to_dna,
    dna_to_bytes,
)
from .ecc import ECCDecodeError, rs_decode, rs_encode

MAGIC = b"OA"
VERSION = 1
COMPACT_VERSION = 2
FLAG_PARITY = 0x01
FLAG_FOUNTAIN = 0x02
_HEADER = struct.Struct(">2sBBIIHI")
_COMPACT_TAG = 0xA0
_COMPACT_TAG_MASK = 0xF0
_COMPACT_FLAG_MASK = 0x0F
_COMPACT_INLINE_TAG = 0xB0
_COMPACT_INLINE_TAG_MASK = 0xF0
_COMPACT_INLINE_MASK_MASK = 0x0F
_LEGACY_MASKS = (0x00, 0x55, 0xAA, 0xFF)


def _inline_selector_for_mask(mask_id: int) -> int:
    """Encode up to 256 masks in one selector byte without changing legacy selectors.

    The original sixteen identifiers use 0xB0..0xBF. The other 240 byte values
    encode masks 16..255. Existing strands remain byte-for-byte decodable.
    """
    if not 0 <= mask_id <= 255:
        raise ValueError("inline mask identifier must fit in one byte")
    if mask_id < 16:
        return _COMPACT_INLINE_TAG | mask_id
    if mask_id < 192:
        return mask_id - 16
    return mask_id


def _inline_mask_from_selector(selector: int) -> int:
    """Inverse of _inline_selector_for_mask, including all legacy 0xB selectors."""
    if selector & _COMPACT_INLINE_TAG_MASK == _COMPACT_INLINE_TAG:
        return selector & _COMPACT_INLINE_MASK_MASK
    if selector < 0xB0:
        return selector + 16
    return selector


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

@dataclass(frozen=True)
class CompactFrameHint:
    """Unverified identity hint recovered from an intact inline compact prefix."""

    index: int
    is_parity: bool
    is_fountain: bool


def compact_inline_frame_hint(
    sequence: str,
    *,
    compact_index_bytes: int = 3,
) -> CompactFrameHint | None:
    """Read only the selector + protected index prefix from a length-shifted strand.

    A single indel occurring after this prefix leaves the prefix byte-aligned and therefore
    gives a cheap identity hint. The hint is never accepted as decoded data: CRC/RS validation
    is still required by decode_frame_resilient before any payload enters recovery.
    """
    sequence = sequence.strip().upper()
    prefix_nt = 4 * (1 + compact_index_bytes)
    if len(sequence) < prefix_nt or set(sequence[:prefix_nt]) - {"A", "C", "G", "T"}:
        return None
    try:
        raw = dna_to_bytes(sequence[:prefix_nt])
    except ValueError:
        return None
    selector = raw[0]
    if selector & _COMPACT_INLINE_TAG_MASK != _COMPACT_INLINE_TAG:
        return None
    mask_id = selector & _COMPACT_INLINE_MASK_MASK
    protected_prefix = _mask(raw[1:], mask_id)
    if len(protected_prefix) != compact_index_bytes:
        return None
    index_bits = 8 * compact_index_bytes - 2
    packed_index = int.from_bytes(protected_prefix, "big")
    flags = packed_index >> index_bits
    if flags & ~(FLAG_PARITY | FLAG_FOUNTAIN):
        return None
    return CompactFrameHint(
        index=packed_index & ((1 << index_bits) - 1),
        is_parity=bool(flags & FLAG_PARITY),
        is_fountain=bool(flags & FLAG_FOUNTAIN),
    )


def frame_overhead_bytes(
    rs_nsym: int,
    *,
    compact_framing: bool = False,
    compact_index_bytes: int = 3,
    inline_mask_framing: bool = False,
) -> int:
    """Return fixed per-strand byte overhead before payload bytes."""
    if compact_framing:
        if not 2 <= compact_index_bytes <= 4:
            raise ValueError("compact_index_bytes must be between 2 and 4")
        if inline_mask_framing:
            # selector/mask byte + protected flags/index + CRC16 + RS parity
            return 1 + compact_index_bytes + 2 + rs_nsym
        # mask byte + control byte + index + CRC16 + Reed-Solomon parity
        return 1 + 1 + compact_index_bytes + 2 + rs_nsym
    if inline_mask_framing:
        raise ValueError("inline_mask_framing requires compact_framing")
    return 1 + _HEADER.size + rs_nsym


def _mask(data: bytes, mask_id: int) -> bytes:
    """Apply a reversible deterministic whitening mask."""
    if not 0 <= mask_id <= 255:
        raise ValueError("mask_id must fit in one byte")
    if mask_id < len(_LEGACY_MASKS):
        value = _LEGACY_MASKS[mask_id]
        if value == 0:
            return data
        return bytes(byte ^ value for byte in data)

    state = (0x9E3779B9 ^ (mask_id * 0x45D9F3B)) & 0xFFFFFFFF
    out = bytearray(len(data))
    for index, byte in enumerate(data):
        state ^= (state << 13) & 0xFFFFFFFF
        state ^= state >> 17
        state ^= (state << 5) & 0xFFFFFFFF
        out[index] = byte ^ (state & 0xFF)
    return bytes(out)


def _flags(*, is_parity: bool, is_fountain: bool) -> int:
    value = 0
    if is_parity:
        value |= FLAG_PARITY
    if is_fountain:
        value |= FLAG_FOUNTAIN
    return value


def _compact_header(payload: bytes, *, index: int, flags: int, index_bytes: int) -> bytes:
    if not 2 <= index_bytes <= 4:
        raise ValueError("compact_index_bytes must be between 2 and 4")
    if not 0 <= index < (1 << (8 * index_bytes)):
        raise ValueError(
            f"compact frame index must fit in {index_bytes} bytes; got index={index}"
        )
    checksum = binascii.crc_hqx(payload, 0xFFFF)
    return (
        bytes([_COMPACT_TAG | (flags & _COMPACT_FLAG_MASK)])
        + index.to_bytes(index_bytes, "big")
        + struct.pack(">H", checksum)
    )


def _compact_inline_header(
    payload: bytes,
    *,
    index: int,
    flags: int,
    index_bytes: int,
) -> bytes:
    """Pack protected frame flags into the high two bits of the compact index."""
    if not 2 <= index_bytes <= 4:
        raise ValueError("compact_index_bytes must be between 2 and 4")
    index_bits = 8 * index_bytes - 2
    if not 0 <= index < (1 << index_bits):
        raise ValueError(
            f"inline compact frame index must fit in {index_bits} bits; got index={index}"
        )
    packed_index = ((flags & 0x03) << index_bits) | index
    checksum = binascii.crc_hqx(payload, 0xFFFF)
    return packed_index.to_bytes(index_bytes, "big") + struct.pack(">H", checksum)


def encode_frame_packed(
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
    compact_framing: bool = False,
    compact_index_bytes: int = 3,
    inline_mask_framing: bool = False,
) -> bytes:
    """Encode one frame directly into compact 2-bit-packed bytes."""
    if not 1 <= mask_search_limit <= 256:
        raise ValueError("mask_search_limit must be between 1 and 256")
    constraints = sequence_constraints or SequenceConstraints()
    constraints.validate()

    flags = _flags(is_parity=is_parity, is_fountain=is_fountain)
    if compact_framing:
        header = (
            _compact_inline_header(
                payload,
                index=index,
                flags=flags,
                index_bytes=compact_index_bytes,
            )
            if inline_mask_framing
            else _compact_header(
                payload,
                index=index,
                flags=flags,
                index_bytes=compact_index_bytes,
            )
        )
    elif inline_mask_framing:
        raise ValueError("inline_mask_framing requires compact_framing")
    else:
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
    # Compact inline framing has a full byte for the selector. Legacy selectors
    # (0xB0..0xBF) remain intact; beyond the first 16, search the remaining 240
    # reversible selectors to avoid rare hard-constraint failures at large scale.
    candidate_ids = range(mask_search_limit) if adaptive_masks else range(1)

    if (
        not adaptive_masks
        and constraints.min_gc_fraction == 0.0
        and constraints.max_gc_fraction == 1.0
        and constraints.max_homopolymer >= 4 * (1 + len(protected))
    ):
        prefix = _COMPACT_INLINE_TAG if inline_mask_framing else 0
        return bytes([prefix]) + _mask(protected, 0)

    # Stop at the first deterministic valid mask. Earlier releases scored every valid
    # candidate and then selected the soft optimum, which multiplied physical-profile
    # encoding cost by the full search budget. Hard GC/homopolymer constraints remain
    # unchanged; this only removes unnecessary work once a valid strand is found.
    for mask_id in candidate_ids:
        prefix = _inline_selector_for_mask(mask_id) if inline_mask_framing else mask_id
        packed = bytes([prefix]) + _mask(protected, mask_id)
        if constraints.accepts(bytes_to_dna(packed)):
            return packed

    raise SequenceConstraintError(
        "No deterministic mask candidate satisfied configured GC/homopolymer constraints"
    )


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
    compact_framing: bool = False,
    compact_index_bytes: int = 3,
    inline_mask_framing: bool = False,
) -> str:
    """Encode one protected strand while enforcing optional hard sequence constraints."""
    packed = encode_frame_packed(
        payload,
        index=index,
        total_data=total_data,
        is_parity=is_parity,
        is_fountain=is_fountain,
        rs_nsym=rs_nsym,
        adaptive_masks=adaptive_masks,
        sequence_constraints=sequence_constraints,
        mask_search_limit=mask_search_limit,
        compact_framing=compact_framing,
        compact_index_bytes=compact_index_bytes,
        inline_mask_framing=inline_mask_framing,
    )
    return bytes_to_dna(packed)


def _decode_legacy_with_mask(raw: bytes, mask_id: int, rs_nsym: int) -> DecodedFrame:
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


def _decode_compact_with_mask(
    raw: bytes,
    mask_id: int,
    rs_nsym: int,
    *,
    expected_total_data: int,
    compact_index_bytes: int,
) -> DecodedFrame:
    protected = _mask(raw[1:], mask_id)
    try:
        inner = rs_decode(protected, rs_nsym)
    except ECCDecodeError as exc:
        raise ValueError("Reed-Solomon recovery failed") from exc

    header_size = 1 + compact_index_bytes + 2
    if len(inner) < header_size:
        raise ValueError("Decoded compact frame is shorter than its header")
    control = inner[0]
    if control & _COMPACT_TAG_MASK != _COMPACT_TAG:
        raise ValueError("Not a supported compact OligoArk strand")
    flags = control & _COMPACT_FLAG_MASK
    if flags & ~(FLAG_PARITY | FLAG_FOUNTAIN):
        raise ValueError("Compact OligoArk frame contains unsupported flags")

    index_start = 1
    index_end = index_start + compact_index_bytes
    index = int.from_bytes(inner[index_start:index_end], "big")
    (checksum,) = struct.unpack(">H", inner[index_end : index_end + 2])
    payload = inner[header_size:]
    if binascii.crc_hqx(payload, 0xFFFF) != checksum:
        raise ValueError("Compact payload CRC16 check failed")

    return DecodedFrame(
        index=index,
        total_data=expected_total_data,
        payload=payload,
        is_parity=bool(flags & FLAG_PARITY),
        is_fountain=bool(flags & FLAG_FOUNTAIN),
    )


def _decode_compact_inline_with_mask(
    raw: bytes,
    mask_id: int,
    rs_nsym: int,
    *,
    expected_total_data: int,
    compact_index_bytes: int,
) -> DecodedFrame:
    # The selector byte is outside RS protection. The complete byte specifies
    # one of 256 masks while retaining the original 0xB0..0xBF mapping. RS and
    # CRC16 still gate any recovered payload; a selector alone proves nothing.
    protected = _mask(raw[1:], mask_id)
    try:
        inner = rs_decode(protected, rs_nsym)
    except ECCDecodeError as exc:
        raise ValueError("Reed-Solomon recovery failed") from exc

    header_size = compact_index_bytes + 2
    if len(inner) < header_size:
        raise ValueError("Decoded inline compact frame is shorter than its header")
    index_bits = 8 * compact_index_bytes - 2
    packed_index = int.from_bytes(inner[:compact_index_bytes], "big")
    flags = packed_index >> index_bits
    index = packed_index & ((1 << index_bits) - 1)
    if flags & ~(FLAG_PARITY | FLAG_FOUNTAIN):
        raise ValueError("Inline compact OligoArk frame contains unsupported flags")
    (checksum,) = struct.unpack(
        ">H",
        inner[compact_index_bytes : compact_index_bytes + 2],
    )
    payload = inner[header_size:]
    if binascii.crc_hqx(payload, 0xFFFF) != checksum:
        raise ValueError("Inline compact payload CRC16 check failed")
    return DecodedFrame(
        index=index,
        total_data=expected_total_data,
        payload=payload,
        is_parity=bool(flags & FLAG_PARITY),
        is_fountain=bool(flags & FLAG_FOUNTAIN),
    )


def _decode_indicated_mask(
    raw: bytes,
    *,
    rs_nsym: int,
    compact_framing: bool,
    compact_index_bytes: int,
    expected_total_data: int | None,
    inline_mask_framing: bool = False,
) -> DecodedFrame:
    """Decode using only the mask identifier carried by this candidate frame."""
    if not raw:
        raise ValueError("empty packed frame")
    mask_id = _inline_mask_from_selector(raw[0]) if inline_mask_framing else raw[0]
    if compact_framing:
        if expected_total_data is None:
            raise ValueError("compact framing requires expected_total_data")
        if inline_mask_framing:
            return _decode_compact_inline_with_mask(
                raw,
                mask_id,
                rs_nsym,
                expected_total_data=expected_total_data,
                compact_index_bytes=compact_index_bytes,
            )
        return _decode_compact_with_mask(
            raw,
            mask_id,
            rs_nsym,
            expected_total_data=expected_total_data,
            compact_index_bytes=compact_index_bytes,
        )
    return _decode_legacy_with_mask(raw, mask_id, rs_nsym)


def decode_frame_packed(
    raw: bytes,
    *,
    rs_nsym: int,
    mask_search_limit: int = 256,
    compact_framing: bool = False,
    compact_index_bytes: int = 3,
    expected_total_data: int | None = None,
    inline_mask_framing: bool = False,
) -> DecodedFrame:
    """Decode a compact 2-bit-packed frame."""
    if not 1 <= mask_search_limit <= 256:
        raise ValueError("mask_search_limit must be between 1 and 256")
    if compact_framing:
        if expected_total_data is None or expected_total_data < 1:
            raise ValueError("compact framing requires expected_total_data")
        minimum = frame_overhead_bytes(
            rs_nsym,
            compact_framing=True,
            compact_index_bytes=compact_index_bytes,
            inline_mask_framing=inline_mask_framing,
        )
    else:
        minimum = 1 + _HEADER.size + rs_nsym

    if len(raw) < minimum:
        raise ValueError("Strand is shorter than the OligoArk frame")

    indicated = _inline_mask_from_selector(raw[0]) if inline_mask_framing else raw[0]
    candidate_limit = mask_search_limit
    candidates = [indicated]
    candidates.extend(mask_id for mask_id in range(candidate_limit) if mask_id != indicated)
    errors: list[Exception] = []
    for mask_id in candidates:
        try:
            if compact_framing:
                if expected_total_data is None:
                    raise ValueError("compact framing requires expected_total_data")
                if inline_mask_framing:
                    return _decode_compact_inline_with_mask(
                        raw,
                        mask_id,
                        rs_nsym,
                        expected_total_data=expected_total_data,
                        compact_index_bytes=compact_index_bytes,
                    )
                return _decode_compact_with_mask(
                    raw,
                    mask_id,
                    rs_nsym,
                    expected_total_data=expected_total_data,
                    compact_index_bytes=compact_index_bytes,
                )
            return _decode_legacy_with_mask(raw, mask_id, rs_nsym)
        except ValueError as exc:
            errors.append(exc)
    raise ValueError("Unable to decode OligoArk strand with any deterministic mask") from errors[-1]


def decode_frame(
    sequence: str,
    *,
    rs_nsym: int,
    mask_search_limit: int = 256,
    compact_framing: bool = False,
    compact_index_bytes: int = 3,
    expected_total_data: int | None = None,
    inline_mask_framing: bool = False,
) -> DecodedFrame:
    """Decode a strand, with a bounded fallback search for a damaged mask byte."""
    return decode_frame_packed(
        dna_to_bytes(sequence),
        rs_nsym=rs_nsym,
        mask_search_limit=mask_search_limit,
        compact_framing=compact_framing,
        compact_index_bytes=compact_index_bytes,
        expected_total_data=expected_total_data,
        inline_mask_framing=inline_mask_framing,
    )


def _single_indel_candidates(sequence: str) -> list[str]:
    """Generate deterministic single-indel repairs that restore byte alignment."""
    sequence = sequence.strip().upper()
    if not sequence or set(sequence) - {"A", "C", "G", "T"}:
        return []
    remainder = len(sequence) % 4
    candidates: list[str] = []
    seen: set[str] = set()

    def add(candidate: str) -> None:
        if candidate not in seen:
            seen.add(candidate)
            candidates.append(candidate)

    if remainder == 1:
        # One inserted base makes an otherwise byte-aligned strand one nucleotide too long.
        for position in range(len(sequence)):
            add(sequence[:position] + sequence[position + 1 :])
    elif remainder == 3:
        # One deleted base makes the observed sequence one nucleotide too short.
        for position in range(len(sequence) + 1):
            for base in "ACGT":
                add(sequence[:position] + base + sequence[position:])
    return candidates


def decode_frame_resilient(
    sequence: str,
    *,
    rs_nsym: int,
    mask_search_limit: int = 256,
    compact_framing: bool = False,
    compact_index_bytes: int = 3,
    expected_total_data: int | None = None,
    max_indel_edits: int = 1,
    inline_mask_framing: bool = False,
) -> DecodedFrame:
    """Decode with a CRC-gated, bounded single-indel realignment fallback.

    Direct decoding remains the fast path. If a single insertion or deletion has shifted the
    2-bit byte boundary, the decoder enumerates only edits that restore byte alignment and
    accepts a candidate only when Reed-Solomon plus the frame checksum validate it. This is
    deliberately bounded: balanced or multi-indel traces still require multi-read consensus.
    """
    try:
        return decode_frame(
            sequence,
            rs_nsym=rs_nsym,
            mask_search_limit=mask_search_limit,
            compact_framing=compact_framing,
            compact_index_bytes=compact_index_bytes,
            expected_total_data=expected_total_data,
            inline_mask_framing=inline_mask_framing,
        )
    except ValueError as direct_error:
        if max_indel_edits < 1:
            raise
        successful: list[DecodedFrame] = []
        identities: set[tuple[int, int, bytes, bool, bool]] = set()
        for candidate in _single_indel_candidates(sequence):
            try:
                raw = dna_to_bytes(candidate)
                frame = _decode_indicated_mask(
                    raw,
                    rs_nsym=rs_nsym,
                    compact_framing=compact_framing,
                    compact_index_bytes=compact_index_bytes,
                    expected_total_data=expected_total_data,
                    inline_mask_framing=inline_mask_framing,
                )
            except ValueError:
                continue
            identity = (
                frame.index,
                frame.total_data,
                frame.payload,
                frame.is_parity,
                frame.is_fountain,
            )
            if identity not in identities:
                identities.add(identity)
                successful.append(frame)
            if len(successful) > 1:
                raise ValueError(
                    "single-indel rescue produced ambiguous valid frames"
                ) from direct_error
        if successful:
            return successful[0]
        raise ValueError(
            "Unable to decode OligoArk strand after single-indel rescue"
        ) from direct_error

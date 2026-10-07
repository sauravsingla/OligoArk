"""Bounded-memory binary archive container for scalable OligoArk storage experiments."""

from __future__ import annotations

import hashlib
import json
import math
import random
import struct
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import BinaryIO, cast

from .archive import ArchiveConfig
from .dna import SequenceConstraints, bytes_to_dna
from .ecc import xor_bytes
from .fountain import indexes_for_seed
from .framing import (
    DecodedFrame,
    decode_frame,
    decode_frame_packed,
    encode_frame_packed,
    frame_overhead_bytes,
)

_MAGIC = b"OAB2"
_HEADER = struct.Struct(">4sI")
_RECORD = struct.Struct(">I")
_HASH_BLOCK = 8 * 1024 * 1024
_QUALITY_SAMPLE = 10_000


@dataclass(frozen=True)
class StreamingArchiveStatistics:
    source_bytes: int
    source_sha256: str
    archive_size_bytes: int
    data_strands: int
    parity_strands: int
    fountain_strands: int
    strand_count: int
    encoded_nucleotides: int
    logical_bits_per_nucleotide: float
    strand_redundancy_ratio: float
    archive_size_overhead_ratio: float
    constraint_pass_rate: float
    constraint_sample_size: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class StreamingFaultProfile:
    dropout_rate: float = 0.0
    substitution_rate: float = 0.0
    insertion_rate: float = 0.0
    deletion_rate: float = 0.0
    seed: int = 7
    controlled_dropout: bool = True

    def validate(self) -> None:
        for name in ("dropout_rate", "substitution_rate", "insertion_rate", "deletion_rate"):
            value = float(getattr(self, name))
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1")


@dataclass(frozen=True)
class StreamingRecoveryReport:
    source_bytes: int
    expected_sha256: str
    output_sha256: str | None
    verified_sha256: bool
    total_data_strands: int
    recovered_data_strands: int
    missing_data_strands: int
    dropped_records: int
    undecodable_records: int
    xor_recovered_strands: int
    fountain_recovered_strands: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _digest(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        while block := handle.read(_HASH_BLOCK):
            size += len(block)
            digest.update(block)
    return size, digest.hexdigest()


def _counts(size: int, config: ArchiveConfig) -> tuple[int, int, int]:
    data = max(1, math.ceil(size / config.chunk_size))
    parity = (
        math.ceil(data / config.parity_group_size)
        if config.redundancy_scheme in {"xor", "hybrid"}
        else 0
    )
    fountain = 0
    if config.redundancy_scheme in {"fountain", "hybrid"} and config.fountain_redundancy > 0:
        fountain = max(1, math.ceil(data * config.fountain_redundancy))
    return data, parity, fountain


def _encoded_nt(size: int, config: ArchiveConfig) -> int:
    data, parity, fountain = _counts(size, config)
    overhead = frame_overhead_bytes(config.rs_nsym)
    if size == 0:
        data_nt = 4 * overhead
    else:
        full, remainder = divmod(size, config.chunk_size)
        data_nt = full * 4 * (overhead + config.chunk_size)
        if remainder:
            data_nt += 4 * (overhead + remainder)
    return data_nt + (parity + fountain) * 4 * (overhead + config.chunk_size)


def _write_record(handle: BinaryIO, packed: bytes) -> None:
    handle.write(_RECORD.pack(len(packed)))
    handle.write(packed)


def _read_record(handle: BinaryIO) -> bytes | None:
    raw = handle.read(_RECORD.size)
    if not raw:
        return None
    if len(raw) != _RECORD.size:
        raise ValueError("truncated streaming record header")
    (length,) = _RECORD.unpack(raw)
    packed = handle.read(length)
    if len(packed) != length:
        raise ValueError("truncated streaming record")
    return packed


def _encode(
    payload: bytes,
    *,
    index: int,
    total: int,
    config: ArchiveConfig,
    parity: bool = False,
    fountain: bool = False,
) -> bytes:
    return encode_frame_packed(
        payload,
        index=index,
        total_data=total,
        is_parity=parity,
        is_fountain=fountain,
        rs_nsym=config.rs_nsym,
        adaptive_masks=config.adaptive_masks,
        sequence_constraints=config.sequence_constraints,
        mask_search_limit=config.mask_search_limit,
    )


def archive_file_streaming(
    source: str | Path,
    destination: str | Path,
    config: ArchiveConfig | None = None,
) -> StreamingArchiveStatistics:
    """Create an oligoark-stream-v2 archive without retaining all strands in RAM."""
    resolved = config or ArchiveConfig()
    resolved.validate()
    source_path, destination_path = Path(source), Path(destination)
    size, sha256 = _digest(source_path)
    data_count, parity_count, fountain_count = _counts(size, resolved)
    nt_count = _encoded_nt(size, resolved)
    metadata: dict[str, object] = {
        "format": "oligoark-stream-v2",
        "original_size": size,
        "sha256": sha256,
        "config": asdict(resolved),
        "data_strands": data_count,
        "parity_strands": parity_count,
        "fountain_strands": fountain_count,
        "strand_count": data_count + parity_count + fountain_count,
        "encoded_nucleotides": nt_count,
        "record_encoding": "length-prefixed 2-bit packed OligoArk v1 frames",
        "note": "software archive container; no wet-lab performance is implied",
    }
    metadata_bytes = json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode()
    reference = SequenceConstraints()
    sampled = passed = 0
    destination_path.parent.mkdir(parents=True, exist_ok=True)

    with destination_path.open("wb") as output, source_path.open("rb") as input_handle:
        output.write(_HEADER.pack(_MAGIC, len(metadata_bytes)))
        output.write(metadata_bytes)
        group: list[bytes] = []
        for index in range(data_count):
            payload = input_handle.read(resolved.chunk_size) if size else b""
            packed = _encode(payload, index=index, total=data_count, config=resolved)
            if sampled < _QUALITY_SAMPLE:
                sampled += 1
                passed += int(reference.accepts(bytes_to_dna(packed)))
            _write_record(output, packed)
            group.append(payload)
            finished = len(group) == resolved.parity_group_size or index == data_count - 1
            if finished and resolved.redundancy_scheme in {"xor", "hybrid"}:
                parity_payload = xor_bytes(group, resolved.chunk_size)
                parity_packed = _encode(
                    parity_payload,
                    index=index // resolved.parity_group_size,
                    total=data_count,
                    config=resolved,
                    parity=True,
                )
                if sampled < _QUALITY_SAMPLE:
                    sampled += 1
                    passed += int(reference.accepts(bytes_to_dna(parity_packed)))
                _write_record(output, parity_packed)
            if finished:
                group.clear()

        for offset in range(fountain_count):
            seed = resolved.fountain_seed + offset
            indexes = indexes_for_seed(data_count, seed, max_degree=resolved.fountain_max_degree)
            parts: list[bytes] = []
            for source_index in indexes:
                input_handle.seek(source_index * resolved.chunk_size)
                parts.append(input_handle.read(resolved.chunk_size))
            packed = _encode(
                xor_bytes(parts, resolved.chunk_size),
                index=seed,
                total=data_count,
                config=resolved,
                fountain=True,
            )
            if sampled < _QUALITY_SAMPLE:
                sampled += 1
                passed += int(reference.accepts(bytes_to_dna(packed)))
            _write_record(output, packed)

    archive_size = destination_path.stat().st_size
    strands = data_count + parity_count + fountain_count
    return StreamingArchiveStatistics(
        source_bytes=size,
        source_sha256=sha256,
        archive_size_bytes=archive_size,
        data_strands=data_count,
        parity_strands=parity_count,
        fountain_strands=fountain_count,
        strand_count=strands,
        encoded_nucleotides=nt_count,
        logical_bits_per_nucleotide=round(size * 8 / nt_count if nt_count else 0.0, 6),
        strand_redundancy_ratio=round((parity_count + fountain_count) / max(1, data_count), 6),
        archive_size_overhead_ratio=round(archive_size / max(1, size), 6),
        constraint_pass_rate=round(passed / sampled if sampled else 0.0, 6),
        constraint_sample_size=sampled,
    )


def _metadata(handle: BinaryIO) -> dict[str, object]:
    raw = handle.read(_HEADER.size)
    if len(raw) != _HEADER.size:
        raise ValueError("truncated streaming archive")
    magic, length = _HEADER.unpack(raw)
    if magic != _MAGIC:
        raise ValueError("not an oligoark-stream-v2 archive")
    value = json.loads(handle.read(length))
    if not isinstance(value, dict) or value.get("format") != "oligoark-stream-v2":
        raise ValueError("invalid streaming metadata")
    return cast(dict[str, object], value)


def _mutate(sequence: str, fault: StreamingFaultProfile, ordinal: int) -> str:
    rng = random.Random((fault.seed << 32) ^ ordinal)
    dna = "ACGT"
    out: list[str] = []
    for base in sequence:
        if rng.random() < fault.deletion_rate:
            continue
        if rng.random() < fault.insertion_rate:
            out.append(rng.choice(dna))
        out.append(
            rng.choice(dna.replace(base, ""))
            if rng.random() < fault.substitution_rate
            else base
        )
    if rng.random() < fault.insertion_rate:
        out.append(rng.choice(dna))
    return "".join(out)


def _drop(
    frame: DecodedFrame,
    config: ArchiveConfig,
    fault: StreamingFaultProfile,
    ordinal: int,
) -> bool:
    if fault.dropout_rate <= 0:
        return False
    if not fault.controlled_dropout:
        return random.Random((fault.seed << 32) ^ ordinal).random() < fault.dropout_rate
    if frame.is_parity or frame.is_fountain:
        return False
    group = max(1, config.parity_group_size)
    if frame.index % group:
        return False
    probability = min(1.0, fault.dropout_rate * group)
    return random.Random((fault.seed << 20) ^ (frame.index // group)).random() < probability


def _decode(
    packed: bytes,
    config: ArchiveConfig,
    fault: StreamingFaultProfile,
    ordinal: int,
) -> tuple[DecodedFrame | None, bool, bool]:
    try:
        clean = decode_frame_packed(
            packed,
            rs_nsym=config.rs_nsym,
            mask_search_limit=config.mask_search_limit,
        )
    except ValueError:
        return None, False, True
    if _drop(clean, config, fault, ordinal):
        return None, True, False
    if fault.substitution_rate == fault.insertion_rate == fault.deletion_rate == 0:
        return clean, False, False
    try:
        return (
            decode_frame(
                _mutate(bytes_to_dna(packed), fault, ordinal),
                rs_nsym=config.rs_nsym,
                mask_search_limit=config.mask_search_limit,
            ),
            False,
            False,
        )
    except ValueError:
        return None, False, True


def _read_chunk(output: BinaryIO, index: int, width: int) -> bytes:
    output.seek(index * width)
    return output.read(width)


def _write_chunk(
    output: BinaryIO,
    known: bytearray,
    index: int,
    payload: bytes,
    width: int,
) -> bool:
    if not 0 <= index < len(known) or known[index]:
        return False
    output.seek(index * width)
    output.write(payload)
    known[index] = 1
    return True


def _recover_fountain_frame(
    frame: DecodedFrame,
    output: BinaryIO,
    known: bytearray,
    total: int,
    config: ArchiveConfig,
) -> bool:
    indexes = indexes_for_seed(
        total,
        frame.index,
        max_degree=config.fountain_max_degree,
    )
    missing_indexes = [index for index in indexes if not known[index]]
    if len(missing_indexes) != 1:
        return False
    parts = [frame.payload] + [
        _read_chunk(output, index, config.chunk_size)
        for index in indexes
        if known[index]
    ]
    payload = xor_bytes(parts, config.chunk_size)
    return _write_chunk(
        output,
        known,
        missing_indexes[0],
        payload,
        config.chunk_size,
    )


def recover_file_streaming(
    archive_path: str | Path,
    output_path: str | Path,
    *,
    fault: StreamingFaultProfile | None = None,
    strict: bool = True,
) -> StreamingRecoveryReport:
    """Recover a compact archive while keeping payload state on disk."""
    resolved_fault = fault or StreamingFaultProfile()
    resolved_fault.validate()
    source, destination = Path(archive_path), Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    dropped = undecodable = xor_recovered = fountain_recovered = 0

    with source.open("rb") as archive_handle:
        metadata = _metadata(archive_handle)
        config_obj = metadata.get("config")
        if not isinstance(config_obj, dict):
            raise ValueError("streaming config metadata must be an object")
        config = ArchiveConfig.from_mapping(cast(dict[str, object], config_obj))
        total = int(cast(int, metadata["data_strands"]))
        records = int(cast(int, metadata["strand_count"]))
        size = int(cast(int, metadata["original_size"]))
        expected = str(metadata["sha256"])
        known = bytearray(total)
        known_count = 0
        fountain_start: int | None = None
        fountain_records = int(cast(int, metadata.get("fountain_strands", 0)))
        fountain_ordinal_start = records - fountain_records

        with destination.open("wb+") as output:
            output.truncate(size)
            for ordinal in range(records):
                record_start = archive_handle.tell()
                if fountain_records and ordinal == fountain_ordinal_start:
                    fountain_start = record_start
                packed = _read_record(archive_handle)
                if packed is None:
                    raise ValueError("archive ended before declared strand count")
                frame, was_dropped, was_bad = _decode(
                    packed,
                    config,
                    resolved_fault,
                    ordinal,
                )
                dropped += int(was_dropped)
                undecodable += int(was_bad)
                if frame is None or frame.total_data != total:
                    continue
                if frame.is_parity:
                    # XOR parity records are emitted directly after their data group. Recover
                    # the group's one controlled erasure immediately instead of retaining all
                    # parity payloads in memory.
                    start = frame.index * config.parity_group_size
                    end = min(start + config.parity_group_size, total)
                    missing_indexes = [
                        index for index in range(start, end) if not known[index]
                    ]
                    if len(missing_indexes) == 1:
                        parts = [frame.payload] + [
                            _read_chunk(output, index, config.chunk_size)
                            for index in range(start, end)
                            if known[index]
                        ]
                        payload = xor_bytes(parts, config.chunk_size)
                        if _write_chunk(
                            output,
                            known,
                            missing_indexes[0],
                            payload,
                            config.chunk_size,
                        ):
                            known_count += 1
                            xor_recovered += 1
                elif frame.is_fountain:
                    # Fountain state stays on disk; opportunistically peel degree-one
                    # equations during the initial pass.
                    if known_count != total and _recover_fountain_frame(
                        frame,
                        output,
                        known,
                        total,
                        config,
                    ):
                        known_count += 1
                        fountain_recovered += 1
                else:
                    if _write_chunk(
                        output,
                        known,
                        frame.index,
                        frame.payload,
                        config.chunk_size,
                    ):
                        known_count += 1

            # If peeling stalled during the initial pass, rescan only the on-disk fountain
            # records. This trades extra sequential I/O for bounded memory instead of retaining
            # every unresolved fountain payload in a Python list.
            changed = True
            while (
                changed
                and fountain_start is not None
                and fountain_records > 0
                and known_count != total
            ):
                changed = False
                archive_handle.seek(fountain_start)
                for offset in range(fountain_records):
                    packed = _read_record(archive_handle)
                    if packed is None:
                        raise ValueError("archive ended inside fountain records")
                    frame, _, _ = _decode(
                        packed,
                        config,
                        resolved_fault,
                        fountain_ordinal_start + offset,
                    )
                    if (
                        frame is not None
                        and frame.is_fountain
                        and frame.total_data == total
                        and _recover_fountain_frame(
                            frame,
                            output,
                            known,
                            total,
                            config,
                        )
                    ):
                        known_count += 1
                        fountain_recovered += 1
                        changed = True
            # Redundancy payloads are padded to chunk width. If the final short chunk was
            # reconstructed, trim any padded tail before hashing the recovered file.
            output.truncate(size)
            output.flush()

    recovered = known_count
    missing = total - recovered
    actual: str | None = None
    if missing == 0:
        _, actual = _digest(destination)
    verified = missing == 0 and actual == expected
    report = StreamingRecoveryReport(
        source_bytes=size,
        expected_sha256=expected,
        output_sha256=actual,
        verified_sha256=verified,
        total_data_strands=total,
        recovered_data_strands=recovered,
        missing_data_strands=missing,
        dropped_records=dropped,
        undecodable_records=undecodable,
        xor_recovered_strands=xor_recovered,
        fountain_recovered_strands=fountain_recovered,
    )
    if strict and not verified:
        raise ValueError(
            "streaming recovery failed SHA-256 verification "
            f"({missing} data strands missing)"
        )
    return report

"""High-level archive, statistics, reconstruction, and recovery pipeline."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import cast

from .dna import sequence_metrics
from .ecc import build_xor_parity, recover_one_missing
from .framing import decode_frame, encode_frame
from .reconstruct import graph_cluster_consensus


@dataclass(frozen=True)
class ArchiveConfig:
    chunk_size: int = 96
    rs_nsym: int = 8
    parity_group_size: int = 8
    adaptive_masks: bool = True

    def validate(self) -> None:
        if self.chunk_size < 8:
            raise ValueError("chunk_size must be at least 8 bytes")
        if self.chunk_size + 18 + self.rs_nsym > 255:
            raise ValueError("chunk_size + protected header + rs_nsym must be <= 255 bytes")
        if not 0 <= self.rs_nsym <= 64:
            raise ValueError("rs_nsym must be between 0 and 64")
        if self.parity_group_size < 2:
            raise ValueError("parity_group_size must be at least 2")

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> ArchiveConfig:
        allowed = {"chunk_size", "rs_nsym", "parity_group_size", "adaptive_masks"}
        unknown = sorted(set(values) - allowed)
        if unknown:
            raise ValueError(f"Unknown archive configuration field(s): {unknown}")

        def integer(name: str, default: int) -> int:
            value = values.get(name, default)
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError(f"{name} must be an integer")
            return value

        adaptive_masks = values.get("adaptive_masks", True)
        if not isinstance(adaptive_masks, bool):
            raise ValueError("adaptive_masks must be a boolean")
        config = cls(
            chunk_size=integer("chunk_size", 96),
            rs_nsym=integer("rs_nsym", 8),
            parity_group_size=integer("parity_group_size", 8),
            adaptive_masks=adaptive_masks,
        )
        config.validate()
        return config


@dataclass(frozen=True)
class ArchiveStatistics:
    strand_count: int
    data_strands: int
    parity_strands: int
    encoded_nucleotides: int
    logical_bits_per_nucleotide: float
    mean_gc_fraction: float
    min_gc_fraction: float
    max_gc_fraction: float
    max_homopolymer: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RecoveryReport:
    direct_attempt_succeeded: bool
    graph_reconstruction_used: bool
    input_reads: int
    consensus_reads: int
    cluster_count: int
    verified_sha256: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass
class DNAArchive:
    metadata: dict[str, object]
    strands: list[str]

    def to_json(self) -> str:
        return json.dumps({"metadata": self.metadata, "strands": self.strands}, indent=2)

    @classmethod
    def from_json(cls, text: str) -> DNAArchive:
        obj = json.loads(text)
        if not isinstance(obj, dict) or not isinstance(obj.get("metadata"), dict):
            raise ValueError("Invalid OligoArk archive JSON")
        strands = obj.get("strands")
        if not isinstance(strands, list) or not all(isinstance(item, str) for item in strands):
            raise ValueError("Archive strands must be a list of DNA strings")
        metadata = cast(dict[str, object], obj["metadata"])
        archive = cls(metadata=metadata, strands=cast(list[str], strands))
        archive.validate()
        return archive

    def validate(self) -> None:
        if self.metadata.get("format") != "oligoark-archive-v1":
            raise ValueError("Unsupported or missing OligoArk archive format")
        required = {"original_size", "sha256", "config", "data_strands", "parity_strands"}
        missing = sorted(required - set(self.metadata))
        if missing:
            raise ValueError(f"Archive metadata is missing required field(s): {missing}")
        config_obj = self.metadata["config"]
        if not isinstance(config_obj, dict):
            raise ValueError("Archive metadata config must be an object")
        ArchiveConfig.from_mapping(config_obj)
        for name in ("original_size", "data_strands", "parity_strands"):
            value = self.metadata[name]
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError(f"{name} metadata must be an integer")
        original_size = cast(int, self.metadata["original_size"])
        data_strands = cast(int, self.metadata["data_strands"])
        parity_strands = cast(int, self.metadata["parity_strands"])
        if original_size < 0:
            raise ValueError("original_size must be non-negative")
        if data_strands <= 0 or parity_strands < 0:
            raise ValueError("data_strands must be positive and parity_strands non-negative")
        if len(self.strands) < data_strands:
            raise ValueError("Archive contains fewer strands than declared data_strands")
        if any(
            not strand or set(strand.upper()) - {"A", "C", "G", "T"}
            for strand in self.strands
        ):
            raise ValueError(
                "Archive strands must be non-empty DNA strings containing only A/C/G/T"
            )
        sha256 = self.metadata["sha256"]
        if not isinstance(sha256, str):
            raise ValueError("sha256 metadata must be a string")
        if len(sha256) != 64 or any(char not in "0123456789abcdef" for char in sha256.lower()):
            raise ValueError("sha256 metadata must be a 64-character hexadecimal digest")

    def save(self, path: str | Path) -> None:
        Path(path).write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> DNAArchive:
        return cls.from_json(Path(path).read_text(encoding="utf-8"))


def archive_statistics(archive: DNAArchive) -> ArchiveStatistics:
    archive.validate()
    metrics = [sequence_metrics(strand) for strand in archive.strands]
    gc_values = [metric.gc_fraction for metric in metrics] or [0.0]
    encoded_nucleotides = sum(len(strand) for strand in archive.strands)
    original_size = int(cast(int, archive.metadata["original_size"]))
    logical_density = (original_size * 8 / encoded_nucleotides) if encoded_nucleotides else 0.0
    return ArchiveStatistics(
        strand_count=len(archive.strands),
        data_strands=int(cast(int, archive.metadata["data_strands"])),
        parity_strands=int(cast(int, archive.metadata["parity_strands"])),
        encoded_nucleotides=encoded_nucleotides,
        logical_bits_per_nucleotide=round(logical_density, 6),
        mean_gc_fraction=round(sum(gc_values) / len(gc_values), 6),
        min_gc_fraction=round(min(gc_values), 6),
        max_gc_fraction=round(max(gc_values), 6),
        max_homopolymer=max((metric.max_homopolymer for metric in metrics), default=0),
    )


def archive_bytes(data: bytes, config: ArchiveConfig | None = None) -> DNAArchive:
    config = config or ArchiveConfig()
    config.validate()
    total = max(1, math.ceil(len(data) / config.chunk_size))
    chunks = [data[i : i + config.chunk_size] for i in range(0, len(data), config.chunk_size)]
    if not chunks:
        chunks = [b""]
    parity = build_xor_parity(chunks, config.parity_group_size, config.chunk_size)
    strands = [
        encode_frame(
            chunk,
            index=index,
            total_data=total,
            is_parity=False,
            rs_nsym=config.rs_nsym,
            adaptive_masks=config.adaptive_masks,
        )
        for index, chunk in enumerate(chunks)
    ]
    strands.extend(
        encode_frame(
            block.payload,
            index=block.group_index,
            total_data=total,
            is_parity=True,
            rs_nsym=config.rs_nsym,
            adaptive_masks=config.adaptive_masks,
        )
        for block in parity
    )
    metadata: dict[str, object] = {
        "format": "oligoark-archive-v1",
        "original_size": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "config": asdict(config),
        "data_strands": total,
        "parity_strands": len(parity),
        "measured": {"encoded_nucleotides": sum(map(len, strands))},
        "note": "DNA strings are software encodings; no wet-lab performance is implied.",
    }
    archive = DNAArchive(metadata, strands)
    stats = archive_statistics(archive)
    archive.metadata["measured"] = stats.to_dict()
    return archive


def _decode_available_frames(
    archive: DNAArchive,
    strands: Iterable[str],
) -> tuple[dict[int, bytes], dict[int, bytes]]:
    config_obj = archive.metadata["config"]
    if not isinstance(config_obj, dict):
        raise ValueError("Archive metadata is missing codec configuration")
    config = ArchiveConfig.from_mapping(config_obj)
    total_data = int(cast(int, archive.metadata["data_strands"]))
    data_chunks: dict[int, bytes] = {}
    parity_chunks: dict[int, bytes] = {}

    for strand in strands:
        try:
            frame = decode_frame(strand, rs_nsym=config.rs_nsym)
        except ValueError:
            continue
        if frame.total_data != total_data:
            continue
        if frame.is_parity:
            parity_chunks.setdefault(frame.index, frame.payload)
        elif 0 <= frame.index < total_data:
            data_chunks.setdefault(frame.index, frame.payload)
    return data_chunks, parity_chunks


def recover_bytes(archive: DNAArchive, strands: Iterable[str] | None = None) -> bytes:
    archive.validate()
    config_obj = archive.metadata["config"]
    if not isinstance(config_obj, dict):
        raise ValueError("Archive metadata is missing codec configuration")
    config = ArchiveConfig.from_mapping(config_obj)
    total_data = int(cast(int, archive.metadata["data_strands"]))
    selected_strands = archive.strands if strands is None else strands
    data_chunks, parity_chunks = _decode_available_frames(archive, selected_strands)
    data_chunks = recover_one_missing(
        data_chunks,
        parity_chunks,
        total_data,
        config.parity_group_size,
        config.chunk_size,
    )
    missing = [index for index in range(total_data) if index not in data_chunks]
    if missing:
        raise ValueError(f"Archive is not recoverable; missing data strand(s): {missing}")
    raw = b"".join(data_chunks[index] for index in range(total_data))
    original_size = int(cast(int, archive.metadata["original_size"]))
    raw = raw[:original_size]
    expected = str(archive.metadata["sha256"])
    actual = hashlib.sha256(raw).hexdigest()
    if actual != expected:
        raise ValueError("Recovered file failed SHA-256 verification")
    return raw


def recover_from_reads(
    archive: DNAArchive,
    reads: Iterable[str],
    *,
    similarity_threshold: float = 0.90,
) -> tuple[bytes, RecoveryReport]:
    """Recover from noisy/duplicate reads, adding graph-consensus reads only when needed."""
    read_list = [read.strip().upper() for read in reads if read.strip()]
    try:
        raw = recover_bytes(archive, read_list)
        return raw, RecoveryReport(True, False, len(read_list), 0, 0, True)
    except ValueError:
        reconstruction = graph_cluster_consensus(read_list, threshold=similarity_threshold)
        augmented_reads = read_list + reconstruction.consensus_reads
        try:
            raw = recover_bytes(archive, augmented_reads)
        except ValueError as reconstructed_error:
            raise ValueError(
                "Recovery failed after direct decode and graph-consensus reconstruction"
            ) from reconstructed_error
        report = RecoveryReport(
            direct_attempt_succeeded=False,
            graph_reconstruction_used=True,
            input_reads=len(read_list),
            consensus_reads=len(reconstruction.consensus_reads),
            cluster_count=len(reconstruction.cluster_sizes),
            verified_sha256=True,
        )
        return raw, report

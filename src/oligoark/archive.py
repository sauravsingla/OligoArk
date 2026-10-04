"""High-level archive, redundancy, statistics, reconstruction, and recovery pipeline."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import cast

from .dna import SequenceConstraints, sequence_metrics
from .ecc import build_xor_parity, recover_one_missing
from .fountain import FountainSymbol, indexes_for_seed, make_symbols, peel_decode
from .framing import decode_frame, encode_frame
from .reconstruct import GraphConsensusReconstructor, ReadReconstructor, TraceConsensusReconstructor

_REDUNDANCY_SCHEMES = {"none", "xor", "fountain", "hybrid"}


@dataclass(frozen=True)
class ArchiveConfig:
    """Configuration that affects generated strands and therefore travels with an archive."""

    chunk_size: int = 96
    rs_nsym: int = 8
    parity_group_size: int = 8
    adaptive_masks: bool = True
    redundancy_scheme: str = "xor"
    fountain_redundancy: float = 0.25
    fountain_seed: int = 1
    fountain_max_degree: int = 4
    min_gc_fraction: float = 0.35
    max_gc_fraction: float = 0.65
    max_homopolymer: int = 4
    mask_search_limit: int = 64

    @property
    def sequence_constraints(self) -> SequenceConstraints:
        return SequenceConstraints(
            min_gc_fraction=self.min_gc_fraction,
            max_gc_fraction=self.max_gc_fraction,
            max_homopolymer=self.max_homopolymer,
        )

    def validate(self) -> None:
        if self.chunk_size < 8:
            raise ValueError("chunk_size must be at least 8 bytes")
        if self.chunk_size + 18 + self.rs_nsym > 255:
            raise ValueError("chunk_size + protected header + rs_nsym must be <= 255 bytes")
        if not 0 <= self.rs_nsym <= 64:
            raise ValueError("rs_nsym must be between 0 and 64")
        if self.parity_group_size < 2:
            raise ValueError("parity_group_size must be at least 2")
        if self.redundancy_scheme not in _REDUNDANCY_SCHEMES:
            raise ValueError(
                f"redundancy_scheme must be one of {sorted(_REDUNDANCY_SCHEMES)}"
            )
        if not 0.0 <= self.fountain_redundancy <= 5.0:
            raise ValueError("fountain_redundancy must be between 0 and 5")
        if not 0 <= self.fountain_seed <= 0xFFFFFFFF:
            raise ValueError("fountain_seed must fit in an unsigned 32-bit integer")
        if not 1 <= self.fountain_max_degree <= 32:
            raise ValueError("fountain_max_degree must be between 1 and 32")
        if not 1 <= self.mask_search_limit <= 256:
            raise ValueError("mask_search_limit must be between 1 and 256")
        self.sequence_constraints.validate()

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> ArchiveConfig:
        """Load current or older archive configuration with strict type validation."""
        allowed = {
            "chunk_size",
            "rs_nsym",
            "parity_group_size",
            "adaptive_masks",
            "redundancy_scheme",
            "fountain_redundancy",
            "fountain_seed",
            "fountain_max_degree",
            "min_gc_fraction",
            "max_gc_fraction",
            "max_homopolymer",
            "mask_search_limit",
        }
        unknown = sorted(set(values) - allowed)
        if unknown:
            raise ValueError(f"Unknown archive configuration field(s): {unknown}")

        def integer(name: str, default: int) -> int:
            value = values.get(name, default)
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError(f"{name} must be an integer")
            return value

        def number(name: str, default: float) -> float:
            value = values.get(name, default)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"{name} must be numeric")
            return float(value)

        adaptive_masks = values.get("adaptive_masks", True)
        if not isinstance(adaptive_masks, bool):
            raise ValueError("adaptive_masks must be a boolean")
        redundancy_scheme = values.get("redundancy_scheme", "xor")
        if not isinstance(redundancy_scheme, str):
            raise ValueError("redundancy_scheme must be a string")

        config = cls(
            chunk_size=integer("chunk_size", 96),
            rs_nsym=integer("rs_nsym", 8),
            parity_group_size=integer("parity_group_size", 8),
            adaptive_masks=adaptive_masks,
            redundancy_scheme=redundancy_scheme,
            fountain_redundancy=number("fountain_redundancy", 0.25),
            fountain_seed=integer("fountain_seed", 1),
            fountain_max_degree=integer("fountain_max_degree", 4),
            min_gc_fraction=number("min_gc_fraction", 0.35),
            max_gc_fraction=number("max_gc_fraction", 0.65),
            max_homopolymer=integer("max_homopolymer", 4),
            mask_search_limit=integer("mask_search_limit", 64),
        )
        config.validate()
        return config


@dataclass(frozen=True)
class ArchiveStatistics:
    strand_count: int
    data_strands: int
    parity_strands: int
    fountain_strands: int
    redundancy_scheme: str
    encoded_nucleotides: int
    logical_bits_per_nucleotide: float
    mean_gc_fraction: float
    min_gc_fraction: float
    max_gc_fraction: float
    max_homopolymer: int
    constraint_pass_rate: float

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
    reconstruction_strategy: str | None = None
    graph_nodes: int = 0
    candidate_pairs: int = 0
    edge_count: int = 0
    component_count: int = 0
    cluster_sizes: tuple[int, ...] = ()
    consensus_lengths: tuple[int, ...] = ()
    reconstruction_runtime_seconds: float = 0.0
    rescue_changed_result: bool = False

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
        fountain_value = self.metadata.get("fountain_strands", 0)
        if isinstance(fountain_value, bool) or not isinstance(fountain_value, int):
            raise ValueError("fountain_strands metadata must be an integer")
        original_size = cast(int, self.metadata["original_size"])
        data_strands = cast(int, self.metadata["data_strands"])
        parity_strands = cast(int, self.metadata["parity_strands"])
        fountain_strands = fountain_value
        if original_size < 0:
            raise ValueError("original_size must be non-negative")
        if data_strands <= 0 or parity_strands < 0 or fountain_strands < 0:
            raise ValueError("strand counts must be non-negative and data_strands positive")
        expected_min = data_strands + parity_strands + fountain_strands
        if len(self.strands) < expected_min:
            raise ValueError("Archive contains fewer strands than declared strand counts")
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
    config_obj = cast(dict[str, object], archive.metadata["config"])
    config = ArchiveConfig.from_mapping(config_obj)
    metrics = [sequence_metrics(strand) for strand in archive.strands]
    gc_values = [metric.gc_fraction for metric in metrics] or [0.0]
    encoded_nucleotides = sum(len(strand) for strand in archive.strands)
    original_size = int(cast(int, archive.metadata["original_size"]))
    logical_density = (original_size * 8 / encoded_nucleotides) if encoded_nucleotides else 0.0
    constraints = config.sequence_constraints
    passes = sum(constraints.accepts(strand) for strand in archive.strands)
    return ArchiveStatistics(
        strand_count=len(archive.strands),
        data_strands=int(cast(int, archive.metadata["data_strands"])),
        parity_strands=int(cast(int, archive.metadata["parity_strands"])),
        fountain_strands=int(cast(int, archive.metadata.get("fountain_strands", 0))),
        redundancy_scheme=config.redundancy_scheme,
        encoded_nucleotides=encoded_nucleotides,
        logical_bits_per_nucleotide=round(logical_density, 6),
        mean_gc_fraction=round(sum(gc_values) / len(gc_values), 6),
        min_gc_fraction=round(min(gc_values), 6),
        max_gc_fraction=round(max(gc_values), 6),
        max_homopolymer=max((metric.max_homopolymer for metric in metrics), default=0),
        constraint_pass_rate=round(passes / max(1, len(archive.strands)), 6),
    )


def _encode_common(
    payload: bytes,
    *,
    index: int,
    total: int,
    config: ArchiveConfig,
    is_parity: bool = False,
    is_fountain: bool = False,
) -> str:
    return encode_frame(
        payload,
        index=index,
        total_data=total,
        is_parity=is_parity,
        is_fountain=is_fountain,
        rs_nsym=config.rs_nsym,
        adaptive_masks=config.adaptive_masks,
        sequence_constraints=config.sequence_constraints,
        mask_search_limit=config.mask_search_limit,
    )


def archive_bytes(data: bytes, config: ArchiveConfig | None = None) -> DNAArchive:
    config = config or ArchiveConfig()
    config.validate()
    total = max(1, math.ceil(len(data) / config.chunk_size))
    chunks = [
        data[index : index + config.chunk_size]
        for index in range(0, len(data), config.chunk_size)
    ]
    if not chunks:
        chunks = [b""]

    parity = []
    if config.redundancy_scheme in {"xor", "hybrid"}:
        parity = build_xor_parity(chunks, config.parity_group_size, config.chunk_size)

    fountain_symbols: list[FountainSymbol] = []
    if config.redundancy_scheme in {"fountain", "hybrid"} and config.fountain_redundancy > 0:
        count = max(1, math.ceil(total * config.fountain_redundancy))
        if config.fountain_seed + count - 1 > 0xFFFFFFFF:
            raise ValueError("fountain seed range exceeds unsigned 32-bit frame index")
        fountain_symbols = make_symbols(
            chunks,
            count=count,
            width=config.chunk_size,
            seed=config.fountain_seed,
            max_degree=config.fountain_max_degree,
        )

    strands = [
        _encode_common(chunk, index=index, total=total, config=config)
        for index, chunk in enumerate(chunks)
    ]
    strands.extend(
        _encode_common(
            block.payload,
            index=block.group_index,
            total=total,
            config=config,
            is_parity=True,
        )
        for block in parity
    )
    strands.extend(
        _encode_common(
            symbol.payload,
            index=symbol.seed,
            total=total,
            config=config,
            is_fountain=True,
        )
        for symbol in fountain_symbols
    )

    metadata: dict[str, object] = {
        "format": "oligoark-archive-v1",
        "original_size": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "config": asdict(config),
        "data_strands": total,
        "parity_strands": len(parity),
        "fountain_strands": len(fountain_symbols),
        "measured": {"encoded_nucleotides": sum(map(len, strands))},
        "note": "DNA strings are software encodings; no wet-lab performance is implied.",
    }
    archive = DNAArchive(metadata, strands)
    archive.metadata["measured"] = archive_statistics(archive).to_dict()
    return archive


def _decode_available_frames(
    archive: DNAArchive,
    strands: Iterable[str],
) -> tuple[dict[int, bytes], dict[int, bytes], list[FountainSymbol]]:
    config_obj = cast(dict[str, object], archive.metadata["config"])
    config = ArchiveConfig.from_mapping(config_obj)
    total_data = int(cast(int, archive.metadata["data_strands"]))
    data_chunks: dict[int, bytes] = {}
    parity_chunks: dict[int, bytes] = {}
    fountain_symbols: list[FountainSymbol] = []

    for strand in strands:
        try:
            frame = decode_frame(
                strand,
                rs_nsym=config.rs_nsym,
                mask_search_limit=config.mask_search_limit,
            )
        except ValueError:
            continue
        if frame.total_data != total_data:
            continue
        if frame.is_fountain:
            indexes = indexes_for_seed(
                total_data,
                frame.index,
                max_degree=config.fountain_max_degree,
            )
            fountain_symbols.append(FountainSymbol(frame.index, indexes, frame.payload))
        elif frame.is_parity:
            parity_chunks.setdefault(frame.index, frame.payload)
        elif 0 <= frame.index < total_data:
            data_chunks.setdefault(frame.index, frame.payload)
    return data_chunks, parity_chunks, fountain_symbols


def _apply_redundancy(
    data_chunks: dict[int, bytes],
    parity_chunks: dict[int, bytes],
    fountain_symbols: list[FountainSymbol],
    total_data: int,
    config: ArchiveConfig,
) -> dict[int, bytes]:
    recovered = dict(data_chunks)
    for _ in range(3):
        before = len(recovered)
        if parity_chunks:
            recovered = recover_one_missing(
                recovered,
                parity_chunks,
                total_data,
                config.parity_group_size,
                config.chunk_size,
            )
        if fountain_symbols:
            recovered = peel_decode(
                recovered,
                fountain_symbols,
                total=total_data,
                width=config.chunk_size,
            )
        if len(recovered) == before:
            break
    return recovered


def recover_bytes(archive: DNAArchive, strands: Iterable[str] | None = None) -> bytes:
    archive.validate()
    config_obj = cast(dict[str, object], archive.metadata["config"])
    config = ArchiveConfig.from_mapping(config_obj)
    total_data = int(cast(int, archive.metadata["data_strands"]))
    selected_strands = archive.strands if strands is None else strands
    data_chunks, parity_chunks, fountain_symbols = _decode_available_frames(
        archive,
        selected_strands,
    )
    data_chunks = _apply_redundancy(
        data_chunks,
        parity_chunks,
        fountain_symbols,
        total_data,
        config,
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
    reconstructor: ReadReconstructor | None = None,
) -> tuple[bytes, RecoveryReport]:
    """Recover noisy/duplicate reads with a pluggable reconstruction fallback."""
    archive.validate()
    read_list = [read.strip().upper() for read in reads if read.strip()]
    try:
        raw = recover_bytes(archive, read_list)
        return raw, RecoveryReport(
            direct_attempt_succeeded=True,
            graph_reconstruction_used=False,
            input_reads=len(read_list),
            consensus_reads=0,
            cluster_count=0,
            verified_sha256=True,
            reconstruction_strategy=None,
            rescue_changed_result=False,
        )
    except ValueError:
        resolved = reconstructor or GraphConsensusReconstructor(
            threshold=similarity_threshold
        )
        reconstruction = resolved.reconstruct(read_list)
        augmented_reads = read_list + reconstruction.consensus_reads
        try:
            raw = recover_bytes(archive, augmented_reads)
        except ValueError as reconstructed_error:
            raise ValueError(
                "Recovery failed after direct decode and reconstruction fallback"
            ) from reconstructed_error

        strategy = type(resolved).__name__
        report = RecoveryReport(
            direct_attempt_succeeded=False,
            graph_reconstruction_used=isinstance(
                resolved,
                (GraphConsensusReconstructor, TraceConsensusReconstructor),
            ),
            input_reads=len(read_list),
            consensus_reads=len(reconstruction.consensus_reads),
            cluster_count=len(reconstruction.cluster_sizes),
            verified_sha256=True,
            reconstruction_strategy=strategy,
            graph_nodes=reconstruction.node_count,
            candidate_pairs=reconstruction.candidate_pairs,
            edge_count=reconstruction.edge_count,
            component_count=reconstruction.component_count,
            cluster_sizes=tuple(reconstruction.cluster_sizes),
            consensus_lengths=reconstruction.consensus_lengths,
            reconstruction_runtime_seconds=reconstruction.runtime_seconds,
            rescue_changed_result=True,
        )
        return raw, report

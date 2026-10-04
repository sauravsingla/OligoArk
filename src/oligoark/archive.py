"""High-level archive and recovery pipeline."""

from __future__ import annotations
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable
from .ecc import build_xor_parity, recover_one_missing
from .framing import decode_frame, encode_frame


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


@dataclass
class DNAArchive:
    metadata: dict[str, object]
    strands: list[str]

    def to_json(self) -> str:
        return json.dumps({"metadata": self.metadata, "strands": self.strands}, indent=2)

    @classmethod
    def from_json(cls, text: str) -> "DNAArchive":
        obj = json.loads(text)
        if not isinstance(obj, dict) or not isinstance(obj.get("metadata"), dict):
            raise ValueError("Invalid OligoArk archive JSON")
        strands = obj.get("strands")
        if not isinstance(strands, list) or not all(isinstance(x, str) for x in strands):
            raise ValueError("Archive strands must be a list of DNA strings")
        return cls(metadata=obj["metadata"], strands=strands)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> "DNAArchive":
        return cls.from_json(Path(path).read_text(encoding="utf-8"))


def archive_bytes(data: bytes, config: ArchiveConfig | None = None) -> DNAArchive:
    config = config or ArchiveConfig()
    config.validate()
    total = max(1, math.ceil(len(data) / config.chunk_size))
    chunks = [data[i:i+config.chunk_size] for i in range(0, len(data), config.chunk_size)]
    if not chunks:
        chunks = [b""]
    parity = build_xor_parity(chunks, config.parity_group_size, config.chunk_size)
    strands = [
        encode_frame(chunk, index=index, total_data=total, is_parity=False,
                     rs_nsym=config.rs_nsym, adaptive_masks=config.adaptive_masks)
        for index, chunk in enumerate(chunks)
    ]
    strands.extend(
        encode_frame(block.payload, index=block.group_index, total_data=total, is_parity=True,
                     rs_nsym=config.rs_nsym, adaptive_masks=config.adaptive_masks)
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
    return DNAArchive(metadata, strands)


def recover_bytes(archive: DNAArchive, strands: Iterable[str] | None = None) -> bytes:
    config_obj = archive.metadata.get("config")
    if not isinstance(config_obj, dict):
        raise ValueError("Archive metadata is missing codec configuration")
    config = ArchiveConfig(**config_obj)
    data_chunks: dict[int, bytes] = {}
    parity_chunks: dict[int, bytes] = {}
    total_data = int(archive.metadata["data_strands"])
    for strand in strands if strands is not None else archive.strands:
        try:
            frame = decode_frame(strand, rs_nsym=config.rs_nsym)
        except ValueError:
            continue
        if frame.total_data != total_data:
            continue
        if frame.is_parity:
            parity_chunks.setdefault(frame.index, frame.payload)
        else:
            data_chunks.setdefault(frame.index, frame.payload)
    data_chunks = recover_one_missing(
        data_chunks, parity_chunks, total_data, config.parity_group_size, config.chunk_size)
    missing = [idx for idx in range(total_data) if idx not in data_chunks]
    if missing:
        raise ValueError(f"Archive is not recoverable; missing data strand(s): {missing}")
    raw = b"".join(data_chunks[idx] for idx in range(total_data))
    raw = raw[:int(archive.metadata["original_size"])]
    expected = str(archive.metadata["sha256"])
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError("Recovered file failed SHA-256 verification")
    return raw

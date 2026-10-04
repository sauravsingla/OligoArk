"""Runtime configuration helpers for OligoArk applications."""

from __future__ import annotations

import json
import os
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class RuntimeConfig:
    """Operational settings that do not change the archive wire format."""

    log_level: str = "INFO"
    reconstruction_threshold: float = 0.90
    max_api_payload_bytes: int = 10 * 1024 * 1024

    def validate(self) -> None:
        if self.log_level.upper() not in {"CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"}:
            raise ValueError("log_level must be a standard Python logging level")
        if not 0.0 <= self.reconstruction_threshold <= 1.0:
            raise ValueError("reconstruction_threshold must be between 0 and 1")
        if self.max_api_payload_bytes <= 0:
            raise ValueError("max_api_payload_bytes must be positive")

    def to_dict(self) -> dict[str, object]:
        return asdict(self)

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> RuntimeConfig:
        allowed = {"log_level", "reconstruction_threshold", "max_api_payload_bytes"}
        unknown = sorted(set(values) - allowed)
        if unknown:
            raise ValueError(f"Unknown runtime configuration field(s): {unknown}")

        log_level = values.get("log_level", "INFO")
        threshold = values.get("reconstruction_threshold", 0.90)
        max_payload = values.get("max_api_payload_bytes", 10 * 1024 * 1024)
        if not isinstance(log_level, str):
            raise ValueError("log_level must be a string")
        if isinstance(threshold, bool) or not isinstance(threshold, (int, float)):
            raise ValueError("reconstruction_threshold must be numeric")
        if isinstance(max_payload, bool) or not isinstance(max_payload, int):
            raise ValueError("max_api_payload_bytes must be an integer")

        config = cls(log_level, float(threshold), max_payload)
        config.validate()
        return config

    @classmethod
    def from_json_file(cls, path: str | Path) -> RuntimeConfig:
        parsed = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(parsed, dict):
            raise ValueError("Runtime configuration JSON must contain an object")
        return cls.from_mapping(parsed)

    @classmethod
    def from_environment(cls, prefix: str = "OLIGOARK_") -> RuntimeConfig:
        values: dict[str, object] = {}
        if value := os.getenv(f"{prefix}LOG_LEVEL"):
            values["log_level"] = value
        if value := os.getenv(f"{prefix}RECONSTRUCTION_THRESHOLD"):
            values["reconstruction_threshold"] = float(value)
        if value := os.getenv(f"{prefix}MAX_API_PAYLOAD_BYTES"):
            values["max_api_payload_bytes"] = int(value)
        return cls.from_mapping(values)

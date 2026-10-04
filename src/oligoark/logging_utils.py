"""Logging configuration for CLI and service entry points."""

from __future__ import annotations

import logging


def configure_logging(level: str = "INFO") -> None:
    """Configure concise structured-enough console logging without external dependencies."""
    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Unknown logging level: {level}")
    logging.basicConfig(
        level=numeric_level,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        force=True,
    )

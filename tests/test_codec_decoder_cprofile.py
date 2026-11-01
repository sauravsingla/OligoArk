"""Opt-in matched-codec decoder profiling preserves outcomes and partial failures."""

from __future__ import annotations

import pstats
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from run_matched_codec_scale import _profiled_decode  # noqa: E402


def test_profiled_decode_disabled(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.delenv("OLIGOARK_CODEC_CPROFILE_DIR", raising=False)
    assert _profiled_decode(lambda reads: reads[::-1], "ACGT", "method", "indel-low", 0) == "TGCA"
    assert not list(tmp_path.iterdir())


def test_profiled_decode_retains_success_and_failure(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setenv("OLIGOARK_CODEC_CPROFILE_DIR", str(tmp_path))
    assert _profiled_decode(lambda reads: len(reads), "ACGT", "method", "indel-low", 0) == 4

    def fail(_reads: str) -> None:
        raise TimeoutError("synthetic deadline")

    with pytest.raises(TimeoutError, match="synthetic deadline"):
        _profiled_decode(fail, "ACGT", "method", "indel-low", 1)

    for trial in (0, 1):
        path = tmp_path / f"method-indel-low-trial-{trial}.prof"
        assert path.is_file()
        stats = pstats.Stats(str(path))
        assert stats.total_calls > 0
        assert stats.stats

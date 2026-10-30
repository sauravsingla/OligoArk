"""Codec timeouts retain a bounded interrupted stack without changing trial outcomes."""

from __future__ import annotations

import inspect
import json
import signal
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

from run_matched_codec_scale import _TrialDeadline, _alarm_handler, _condition_row

pytestmark = pytest.mark.skipif(
    not hasattr(signal, "SIGALRM") or not hasattr(signal, "setitimer"),
    reason="SIGALRM trial diagnostics require POSIX timers",
)


def _run_trial(decode: object) -> dict[str, object]:
    return _condition_row(
        method="oligoark-compact-hybrid",
        payload=b"x",
        sequences=("ACGT" * 38,),
        decode=decode,
        condition_name="clean",
        trials=1,
        redundancy=0.25,
        target_max_strand_nt=152,
        data_units=1,
        encoded_nucleotides=152,
        logical_bits_per_nucleotide=0.05,
        encode_seconds=0.01,
        condition_timeout_seconds=1,
    )


def test_timeout_retains_interrupted_stack_and_failed_sha256() -> None:
    def slow_decode(_reads: object) -> bytes:
        _alarm_handler(signal.SIGALRM, inspect.currentframe())
        raise AssertionError("unreachable after the deadline")

    row = _run_trial(slow_decode)
    assert row["timeout_trials"] == 1
    assert row["failed_trials"] == 1
    assert row["successes"] == 0
    trial = row["trial_results"][0]
    assert trial["timed_out"] is True
    assert trial["timeout_stage"] == "decode"
    assert trial["sha256_verified"] is False
    stack = trial["timeout_stack"]
    assert stack
    assert stack[-1]["function"] == "slow_decode"
    assert all(not frame["file"].startswith("/") for frame in stack)
    assert all(set(frame) == {"file", "line", "function"} for frame in stack)
    json.dumps(row)


def test_successful_trial_does_not_invent_timeout_evidence() -> None:
    row = _run_trial(lambda _reads: b"x")
    assert row["timeout_trials"] == 0
    assert row["successes"] == 1
    trial = row["trial_results"][0]
    assert trial["sha256_verified"] is True
    assert "timeout_stack" not in trial


def test_deadline_with_missing_frame_retains_failure() -> None:
    with pytest.raises(_TrialDeadline) as exc:
        _alarm_handler(signal.SIGALRM, None)
    assert exc.value.frames == []

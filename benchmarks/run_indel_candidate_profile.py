"""Controlled reference-vs-streamed single-indel candidate memory profile.

This microbenchmark measures enumeration only. It does not measure end-to-end
indel rescue, codec throughput, 10 MiB timeout recovery or physical DNA.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import random
import statistics
import sys
import time
import tracemalloc
from collections.abc import Callable, Iterable
from pathlib import Path

from oligoark.framing import _single_indel_candidates


def _reference_candidates(sequence: str) -> list[str]:
    """Original list-and-set implementation retained as a comparison baseline."""
    sequence = sequence.strip().upper()
    if not sequence or set(sequence) - {"A", "C", "G", "T"}:
        return []
    candidates: list[str] = []
    seen: set[str] = set()

    def add(candidate: str) -> None:
        if candidate not in seen:
            seen.add(candidate)
            candidates.append(candidate)

    if len(sequence) % 4 == 1:
        for position in range(len(sequence)):
            add(sequence[:position] + sequence[position + 1 :])
    elif len(sequence) % 4 == 3:
        for position in range(len(sequence) + 1):
            for base in "ACGT":
                add(sequence[:position] + base + sequence[position:])
    return candidates


def _measure(
    candidate_fn: Callable[[str], Iterable[str]],
    sequence: str,
    repeats: int,
) -> dict[str, object]:
    times_ns: list[int] = []
    peak_bytes: list[int] = []
    fingerprints: list[str] = []
    counts: list[int] = []
    for _ in range(repeats):
        tracemalloc.start()
        started = time.perf_counter_ns()
        digest = hashlib.sha256()
        count = 0
        for candidate in candidate_fn(sequence):
            digest.update(candidate.encode("ascii"))
            digest.update(bytes([0]))
            count += 1
        times_ns.append(time.perf_counter_ns() - started)
        peak_bytes.append(tracemalloc.get_traced_memory()[1])
        tracemalloc.stop()
        fingerprints.append(digest.hexdigest())
        counts.append(count)
    if len(set(fingerprints)) != 1 or len(set(counts)) != 1:
        raise RuntimeError("nondeterministic candidate enumeration")
    return {
        "count": counts[0],
        "sha256": fingerprints[0],
        "timings_ns": times_ns,
        "peak_traced_bytes": peak_bytes,
        "median_time_ns": statistics.median(times_ns),
        "median_peak_traced_bytes": statistics.median(peak_bytes),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument(
        "--output", type=Path, default=Path("performance-results/indel-candidates.json")
    )
    args = parser.parse_args()
    if args.repeats < 3:
        parser.error("--repeats must be at least 3")

    rng = random.Random(args.seed)
    cases: list[dict[str, object]] = []
    for length in (599, 601):
        sequence = "".join(rng.choices("ACGT", k=length))
        reference = _measure(_reference_candidates, sequence, args.repeats)
        streamed = _measure(_single_indel_candidates, sequence, args.repeats)
        if (reference["count"], reference["sha256"]) != (
            streamed["count"], streamed["sha256"]
        ):
            raise SystemExit(f"candidate count/order mismatch at length {length}")
        reference_peak = float(reference["median_peak_traced_bytes"])
        streamed_peak = float(streamed["median_peak_traced_bytes"])
        if streamed_peak > reference_peak * 0.5:
            raise SystemExit(
                f"streamed candidate memory regression at length {length}: "
                f"{streamed_peak} vs {reference_peak} traced bytes"
            )
        cases.append(
            {
                "length_nt": length,
                "input_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "reference": reference,
                "streamed": streamed,
                "exact_candidate_order_parity": True,
                "reference_over_streamed_median_peak_traced_ratio": round(
                    reference_peak / max(streamed_peak, 1.0), 6
                ),
                "reference_over_streamed_median_time_ratio": round(
                    float(reference["median_time_ns"])
                    / max(float(streamed["median_time_ns"]), 1.0),
                    6,
                ),
            }
        )
    report = {
        "seed": args.seed,
        "repeats": args.repeats,
        "python": sys.version,
        "platform": platform.platform(),
        "scope": "candidate enumeration only; not end-to-end codec performance",
        "memory_gate": "streamed median tracemalloc peak <= 50% of reference",
        "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + chr(10))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

"""Matched DNA-to-byte decoding microbenchmark (not end-to-end codec speed).

Run: python benchmarks/run_dna_decode_profile.py --seed 20261031 --repeats 7
Reports per-call wall time and traced Python allocation peaks. The memory gate
compares both implementations in the same interpreter; timing is report-only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import statistics
import sys
import time
import tracemalloc
from collections.abc import Callable
from pathlib import Path

from oligoark.dna import BASE_TO_BITS, dna_to_bytes


def reference(sequence: str) -> bytes:
    """Original list-based implementation."""
    sequence = sequence.strip().upper()
    if len(sequence) % 4:
        raise ValueError("DNA sequence length must be divisible by 4")
    try:
        vals = [BASE_TO_BITS[base] for base in sequence]
    except KeyError as exc:
        raise ValueError(f"Invalid DNA base: {exc.args[0]!r}") from exc
    return bytes(
        (vals[i] << 6) | (vals[i + 1] << 4) | (vals[i + 2] << 2) | vals[i + 3]
        for i in range(0, len(vals), 4)
    )


def measure(
    fn: Callable[[str], bytes], sequence: str, *, repeats: int, loops: int
) -> dict[str, object]:
    timing_ns: list[float] = []
    peaks: list[int] = []
    fingerprint = hashlib.sha256(fn(sequence)).hexdigest()
    for _ in range(repeats):
        started = time.perf_counter_ns()
        for _ in range(loops):
            if hashlib.sha256(fn(sequence)).hexdigest() != fingerprint:
                raise RuntimeError("nondeterministic DNA decoding")
        timing_ns.append((time.perf_counter_ns() - started) / loops)
        tracemalloc.start()
        output = fn(sequence)
        peaks.append(tracemalloc.get_traced_memory()[1])
        tracemalloc.stop()
        if hashlib.sha256(output).hexdigest() != fingerprint:
            raise RuntimeError("traced DNA decoding differs")
    return {
        "output_sha256": fingerprint,
        "per_call_timing_ns": timing_ns,
        "traced_peak_bytes": peaks,
        "median_per_call_ns": statistics.median(timing_ns),
        "median_traced_peak_bytes": statistics.median(peaks),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument(
        "--output", type=Path, default=Path("performance-results/dna-decode.json")
    )
    args = parser.parse_args()
    if args.repeats < 3:
        parser.error("--repeats must be at least 3")

    rng = random.Random(args.seed)
    cases: list[dict[str, object]] = []
    for length, loops in ((152, 1000), (600, 500), (1024, 250)):
        sequence = "".join(rng.choices("ACGT", k=length))
        before = measure(reference, sequence, repeats=args.repeats, loops=loops)
        after = measure(dna_to_bytes, sequence, repeats=args.repeats, loops=loops)
        if before["output_sha256"] != after["output_sha256"]:
            raise SystemExit(f"exact DNA decode parity failure at {length} nt")
        before_peak = float(before["median_traced_peak_bytes"])
        after_peak = float(after["median_traced_peak_bytes"])
        if after_peak > 0.6 * before_peak:
            raise SystemExit(
                f"DNA decode memory regression at {length} nt: "
                f"{after_peak} vs {before_peak} traced bytes"
            )
        cases.append(
            {
                "length_nt": length,
                "loops_per_repeat": loops,
                "input_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "reference": before,
                "optimized": after,
                "exact_output_parity": True,
                "reference_over_optimized_traced_peak_ratio": round(
                    before_peak / max(after_peak, 1.0), 6
                ),
                "reference_over_optimized_median_time_ratio": round(
                    float(before["median_per_call_ns"])
                    / max(float(after["median_per_call_ns"]), 1.0),
                    6,
                ),
            }
        )
    report = {
        "seed": args.seed,
        "repeats": args.repeats,
        "python": sys.version,
        "platform": platform.platform(),
        "github_sha": os.environ.get("GITHUB_SHA", "local"),
        "scope": "DNA-to-byte conversion only, not end-to-end codec decoding",
        "memory_gate": "optimized median traced peak <= 60% of reference",
        "timing_gate": "none; shared CI wall-clock measurements are descriptive",
        "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

"""Matched randomized reference/lookup RS syndrome microbenchmark.

Benchmark arithmetic only, not full DNA reconstruction or archive recovery.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import random
import resource
import statistics
import sys
import time
from pathlib import Path

from oligoark.ecc import _gf_pow, _poly_eval, _syndromes


def reference_syndromes(msg: bytes, nsym: int) -> list[int]:
    """Pre-change syndrome calculation, used as a parity oracle."""
    return [0] + [_poly_eval(msg, _gf_pow(2, i)) for i in range(nsym)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261102)
    parser.add_argument("--cases", type=int, default=400)
    parser.add_argument("--message-length", type=int, default=96)
    parser.add_argument("--nsym", type=int, default=12)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("performance-results/rs-syndrome-lookup.json"),
    )
    args = parser.parse_args()
    if not (1 <= args.message_length <= 255):
        parser.error("message-length must be between 1 and 255")
    if not (1 <= args.nsym < 255):
        parser.error("nsym must be between 1 and 254")
    if args.cases < 1 or args.repeats < 1:
        parser.error("cases and repeats must be positive")

    rng = random.Random(args.seed)
    messages = [rng.randbytes(args.message_length) for _ in range(args.cases)]
    expected = [reference_syndromes(item, args.nsym) for item in messages]
    observed = [_syndromes(item, args.nsym) for item in messages]
    if observed != expected:
        raise RuntimeError("Syndrome lookup differs from reference")
    fingerprint = hashlib.sha256(repr(expected).encode()).hexdigest()

    functions = {"reference": reference_syndromes, "lookup": _syndromes}
    timings: dict[str, list[float]] = {"reference": [], "lookup": []}
    cpu_timings: dict[str, list[float]] = {"reference": [], "lookup": []}
    for fn in functions.values():
        for item in messages:
            fn(item, args.nsym)

    for repeat in range(args.repeats):
        order = ("reference", "lookup") if repeat % 2 == 0 else ("lookup", "reference")
        for name in order:
            fn = functions[name]
            wall_start = time.perf_counter()
            cpu_start = time.process_time()
            for item in messages:
                fn(item, args.nsym)
            cpu_timings[name].append(time.process_time() - cpu_start)
            timings[name].append(time.perf_counter() - wall_start)

    baseline = statistics.median(timings["reference"])
    optimized = statistics.median(timings["lookup"])
    raw_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_mib = raw_rss / (1024 * 1024 if sys.platform == "darwin" else 1024)
    result = {
        "parameters": {
            "seed": args.seed,
            "cases": args.cases,
            "message_length": args.message_length,
            "nsym": args.nsym,
            "repeats": args.repeats,
            "alternating_order": True,
        },
        "environment": {"python": sys.version, "platform": platform.platform()},
        "exact_syndrome_parity": True,
        "results_sha256": fingerprint,
        "wall_seconds": timings,
        "cpu_seconds": cpu_timings,
        "median_wall_seconds": {"reference": baseline, "lookup": optimized},
        "median_cpu_seconds": {
            name: statistics.median(samples) for name, samples in cpu_timings.items()
        },
        "reference_over_lookup_wall_ratio": baseline / optimized,
        "peak_process_rss_mib": round(rss_mib, 3),
        "limitations": (
            "Isolated warm-cache Python RS syndrome arithmetic; process RSS is shared "
            "across methods. Does not prove improved noisy-channel exact recovery, "
            "end-to-end throughput or indel timeout resolution."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

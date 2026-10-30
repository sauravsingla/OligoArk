"""Matched, deterministic Chien-search reference vs optimized microbenchmark.

Measures only error-location evaluation; not end-to-end archive throughput.
Both implementations run on the same interpreter and randomized valid locators.
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
from pathlib import Path

from oligoark.ecc import _find_errors, _gf_pow, _poly_eval, _poly_mul


def reference_find_errors(locator: list[int], message_len: int) -> list[int]:
    """Generic evaluator from the pre-optimization implementation."""
    errors = len(locator) - 1
    positions: list[int] = []
    for i in range(message_len):
        if _poly_eval(locator, _gf_pow(2, i)) == 0:
            positions.append(message_len - 1 - i)
    if len(positions) != errors:
        raise ValueError("could not locate all Reed-Solomon symbol errors")
    return positions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--cases-per-degree", type=int, default=40)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output", type=Path, default=Path("performance-results/rs-locator.json"))
    args = parser.parse_args()
    if args.cases_per_degree < 1 or args.repeats < 1:
        parser.error("cases-per-degree and repeats must be positive")

    rng = random.Random(args.seed)
    cases: list[list[int]] = []
    for degree in (1, 2, 4, 8, 16):
        for _ in range(args.cases_per_degree):
            locator = [1]
            for root in rng.sample(range(255), degree):
                locator = _poly_mul(locator, [1, _gf_pow(2, root)])
            cases.append(locator)

    expected = [reference_find_errors(locator, 255) for locator in cases]
    observed = [_find_errors(locator, 255) for locator in cases]
    if expected != observed:
        raise RuntimeError("optimized error positions differ from reference")
    fingerprint = hashlib.sha256(repr(expected).encode()).hexdigest()
    timings: dict[str, list[float]] = {}
    for name, fn in (("reference", reference_find_errors), ("optimized", _find_errors)):
        for _ in range(2):
            for locator in cases:
                fn(locator, 255)
        samples = []
        for _ in range(args.repeats):
            start = time.perf_counter()
            for locator in cases:
                fn(locator, 255)
            samples.append(time.perf_counter() - start)
        timings[name] = samples

    before = statistics.median(timings["reference"])
    after = statistics.median(timings["optimized"])
    payload = {
        "parameters": {
            "seed": args.seed,
            "cases_per_degree": args.cases_per_degree,
            "repeats": args.repeats,
            "message_length": 255,
            "degrees": [1, 2, 4, 8, 16],
        },
        "environment": {"python": sys.version, "platform": platform.platform()},
        "exact_position_parity": True,
        "position_results_sha256": fingerprint,
        "timings_seconds": timings,
        "median_seconds": {"reference": before, "optimized": after},
        "median_reference_over_optimized": before / after,
        "limitations": (
            "Single-process matched microbenchmark; no end-to-end codec speedup, "
            "statistical significance, total RSS, or indel timeout resolution is claimed."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

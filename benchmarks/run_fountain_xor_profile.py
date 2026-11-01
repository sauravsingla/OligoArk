"""Matched, deterministic microbenchmark of DNA Fountain droplet XOR alternatives.

This is NOT a full 100 MiB encoding benchmark or evidence of exact archive recovery.
It measures only XOR payload assembly; the existing codec/recovery gates are separate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import resource
import statistics
import time
import tracemalloc
from collections.abc import Callable
from pathlib import Path

from oligoark.ecc import xor_bytes

XorFunction = Callable[[list[bytes], int], bytes]


def integer_droplet_xor(parts: list[bytes], width: int) -> bytes:
    """Candidate implementation, with the reference's right-zero padding."""
    value = 0
    for chunk in parts:
        value ^= int.from_bytes(chunk, "big") << (8 * (width - len(chunk)))
    return value.to_bytes(width, "big")


def _cases(seed: int, operations: int, width: int) -> list[list[bytes]]:
    rng = random.Random(seed)
    chunks = [rng.randbytes(width) for _ in range(512)]
    chunks.extend((b"", rng.randbytes(width // 2), rng.randbytes(width - 1)))
    degrees = (1, 2, 3, 5, 8, 13, 21, 34)
    cases = [
        rng.sample(chunks, degrees[index % len(degrees)])
        for index in range(operations)
    ]
    # Explicitly exercise empty inputs and truncated final chunks.
    cases[:4] = [[], [b""], [chunks[-1]], [chunks[-2], chunks[0]]]
    return cases


def _apply(fn: XorFunction, cases: list[list[bytes]], width: int) -> list[bytes]:
    return [fn(parts, width) for parts in cases]


def _traced_peak(fn: XorFunction, cases: list[list[bytes]], width: int) -> int:
    tracemalloc.start()
    try:
        _apply(fn, cases, width)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return peak


def profile(seed: int, operations: int, width: int, repeats: int) -> dict[str, object]:
    if operations < 64 or width < 1 or width > 128 or repeats < 3:
        raise ValueError("require operations >= 64, 1 <= width <= 128, repeats >= 3")

    cases = _cases(seed, operations, width)
    expected = _apply(xor_bytes, cases, width)
    actual = _apply(integer_droplet_xor, cases, width)
    if actual != expected:
        raise AssertionError("integer XOR disagrees with the reference")

    functions: dict[str, XorFunction] = {
        "reference": xor_bytes,
        "integer": integer_droplet_xor,
    }
    # Warm both paths before timing; alternate order to limit first-run bias.
    for fn in functions.values():
        if _apply(fn, cases[:64], width) != expected[:64]:
            raise AssertionError("warm-up parity mismatch")
    samples: dict[str, list[float]] = {"reference": [], "integer": []}
    for repeat in range(repeats):
        order = ("reference", "integer") if repeat % 2 == 0 else ("integer", "reference")
        for name in order:
            start = time.perf_counter()
            output = _apply(functions[name], cases, width)
            elapsed = time.perf_counter() - start
            if output != expected:
                raise AssertionError(f"{name} changed outputs during timing")
            samples[name].append(elapsed)

    reference_median = statistics.median(samples["reference"])
    integer_median = statistics.median(samples["integer"])
    peak_reference = _traced_peak(xor_bytes, cases, width)
    peak_integer = _traced_peak(integer_droplet_xor, cases, width)
    peak_rss_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {
        "schema_version": 1,
        "seed": seed,
        "operations": operations,
        "width_bytes": width,
        "repeats": repeats,
        "case_degrees": [1, 2, 3, 5, 8, 13, 21, 34],
        "exact_xor_parity": True,
        "payload_sha256": hashlib.sha256(b"".join(expected)).hexdigest(),
        "seconds_by_repeat": samples,
        "reference_median_seconds": reference_median,
        "integer_median_seconds": integer_median,
        "reference_over_integer_median_ratio": reference_median / integer_median,
        "traced_peak_bytes": {"reference": peak_reference, "integer": peak_integer},
        "process_peak_rss_mib": round(peak_rss_kib / 1024, 3),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "git_sha": os.environ.get("GITHUB_SHA", "unknown"),
        "scope": (
            "Synthetic 32-byte-class XOR microbenchmark only; memory is traced "
            "Python allocation and whole-process RSS high-water, not isolated "
            "codec-worker peak RSS. No full encode/decode, throughput, 100 MiB "
            "recovery, 1 GiB milestone or physical DNA validation."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261101)
    parser.add_argument("--operations", type=int, default=4096)
    parser.add_argument("--width", type=int, default=32)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--min-relative-ratio", type=float, default=0.75)
    parser.add_argument(
        "--output", type=Path, default=Path("performance-results/fountain-xor.json")
    )
    args = parser.parse_args()

    result = profile(args.seed, args.operations, args.width, args.repeats)
    result["minimum_allowed_reference_over_integer_ratio"] = args.min_relative_ratio
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    ratio = float(result["reference_over_integer_median_ratio"])
    print(f"Exact XOR parity passed; reference/integer median ratio={ratio:.3f}x")
    print(f"Evidence: {args.output}")
    if ratio < args.min_relative_ratio:
        raise SystemExit(
            f"Fountain XOR relative microbenchmark regression: {ratio:.3f}x "
            f"< {args.min_relative_ratio:.3f}x"
        )


if __name__ == "__main__":
    main()

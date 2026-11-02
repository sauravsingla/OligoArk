"""Matched reference-vs-cached DNA whitening benchmark, not an archive speedup claim."""

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

from oligoark.framing import _mask, _mask_pad_int


def reference_mask(data: bytes, mask_id: int) -> bytes:
    """Byte-for-byte reference copied from the previous production implementation."""
    if not 0 <= mask_id <= 255:
        raise ValueError("mask_id must fit in one byte")
    legacy = (0, 0x55, 0xAA, 0xFF)
    if mask_id < len(legacy):
        value = legacy[mask_id]
        if value == 0:
            return data
        return bytes(byte ^ value for byte in data)
    state = (0x9E3779B9 ^ (mask_id * 0x45D9F3B)) & 0xFFFFFFFF
    out = bytearray(len(data))
    for index, byte in enumerate(data):
        state ^= (state << 13) & 0xFFFFFFFF
        state ^= state >> 17
        state ^= (state << 5) & 0xFFFFFFFF
        out[index] = byte ^ (state & 0xFF)
    return bytes(out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261102)
    parser.add_argument("--cases-per-length", type=int, default=512)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("performance-results/whitening-mask.json"),
    )
    args = parser.parse_args()
    if args.cases_per_length < 256 or args.repeats < 3:
        parser.error("cases-per-length >= 256 and repeats >= 3 are required")
    rng = random.Random(args.seed)
    # Covers all 256 masks at each width plus independent fixed-seed repeats.
    widths = (36, 96, 255)
    cases: list[tuple[bytes, int]] = []
    for width in widths:
        cases.extend(
            (rng.randbytes(width), index % 256)
            for index in range(args.cases_per_length)
        )
    expected = [reference_mask(data, mask) for data, mask in cases]
    actual = [_mask(data, mask) for data, mask in cases]
    if expected != actual:
        raise RuntimeError("cached whitening does not match reference bytes")
    sha256 = hashlib.sha256(b"".join(expected)).hexdigest()
    # Both functions run in the same interpreter on the same inputs with
    # alternating order; optimized caches are warm before timing.
    funcs = {"reference": reference_mask, "cached": _mask}
    for func in funcs.values():
        for data, mask in cases:
            func(data, mask)
    wall: dict[str, list[float]] = {"reference": [], "cached": []}
    cpu: dict[str, list[float]] = {"reference": [], "cached": []}
    for rep in range(args.repeats):
        order = ("reference", "cached") if rep % 2 == 0 else ("cached", "reference")
        for name in order:
            start_wall = time.perf_counter()
            start_cpu = time.process_time()
            for data, mask in cases:
                funcs[name](data, mask)
            cpu[name].append(time.process_time() - start_cpu)
            wall[name].append(time.perf_counter() - start_wall)

    median_wall = {name: statistics.median(vals) for name, vals in wall.items()}
    median_cpu = {name: statistics.median(vals) for name, vals in cpu.items()}
    raw_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_mib = raw_rss / (1024 * 1024 if sys.platform == "darwin" else 1024)
    result = {
        "parameters": {
            "seed": args.seed,
            "cases_per_length": args.cases_per_length,
            "widths": widths,
            "repeats": args.repeats,
            "all_256_mask_ids_per_width": True,
        },
        "environment": {"python": sys.version, "platform": platform.platform()},
        "exact_byte_parity": True,
        "reference_output_sha256": sha256,
        "wall_seconds": wall,
        "cpu_seconds": cpu,
        "median_wall_seconds": median_wall,
        "median_cpu_seconds": median_cpu,
        "reference_over_cached_wall_ratio": (
            median_wall["reference"] / median_wall["cached"]
        ),
        "reference_over_cached_cpu_ratio": (
            median_cpu["reference"] / median_cpu["cached"]
        ),
        "cache_info": {
            "maxsize": _mask_pad_int.cache_info().maxsize,
            "currsize": _mask_pad_int.cache_info().currsize,
            "hits": _mask_pad_int.cache_info().hits,
            "misses": _mask_pad_int.cache_info().misses,
        },
        "whole_process_peak_rss_mib": round(rss_mib, 3),
        "limitations": (
            "Warm-cache isolated whitening arithmetic under fixed-seed synthetic "
            "inputs, no independent noisy read traces. Whole-process RSS is shared "
            "between both implementations and not a differential RAM claim. "
            "No measured 10 MiB/100 MiB archive recovery speedup or timeout fix."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()

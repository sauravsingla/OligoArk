"""Reproducible graph and consensus microbenchmark with measured peak Python memory.

Run identical parameters against two checked-out commits and compare the JSON
artifacts; do not infer before/after improvements from a single invocation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import statistics
import time
import tracemalloc
from pathlib import Path

from oligoark.reconstruct import build_similarity_graph, multistart_trace_consensus


def reads_for(seed: int, strands: int, copies: int, length: int) -> list[str]:
    rng = random.Random(seed)
    reads: list[str] = []
    for _ in range(strands):
        original = "".join(rng.choices("ACGT", k=length))
        for _ in range(copies):
            noisy = list(original)
            position = rng.randrange(length)
            noisy[position] = rng.choice("ACGT")
            reads.append("".join(noisy))
    return reads


def measure(fn: object) -> tuple[float, int, object]:
    tracemalloc.start()
    start = time.perf_counter()
    try:
        result = fn()  # type: ignore[operator]
        elapsed = time.perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return elapsed, peak, result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--strands", type=int, default=12)
    parser.add_argument("--copies", type=int, default=4)
    parser.add_argument("--length", type=int, default=80)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, default=Path("performance-results/graph.json"))
    args = parser.parse_args()
    if min(args.strands, args.copies, args.length, args.repeats) < 1:
        parser.error("strands, copies, length and repeats must be positive")
    reads = reads_for(args.seed, args.strands, args.copies, args.length)
    cluster = reads[: args.copies]
    cases = {
        "graph": lambda: build_similarity_graph(reads, threshold=0.90),
        "consensus": lambda: multistart_trace_consensus(cluster),
    }
    results: dict[str, object] = {}
    for name, fn in cases.items():
        timings: list[float] = []
        peaks: list[int] = []
        fingerprints: set[str] = set()
        for _ in range(args.repeats):
            elapsed, peak, result = measure(fn)
            timings.append(elapsed)
            peaks.append(peak)
            fingerprints.add(hashlib.sha256(repr(result).encode()).hexdigest())
        if len(fingerprints) != 1:
            raise RuntimeError(f"{name} produced inconsistent results")
        results[name] = {
            "median_seconds": statistics.median(timings),
            "timings_seconds": timings,
            "peak_tracemalloc_bytes": max(peaks),
            "result_sha256": next(iter(fingerprints)),
        }
    payload = {
        "parameters": {
            "seed": args.seed, "strands": args.strands, "copies": args.copies,
            "length": args.length, "repeats": args.repeats,
        },
        "note": "tracemalloc reports traced Python allocations, not total process RSS",
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

# DNA-to-byte conversion performance validation

The 10 MiB codec research run [41170243308](https://github.com/sauravsingla/OligoArk/actions/runs/41170243308)
on October 31, 2026 retains 10/10 indel-low decode timeouts for both compact
hybrids and 10/10 mixed-noise decode timeouts for the original compact hybrid,
under the unchanged 75-second per-trial deadline. Deadline stack snapshots
include `dna_to_bytes`, `_mask`, and Reed–Solomon functions. These snapshots
are **not** representative CPU profiles.

The prior DNA-to-byte decoder materialized one Python integer per nucleotide
before constructing output bytes. The direct-group decoder reads four bases
at a time without that intermediate list. It preserves normalization,
left-to-right invalid-base reporting, byte output, and SHA-256 exact recovery.
This change does not alter RS, CRC, rescue ordering, or the channel deadlines.

## Reproduce

```bash
python -m pip install -e ".[dev,bench]"
pytest tests/test_dna_decode_direct.py tests/test_streaming.py
python benchmarks/run_dna_decode_profile.py --seed 20261031 --repeats 7
```

The benchmark retains the original decoder as an independent reference and
runs both on identical deterministic 152-, 600- and 1024-nucleotide inputs.
It records raw per-call timings, median timings, `tracemalloc` peak
allocations, Python/platform information and output SHA-256 fingerprints.
The quality CI job checks exact output parity and requires the direct
decoder's median traced peak to be at most 60% of the reference. Wall-clock
timing has **no CI pass/fail gate** due to shared-runner variability.

An initial local Python 3.13.5 controlled microprobe (seed 20261031,
nine alternating-order repetitions, 8000/3000/1500 calls per repetition)
observed reference/direct median time ratios of approximately 1.37x
(152 nt), 1.36x (600 nt) and 1.27x (1024 nt), with traced peak ratios of
approximately 2.94x, 5.41x and 6.28x. This was a local conversion-only
probe, not the CI script or a codec throughput measurement; use the
versioned CI artifact for reproducible confirmation.

**Limitations:** `tracemalloc` does not measure total process RSS. This
microbenchmark cannot establish end-to-end decode throughput, 10 MiB
indel/mixed timeout resolution, storage scale, or wet-lab recovery.
Any claims about those require separately controlled complete trials
with exact SHA-256 outcomes and all failures/timeouts retained.

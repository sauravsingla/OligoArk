# Reed–Solomon syndrome lookup experiment

The [10 MiB decoder CPU profile](results/decoder-cpu-profile-10mib-2026-11-02.md)
recorded hundreds of thousands of repeated RS syndrome calculations while
indel and mixed-noise candidate masks were tested. This change replaces the
per-symbol Python `_gf_mul()` calls **inside `_syndromes()` only** with small
cached immutable GF multiplication lookup tables. The evaluation points and
Horner polynomial arithmetic are unchanged, including the zero symbol.

## Reproduce and validate

```bash
python -m pip install -e ".[dev,bench]"
pytest -q tests/test_rs_syndrome_lookup.py tests/test_ecc_and_simulator.py
python benchmarks/run_rs_syndrome_profile.py \
    --cases 400 --message-length 96 --nsym 12 --repeats 7 \
    --output performance-results/rs-syndrome-lookup.json
```

The benchmark runs the previous evaluator and the new evaluator over identical
fixed-seed byte strings, checks **exact syndrome parity**, warms both paths,
alternates execution order, and retains individual wall- and CPU-time samples,
medians, process RSS and environment metadata. Repeat at other relevant
message widths and RS parity sizes before claiming a general result.
`tests/test_rs_syndrome_lookup.py` also checks randomized bytearray/bytes
parity and corrects randomized RS symbol substitutions. Existing
archive/codec tests remain required.

## Scientific limits

- The lookup tables are bounded by an LRU cache of 32 distinct `nsym`
  values. Each entry holds `nsym × 256` bytes, plus Python object overhead.
- A measured RS-syndrome microbenchmark gain **does not** imply an equal
  end-to-end reconstruction speedup. Decoder performance includes mask search,
  indel enumeration, DNA conversion and other work.
- The run does **not** demonstrate 10 MiB/100 MiB indel recovery, physical DNA
  storage or a 1 GiB matched 152-nt codec milestone.
- For the decisive experiment, alternate unprofiled base and optimized
  commits on the same runner and seeds, preserve unchanged 75 s (10 MiB) and
  300 s (100 MiB) per-trial limits, and retain every success, failure,
  timeout and exact SHA-256 result. Do not relax recovery gates.

No unmeasured whole-decoder speedup or noisy-channel recovery claim is made.

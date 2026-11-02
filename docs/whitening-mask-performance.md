# Bounded deterministic whitening-pad reuse: indel/mixed-noise experiment

## Why this path

The repeated 10 MiB decoder CPU profiles ([2026-11-02 evidence](results/decoder-cpu-profile-10mib-2026-11-02.md)) implicate mask candidate attempts and Reed–Solomon operations. PR #80 accelerated isolated RS syndrome arithmetic, but paired unprofiled [10 MiB recovery evidence](rs-lookup-e2e-validation.md) retained **all indel-low timeouts**. This experiment focuses on another repeated candidate cost: `framing._mask`, whose original Python xorshift loop rebuilt identical deterministic whitening masks every time a candidate was attempted.

The implementation retains the exact 256 IDs, xorshift constants, legacy `0x00, 0x55, 0xAA, 0xFF` byte masks, seed arithmetic, selector encoding, CRC, RS checks, and all decoding and ambiguity behavior. A whitening pad for a given `(mask_id, frame_length)` is reusable because it depends on neither payload content nor frame index. The algorithm now constructs that pad once, caches an immutable integer and XORs packed bytes using integer arithmetic. `lru_cache(maxsize=1024)` bounds the number of small cached entries; input lengths above 255 bytes **bypass the cache**. Byte outputs must be bit-for-bit identical, so existing archives remain compatible.

## Tests and reproducible measurement

```bash
python -m pip install -e ".[dev,bench]"
pytest -q tests/test_whitening_mask_cache.py tests/test_extended_inline_masks.py \
  tests/test_indel_candidate_streaming.py
python benchmarks/run_whitening_mask_profile.py \
  --seed 20261102 --cases-per-length 512 --repeats 7 \
  --output performance-results/whitening-mask.json
```

The randomized parity gate exercises all 256 masks across 13 byte widths, empty inputs, 255/256-byte boundary, double-whitening reversibility, cache bounds and full SHA-256 archive recovery. The dedicated [whitening microbenchmark workflow](../.github/workflows/whitening-mask-profile.yml) verifies exact output bytes for all IDs at 36/96/255 bytes, uses fixed-seed inputs and seven alternating-order warm-cache repetitions, records wall/CPU samples and cache occupancy, and retains evidence as a research artifact. A generous 0.75× relative performance *regression* bound is not a claim of improved real-world decoding.

The [paired 10 MiB workflow](../.github/workflows/whitening-mask-e2e.yml) independently compares:

- Baseline `bce161496bdee8465d9805091a53ef28c5c177ad` (before this whitening change).
- The exact candidate pull-request head SHA (or dispatch SHA when manually triggered).
- Two isolated Python environments on the same runner, alternating baseline/optimized order on two repetitions, for both compact-hybrid codecs.
- The same deterministic 10 MiB payload, 152-nt ceiling, 25% nominal redundancy, two fixed seeds (20260000/20260001) per repetition and unchanged **75-second per-trial timeout** under indel-low and mixed conditions.
- Complete worker JSON/stderr/exit codes/progress, SHA-256 exact verdicts, timeouts, encoding/decoding times and peak RSS. The existing strict `compare_rs_lookup_e2e.py` validator checks that the encoded output, trial accounting and payload fingerprints match between revisions and that missing workers are not misreported as recovery failures.

This is a software simulation and a small experiment. Repeated seeds are **not** independent observations. An isolated mask calculation ratio must never be presented as end-to-end reconstruction speedup, an indel recovery fix, improved physical DNA reliability or 1 GiB matched-codec robustness. The original 100 MiB DNA Fountain encoding-watchdog failure remains recorded and unrelated. If paired outcomes remain timed out, document that negative result rather than relaxing deadlines or hiding failures.


## Verified measurements — 2 November 2026

### Isolated whitening microbenchmark

[GitHub run 41380962191](https://github.com/sauravsingla/OligoArk/actions/runs/41380962191), [retained artifact 13112772694](https://github.com/sauravsingla/OligoArk/actions/runs/41380962191/artifacts/13112772694), Python 3.12.15 on one Linux GitHub runner, seed 20261102. Both versions processed the same 512 cases per input length (36, 96, 255 bytes), covering all 256 masks at each length; seven alternating-order repetitions ran with the optimized pad cache warm. Exact byte parity passed.

| Isolated metric | Frozen original | Cached optimized |
| --- | ---: | ---: |
| Median wall time for the matched input batch | 65.776 ms | 1.026 ms |
| Median CPU time for the matched input batch | 65.774 ms | 1.025 ms |

The observed reference/optimized median ratio was **64.079× wall time** and **64.172× CPU time**. The bounded immutable cache used **765 of 1,024 entries** with 13,005 hits. Whole-process sampled peak RSS was **23.457 MiB**, shared between both implementations and **not** a differential memory comparison. These are *whitening operations only*, not archive throughput, indel correction, statistically significant improvements or wet-lab DNA storage results.

### Matched unprofiled 10 MiB before/after recovery

[GitHub run 41380962156](https://github.com/sauravsingla/OligoArk/actions/runs/41380962156) passed both paired research jobs. I independently inspected the raw [compact artifact](https://github.com/sauravsingla/OligoArk/actions/runs/41380962156/artifacts/13115536246) and [compact-v3 artifact](https://github.com/sauravsingla/OligoArk/actions/runs/41380962156/artifacts/13114469643), including all worker JSON, stderr, exit codes, progress streams and paired reports. All eight workers exited 0, emitted exactly two conditions with two trials each and ended with a `worker-complete` marker; stderr was empty. Baseline `bce161496bdee8465d9805091a53ef28c5c177ad` and optimized PR #83 head `f76fb96caab57d7fb03173143ac711208b6a3ea8` produced the same deterministic payload SHA-256 `74e1f0601eb2381d8a6669aadea726dcbfd44cbdce7e0873ae53642e773568bd` and matched encoded-strand fingerprints. The optimized production change was merged as [PR #83](https://github.com/sauravsingla/OligoArk/pull/83), merge commit `e8a263062a025acf22f151e8305b02e81c543f7b`.

| Codec | Software noise | Baseline exact recoveries (timeouts) | Optimized exact recoveries (timeouts) | Median successful-only decode (before → after) |
| --- | --- | ---: | ---: | --- |
| Compact hybrid | Indel-low | 0/4 (4) | 0/4 (4) | Not available |
| Compact hybrid | Mixed | 0/4 (4) | 0/4 (4) | Not available |
| Compact-v3 hybrid | Indel-low | 0/4 (4) | 0/4 (4) | Not available |
| Compact-v3 hybrid | Mixed | **4/4 (0)** | **4/4 (0)** | **16.100 s → 13.960 s** |

The experiment used two distinct trial seeds (20260000/20260001) observed twice each. The four observations per cell are **not four independent seeds**. The per-trial timeout remained 75 seconds, with 152-nt strand ceiling and nominal 25% redundancy unchanged. A green research workflow denotes complete accounting, **not** universal noisy recovery.

Two matched compact-hybrid encoding times were **139.424 and 139.521 s baseline**, compared with **120.002 and 117.625 s optimized**. Compact-v3 encoded in **67.392 and 66.771 s baseline**, compared with **64.756 and 63.804 s optimized**. Successful compact-v3 mixed-noise median decode was about 13.3% lower in the optimized observations. Those descriptive differences are not a statistically established end-to-end improvement on independent seeds or hardware. Worker RSS was roughly 300–357 MiB for compact hybrid and 372–411 MiB for compact-v3 across channels; no memory reduction is claimed.

**Negative results remain:** every indel-low observation timed out for both codecs, and compact-hybrid mixed also timed out. PR #83 did not resolve the indel bottleneck. Subsequent work should profile and measure other candidate/RS/realignment stages using independent trial seeds, compare with unchanged deadlines and preserve every exact SHA-256 verdict and timeout. Do not infer physical DNA reliability or 1 GiB matched-noisy-codec recovery from these software experiments.

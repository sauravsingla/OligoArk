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

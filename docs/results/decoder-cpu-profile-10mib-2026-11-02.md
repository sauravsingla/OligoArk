# Matched 10 MiB decoder CPU profiles: retained negative results (2026-11-02)

## Provenance

This report transcribes the **2026-11-01** PR #78 [focused profiling run](https://github.com/sauravsingla/OligoArk/actions/runs/41344124815) and inspects all eight archived `cProfile` files. The source PR was merged as [#78](https://github.com/sauravsingla/OligoArk/pull/78) after its 13 GitHub check runs passed. The experiment ran at PR head `426bdfb7ad6f8f367900614875db2d7f828118f9` (GitHub Actions checkout records a synthetic PR merge SHA `d20d55182e82c90977ef6b6ec22e07c9aeede1b2` in the report). The codec and noise model were unchanged.

- [Original compact artifact](https://github.com/sauravsingla/OligoArk/actions/runs/41344124815/artifacts/13098781747), ZIP SHA-256 `26ca3a30fe294f47021a13692f1c7a83d79d3588b2164afa96c4ea8ed1f0faf2`.
- [Compact-v3 artifact](https://github.com/sauravsingla/OligoArk/actions/runs/41344124815/artifacts/13097474938), ZIP SHA-256 `01d6fe18fc6e688651d2036ddaa99a51be7ec2c07a465a8ca96d89ac13163305`.

The `summary.json` in each artifact records individual `.prof` SHA-256 hashes; all eight downloaded profiles were independently checked against their recorded hashes. The ZIPs also retain per-trial `top30.txt`, `rows.json`, `progress.jsonl`, and summaries.

## Matched experiment and observed outcomes

Both methods used deterministic 10,485,760-byte input, trial seeds 20260000 and 20260001, the same 152-nt strand ceiling, nominal 25% redundancy, `indel-low` and `mixed` software channels, and an unchanged **75-second per-trial deadline**. The workflow preserves each failed trial; a green *profiling workflow* means complete accounting and artifact retention, **not successful DNA recovery**.

| Method | Channel | SHA-256 verified successes | Timeouts | Failures | Encode seconds (shared per method) | Peak worker RSS MiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Original compact hybrid | indel-low | 0/2 | 2/2 | 2/2 | 137.933 | 299.242 |
| Original compact hybrid | mixed | 0/2 | 2/2 | 2/2 | 137.933 | 345.930 |
| Compact-v3 hybrid | indel-low | 0/2 | 2/2 | 2/2 | 89.512 | 372.285 |
| Compact-v3 hybrid | mixed | 2/2 | 0/2 | 0/2 | 89.512 | 410.484 |

The two successful *profiled* compact-v3 mixed trials reported a mean decode time of 55.053 seconds. The `mean_decode_seconds=0` field in timeout-only rows is **not a zero-latency decode**; no completed decode time is available for those rows. Encoding seconds are reused across the two conditions of each method and must not be added as independent encodes. These are **profiled** measurements, not unprofiled speed or recovery comparisons.

## Repeated CPU hotspots

Across the four original-compact profiles, `framing.decode_frame`, `framing._decode_compact_with_mask`, `ecc.rs_decode`, `ecc._syndromes` and `ecc._poly_eval` recur among the highest cumulative-time functions. For example, original-compact indel-low trial 0 recorded ~779,601 calls to `_decode_compact_with_mask`, ~779,601 calls to `rs_decode`, ~844,113 calls to `_syndromes`, and ~1.75 million calls to `_poly_eval`. `framing._mask`, `dna.dna_to_bytes`, and `archive.rescue` also recur.

Compact-v3 indel-low trial 0 recorded ~610,780 `_decode_compact_inline_with_mask` calls, ~610,780 `rs_decode` calls, ~650,630 `_syndromes` calls, and ~377,034 `_inline_decode_mask_candidates` calls. Both compact-v3 mixed profiles similarly show repeated `decode_frame_packed`, `_inline_decode_mask_candidates`, `rs_decode`, `_syndromes` and `_poly_eval` activity. These data point to **repeated candidate/mask and Reed–Solomon decoding work** as an experiment-worthy cost center; they do not prove a particular safe algorithm change.

The `cProfile` cumulative totals are **inclusive** (nested callers overlap); do not add them to estimate total CPU. The interrupted profiles may contain partially flushed or instrumented timing. The counts above are observations from individual fixed-seed traces, not general throughput estimates.

## Interpretation, limitations and next validation

1. **No speedup or timeout fix is claimed.** `cProfile` adds overhead and can change deadline outcomes. The earlier [unprofiled 10 MiB codec run 41311007188](https://github.com/sauravsingla/OligoArk/actions/runs/41311007188) separately recorded 10/10 original-compact indel-low timeouts, 10/10 original-compact mixed timeouts and 6/10 compact-v3 indel-low timeouts. Its outcomes are not numerically interchangeable with these two-trial profiled runs.
2. First isolate the repeated mask/RS attempts in a small deterministic, byte-for-byte parity microbenchmark with fixed noisy inputs. Evaluate a candidate only if it preserves successful and failed decode semantics, candidate order where relevant, exception behavior, SHA-256 archive recovery, and randomized parity. Run alternating-order, unprofiled before/after repeats on the same runner, with median latency, throughput, traced allocations and peak RSS. Reject changes that do not measurably help.
3. Follow with **unprofiled** matched 10 MiB indel-low and mixed trials, unchanged seeds/channel rates, 75-second deadlines, full timeout and SHA-256 accounting, then larger 100 MiB trials only within available runner resources. Do not relax deadlines or suppress failed trials to obtain green scientific results.
4. Parallel reconstruction, caching, graph search and adaptive decoding are **not yet justified by these profiles** without isolated parity and resource measurements. No production-code optimization is included in this report.
5. The separate [100 MiB research run 41270450630](https://github.com/sauravsingla/OligoArk/actions/runs/41270450630) retains a 2,400-second DNA Fountain encoding watchdog failure. The [1 GiB software streaming acceptance](storage-scale-1gib-2026-11-01.md) is a different relaxed-constraint workload and is not evidence of 152-nt noisy matched-codec or wet-lab recovery.

### Reproduction

```bash
python -m pip install -e ".[dev,bench]"
pytest -q tests/test_codec_decoder_cprofile.py tests/test_matched_codec_runtime.py tests/test_streaming.py
export OLIGOARK_CODEC_CPROFILE_DIR="$PWD/decoder-profile-results/profiles"
export OLIGOARK_CODEC_PROGRESS_FILE="$PWD/decoder-profile-results/progress.jsonl"
python benchmarks/run_matched_codec_scale.py --worker --profile research-10 \
  --method oligoark-compact-hybrid --size 10485760 \
  --conditions indel-low,mixed --trials 2 \
  --condition-timeout-seconds 75 > decoder-profile-results/rows.json
```

Repeat with `--method oligoark-compact-v3-hybrid` in a separate output directory. Preserve `rows.json`, partial `.prof` files, progress and RSS on deadline exceptions. For the canonical two-method matrix and evidence checks, use [the merged workflow](../../.github/workflows/decoder-cpu-profile.yml). CPU profiles are diagnostic and do not certify exact recovery unless the corresponding trial has `sha256_verified=true`.

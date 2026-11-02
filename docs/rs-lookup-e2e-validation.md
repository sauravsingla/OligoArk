# Matched unprofiled 10 MiB RS lookup end-to-end research

The syndrome arithmetic experiment in [PR #80](https://github.com/sauravsingla/OligoArk/pull/80) found an isolated warm-cache CPU/wall improvement. That result is **not** an end-to-end decoding improvement. This workflow measures exact archive recovery under indel-low and mixed noise without cProfile, holding the input and channel conditions fixed.

## Reproducible design

- **Before**: f2dc56e2ff0a845794271a46345976ae75af8da6 (immediate parent of the PR #80 squash merge).
- **After**: 8af0391fe0d0988b22d13be564a791ef4721dad0 (merged PR #80). The only production-code change between these revisions is the RS syndrome lookup in src/oligoark/ecc.py; the other changes are tests, workflow, benchmark and documentation.
- Two isolated Python virtual environments on the **same GitHub runner** per codec, Python 3.13, and two alternating-order repetitions: before/after then after/before. Separate worker processes prevent cross-revision Python caches from contaminating measurements.
- 10,485,760 deterministic input bytes, same SHA-256 payload, 152-nt strand ceiling, nominal 25% redundancy, indel-low and mixed software channels, and trial seeds 20260000 and 20260001. **These seeds repeat across timing repetitions**, so four observations per condition/version are *not four independent recovery trials*.
- Every trial has an unchanged **75-second hard deadline**, including channel simulation and decoding. The outer 900-second watchdog protects the runner during encoding/worker orchestration and does **not** change the per-trial deadline. A 900-second worker interruption is an *incomplete experiment*, not a verified negative recovery result.
- Raw worker JSON, stderr, exit codes and worker-side progress/RSS diagnostics are retained even if a worker fails. A validator rejects incomplete, missing, duplicate or unavailable outcomes; true decode failures and per-trial timeouts remain in the denominators. Green workflow status indicates complete evidence, not successful DNA recovery.
- Reports include exact recoveries, failures, timeouts, encoding seconds, successful-only decoding times, and worker peak RSS. The successful-only median **cannot** be used as overall throughput when success rates differ. The workflow is report-only: no speedup gate, no weakened SHA-256 acceptance rule, and no changes to production codec logic.

## Running and interpreting

The workflow [.github/workflows/rs-lookup-e2e-comparison.yml](../.github/workflows/rs-lookup-e2e-comparison.yml) runs on changes to itself or its validator, and can also be dispatched manually. Each method produces a 90-day artifact containing paired results or incomplete-worker diagnostics. The validator is benchmarks/compare_rs_lookup_e2e.py, with negative-evidence tests in tests/test_compare_rs_lookup_e2e.py.

A small fixed-seed experiment cannot establish statistical significance or generalize to all noise distributions. Review raw results and run additional independent seeds and runner repetitions before making a performance claim. Follow with resource-bounded **100 MiB** work using the existing **300-second per-trial** deadline, only after the 10 MiB evidence is evaluated. The separate 1 GiB software-streaming milestone is not evidence of realistic noisy DNA recovery. External CNR comparisons and wet-lab validation remain separate research questions.

## Observed 10 MiB results (2026-11-02, PR #81)

[Workflow run 41370771730](https://github.com/sauravsingla/OligoArk/actions/runs/41370771730) completed both unprofiled, alternating-order matrix jobs. We inspected the two raw artifact ZIPs, not only the green Actions status: all eight workers exited 0, each recorded two condition rows and four trial verdicts, and every progress stream ended with `worker-complete`. All eight stderr files were empty. Each result uses the same deterministic 10 MiB payload SHA-256, `74e1f0601eb2381d8a6669aadea726dcbfd44cbdce7e0873ae53642e773568bd`. The 75 s per-trial limit, 152 nt ceiling and 25% nominal redundancy were unchanged.

| Codec | Channel | Before #80 exact / observations (timeouts) | After #80 exact / observations (timeouts) | Successful-only decode median, before → after |
| --- | --- | ---: | ---: | ---: |
| compact-hybrid | indel-low | 0/4 (4) | 0/4 (4) | Not available → not available |
| compact-hybrid | mixed | 0/4 (4) | 0/4 (4) | Not available → not available |
| compact-v3-hybrid | indel-low | 0/4 (4) | 0/4 (4) | Not available → not available |
| compact-v3-hybrid | mixed | 4/4 (0) | 4/4 (0) | 24.450 s → 22.484 s |

Each row contains **two fixed seeds observed twice**, not four independent samples. The 4/4 exact recoveries in compact-v3 mixed are SHA-256-verified. All other rows were genuine per-trial timeouts, **not** successful reconstructions or missing worker output. Green research jobs mean complete, validated evidence, **not** that noisy recovery passed.

The compact-hybrid encoder took 104.016–107.456 s across the before/after repetitions; compact-v3-hybrid took 86.568–89.466 s. Observed worker peak RSS ranged 305.031–306.969 MiB (compact indel-low), 356.023–356.973 MiB (compact mixed), 371.871–373.551 MiB (v3 indel-low), and 409.973–411.672 MiB (v3 mixed) across both revisions. These ranges are descriptive, not paired memory-reduction claims.

**Interpretation:** The observed lower successful-only median for v3 mixed is not evidence of a statistically established end-to-end speedup. Two independent seeds and two repetitions on one runner are insufficient to establish generalizable throughput or latency improvement; all indel-low cases and compact mixed remain limited by the original 75 s deadline. The isolated ~2.2× syndrome microbenchmark in #80 must not be reported as an end-to-end gain. Before attempting another production optimization, use targeted profiling of the timed-out decode path and more independent seeds, without increasing the original deadline or concealing failures.

Raw data, including worker JSON, progress, exit codes, stderr, RSS and paired reports:
- [Compact-hybrid artifact](https://github.com/sauravsingla/OligoArk/actions/runs/41370771730/artifacts/13110076878)
- [Compact-v3-hybrid artifact](https://github.com/sauravsingla/OligoArk/actions/runs/41370771730/artifacts/13110131695)

The artifacts are configured for 90-day retention. The [PR #81](https://github.com/sauravsingla/OligoArk/pull/81) harness and its validator remain reproducible via manual workflow dispatch. Separate [100 MiB research](https://github.com/sauravsingla/OligoArk/actions/runs/41270450630) retained a DNA Fountain 2400 s watchdog failure and must not be represented as passing. The CI storage profile is not a 1 GiB exact-recovery milestone.

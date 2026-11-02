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

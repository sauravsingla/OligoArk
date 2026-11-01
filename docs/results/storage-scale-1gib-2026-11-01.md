# 1 GiB software streaming-storage acceptance — 2026-11-01

## Provenance and reproducibility

The scheduled **Storage scale acceptance** run [41292938682](https://github.com/sauravsingla/OligoArk/actions/runs/41292938682) completed successfully on 2026-11-01 at commit [`a6463d9930266d94dfc15c4102dfe97c6f1fa8dd`](https://github.com/sauravsingla/OligoArk/commit/a6463d9930266d94dfc15c4102dfe97c6f1fa8dd). The run used Python 3.13.17 on `Linux-7.0.0-1017-azure-x86_64-with-glibc2.39` and the `storage` profile:

```bash
python -m pip install -e ".[bench]"
python benchmarks/run_storage_scale.py --profile storage
```

Source: [storage-scale-41292938682 artifact](https://github.com/sauravsingla/OligoArk/actions/runs/41292938682/artifacts/13078815873), containing `metadata.json`, `results.json`, `results.csv`, `summary.json` and plots. SHA-256 of the downloaded artifact ZIP: `42cdde72beeb610369a4a613ed030ae9e8fffa078f96509229e7cd90d2ddf027`; SHA-256 of its uncompressed `results.json`: `564c1dbf4e25c806365b41b56810a0f58a93613b2611ff41f64a49ff3fac5b92`. Check the original artifact and these hashes before reusing the measurements. This is one environment and one run, **not** a controlled before/after speedup comparison.

## Verified 1 GiB exact recovery

The `scale-1024` software streaming profile used a deterministic heterogeneous **1,073,741,824-byte** payload and XOR parity. All three channel cases returned the **same SHA-256** `cc6286341b8650694bdb4f565a75372a981fde3935d1577eb415031d9fb12b5e` and set `sha256_verified=true`. The reported 1 GiB milestone is therefore **achieved for software archive streaming under clean and controlled strand erasure**.

| Channel | Exact SHA-256 | Erased/recovered data strands | Encode s | Decode s | Encode MiB/s | Decode MiB/s | Peak RSS MiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Clean | PASS | 0 / 0 | 98.562 | 44.096 | 10.389 | 23.222 | 43.996 |
| Controlled 1% dropout | PASS | 45,338 / 45,338 | 97.352 | 57.522 | 10.518 | 17.802 | 44.031 |
| Controlled 5% dropout | PASS | 226,705 / 226,705 | 99.102 | 85.371 | 10.333 | 11.995 | 43.992 |

The 1 GiB container had **4,530,557 data strands**, **566,320 parity strands**, **5,096,877 total strands**, **5,219,201,308 encoded nucleotides**, **1.645833 logical bits/nt**, **12.5% strand redundancy**, and **1,325,188,633 archive bytes** (1.234178x source size). The measured interpreter baseline RSS was 25.762 MiB, with approximately 18.23–18.27 MiB additional peak RSS across the three 1 GiB cases. These values describe the streaming worker, not total filesystem cache or aggregate storage footprint.

Across the fixed clean XOR scale sweep (1 KiB, 64 KiB, 1 MiB, 10 MiB, 100 MiB, 1 GiB), `summary.json` reports a **0.043306** log-log peak-RSS slope (classified `bounded`) and a **0.776126** runtime slope (classified `linear-or-better`). Those are measurements on the selected workloads, not universal complexity proofs or comparisons with another revision.

## Negative results and scope boundaries

The separate 64 KiB `oligoark-248-compact` software fault matrix retained **8 failures out of 24 rows**, including unprotected `none` and `fountain` dropout/substitution/mixed cases. The green workflow does **not** mean every fault-matrix condition recovered exactly. The 1 GiB acceptance rows use `scale-1024`, `rs_nsym=0`, with biochemical constraint checks intentionally disabled (reported sampled constraint pass rate 0); they do **not** demonstrate realistic 152-nt physical-strand constraints, insertion/deletion tolerance, mixed-noise robustness, DNA synthesis/sequencing, or end-to-end wet-lab file recovery.

This 1 GiB **streaming storage** success is distinct from the matched **152-nt codec** research benchmark. The [2026-11-01 codec run 41270450630](https://github.com/sauravsingla/OligoArk/actions/runs/41270450630) still has unresolved indel-low/mixed-noise recovery failures and a clean-room DNA Fountain 100 MiB encoding watchdog failure. See [codec worker evidence](../codec-worker-progress.md) and [decoder CPU profiling instructions](../decoder-cpu-profiling.md). Do not infer a 1 GiB noisy matched-codec result or an end-to-end codec speedup from this streaming-storage milestone.

# Fast external physical-read benchmark

This benchmark evaluates the current OligoArk **confidence-fusion** reconstruction method on Microsoft's public **Clustered Nanopore Reads (CNR)** physical DNA-storage dataset and compares it with the pinned external **Bidirectional Beam Search (BBS)** baseline.

## Dataset and protocol

The benchmark uses the CNR dataset released with *Trellis BMA: coded trace reconstruction on IDS channels for DNA storage*. It contains 10,000 explicit 110-base reference strands and 269,709 clustered Oxford Nanopore MinION reads. `Clusters.txt` maps in order to `Centers.txt`. The workflow pins dataset commit `6938f44796185902a08381943c2895782886c5c3` and verifies Git blob and SHA-256 checksums.

The external baseline is the official BBS implementation pinned at commit `3e4ab46871929819e4f3e34a831c57cac88bb456`, built with locked dependencies and run on one CPU thread.

The authoritative held-out set is the unchanged deterministic **96-cluster** subset selected from clusters with at least 10 reads using SHA-256 ranking with seed `20261005`. Confidence-fusion development used a separate deterministic split with seed `20261007`, after excluding both the 96 held-out clusters and the earlier 48 calibration clusters. The final configuration was frozen before the held-out run.

## Current held-out result

The authoritative one-shot validation is workflow run `37456317557`, artifact `11410090111`.

| Reads / strand | **OligoArk confidence fusion** | External BBS |
| ---: | ---: | ---: |
| 1 | 1/96 (1.0%) | 1–2/96 across five repeats |
| 5 | **73/96 (76.0%)** | 72–74/96 across five repeats |
| 10 | **93/96 (96.9%)** | **93/96 (96.9%)** in all five repeats |

At 10 reads, OligoArk mean edit distance was **0.0417** versus **0.3542** for BBS. OligoArk had **2 one-edit failures, 1 two-edit failure, and no 3+ edit failures**. BBS had 1 one-edit failure and 2 failures with 3+ edits.

At 5 reads, OligoArk mean edit distance was **0.4896**. BBS ranged from 72/96 to 74/96 exact recovery across its five upstream-code repeats because equal-score candidates can vary with Rust HashMap iteration.

Using the first pinned BBS repeat for paired testing, the 5-read comparison had **11 OligoArk-only exact recoveries versus 10 BBS-only recoveries** (two-sided exact McNemar p = **1.0**). At 10 reads, each method had **2 exclusive exact recoveries** (p = **1.0**). Therefore this benchmark shows **no statistically significant exact-recovery difference between OligoArk confidence fusion and BBS at either 5 or 10 reads**.

## Runtime

Confidence fusion took about **22.6 seconds** for the 96-cluster 5-read evaluation and **52.6 seconds** at 10 reads, with about **25.3 MB peak RSS**. BBS remained much faster at roughly **0.18–0.23 seconds** per 96-cluster multi-read evaluation and about **4 MB peak RSS**.

The complete one-shot GitHub Actions workflow—including focused tests, dataset verification, pinned BBS build, all required reconstruction comparisons, five BBS repeats per coverage, acceptance checks, and artifact upload—completed in about **6 minutes 31 seconds**, below the 10-minute budget.

## Development provenance

The protected development split contained **48 clusters** with zero overlap with either the held-out set or the historical calibration set. The frozen configuration `fusion-fast-t2-c8-r0-q025-g005` improved 10-read development recovery from **44/48 to 46/48**, with 3 rescues and 1 regression, while 5-read recovery improved from 34/48 to 35/48. No held-out reference was used for parameter selection.

## Claim boundary

This is a **physical-read reference-reconstruction benchmark**, not an end-to-end OligoArk archive decode. The CNR strands were not encoded with OligoArk framing, ECC, fountain/hybrid redundancy, or the combined optimizer, so those archive-level methods and original-file SHA-256 recovery are not applicable to this comparison.

The current conclusion is: **on this fixed 96-cluster physical CNR benchmark, OligoArk confidence fusion matches the pinned BBS baseline at 93/96 exact recovery with 10 reads per strand and shows no statistically significant exact-recovery difference from BBS at 5 or 10 reads, while remaining CPU-only and below the 10-minute workflow budget. BBS is still substantially faster, and this single benchmark does not establish a general state-of-the-art claim.**

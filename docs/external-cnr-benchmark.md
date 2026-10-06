# Fast external physical-read benchmark

This benchmark evaluates the current OligoArk targeted low-compute reconstruction method on Microsoft's public **Clustered Nanopore Reads (CNR)** physical DNA-storage dataset and compares it with the pinned external **Bidirectional Beam Search (BBS)** baseline.

## Dataset and protocol

The benchmark uses the CNR dataset released with *Trellis BMA: coded trace reconstruction on IDS channels for DNA storage*. It contains 10,000 explicit 110-base reference strands and 269,709 clustered Oxford Nanopore MinION reads. `Clusters.txt` maps in order to `Centers.txt`. The workflow pins dataset commit `6938f44796185902a08381943c2895782886c5c3` and verifies Git blob and SHA-256 checksums.

The external baseline is the official BBS implementation pinned at commit `3e4ab46871929819e4f3e34a831c57cac88bb456`, built with locked dependencies and run on one CPU thread.

The held-out evaluation uses a fixed deterministic **96-cluster** subset selected from clusters with at least 10 reads using SHA-256 ranking with seed `20261005`. A separate deterministic **48-cluster calibration set** uses seed `20261006` after excluding all held-out cluster IDs. The targeted repair configuration was frozen from calibration before the held-out run.

## Current held-out result

The authoritative one-shot held-out validation is workflow run `37412327238`, artifact `11390090747`.

| Reads / strand | **OligoArk targeted repair** | External BBS |
| ---: | ---: | ---: |
| 5 | **68/96 (70.8%)** | 73/96 (76.0%) |
| 10 | **84/96 (87.5%)** | 93/96 (96.9%) |

At 10 reads, OligoArk achieved mean edit distance **0.167**, with **8 one-edit failures, 4 two-edit failures, and no 3+ edit failures**. At 5 reads, mean edit distance was **0.594**.

Against BBS, the 5-read exact-recovery difference was not statistically significant on this subset (two-sided exact McNemar p = **0.424**). At 10 reads, BBS remained significantly stronger on exact recovery (p = **0.0117**).

## Runtime

The targeted method took about **17.1 seconds** for the 96-cluster 5-read evaluation and **41.2 seconds** at 10 reads, with about **24.4 MB peak RSS**. The complete one-shot GitHub Actions workflow—including dataset verification, pinned BBS build, all required comparisons, repeated BBS measurements, acceptance checks, and artifact upload—finished in roughly **5 minutes**, within the low-compute target.

## Calibration result

The frozen targeted configuration reached **34/48 (70.8%)** exact recovery at 5 reads and **47/48 (97.9%)** at 10 reads on the separate calibration set. The calibration and held-out cluster sets have zero overlap.

## Claim boundary

This is a **physical-read reference-reconstruction benchmark**, not an end-to-end OligoArk archive decode. The CNR strands were not encoded with OligoArk framing, ECC, fountain/hybrid redundancy, or the combined optimizer, so those archive-level methods and original-file SHA-256 recovery are not applicable to this comparison.

The current conclusion is: **OligoArk's targeted low-compute repair reaches 87.5% exact reconstruction at 10 reads on the held-out CNR subset while remaining practical on a standard CPU runner, but BBS remains stronger at 96.9% exact recovery.**

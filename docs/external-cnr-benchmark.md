# Fast external physical-read benchmark

This benchmark evaluates OligoArk reconstruction on one public physical DNA-storage dataset and compares it with a recent external trace-reconstruction baseline while staying practical on a standard GitHub-hosted CPU runner.

## Dataset and external baseline

The benchmark uses Microsoft's **Clustered Nanopore Reads (CNR)** dataset released with *Trellis BMA: coded trace reconstruction on IDS channels for DNA storage* (ISIT 2021). It contains 10,000 explicit 110-base reference strands and 269,709 clustered Oxford Nanopore MinION reads. The clusters in `Clusters.txt` map in order to `Centers.txt`. The workflow pins dataset commit `6938f44796185902a08381943c2895782886c5c3` and verifies both Git blob and SHA-256 checksums.

The external baseline is **Bidirectional Beam Search (BBS)** from Gu et al., *Efficient trace reconstruction in DNA storage systems using bidirectional beam search* (iScience, 2025). The workflow builds the official Rust implementation at commit `3e4ab46871929819e4f3e34a831c57cac88bb456` with locked dependencies and uses one CPU thread. HEDGES is not included because the CNR strands were not encoded with HEDGES; retrofitting a codec would change the physical experiment rather than compare reconstruction on identical reads.

The CNR maintainers note that the generated centers contain unintended long-range dependencies and that some recovered clusters may be malformed. This limitation is retained.

## Leakage-controlled low-compute improvement

The original **96-cluster held-out set is unchanged**: clusters with at least 10 reads are SHA-256 ranked with seed `20261005`, and the same nested 1-, 5-, and 10-read subsets are used for every method.

A separate **48-cluster calibration set** uses seed `20261006` after all 96 held-out cluster IDs have been excluded. Five fixed low-compute multi-start configurations are evaluated only on calibration references at 5 and 10 reads. Candidate runtime is limited to at most `max(4 × iterative baseline, baseline + 2 s)`. The winner is frozen before the held-out references are scored.

The new `multistart_trace_consensus()` method ranks a small number of observed reads by agreement with the other reads, performs the existing alignment refinement from those anchors, optionally repeats refinement on reversed traces, and selects the final candidate by observed-read edit distance plus a known-length penalty. It does **not** use the unknown reference, deep learning, a GPU, BBS source code, or a new runtime dependency.

Calibration selected **3 anchors, 1 refinement round, bidirectional mode, length penalty 1.0**. Across the 48 calibration clusters at 5 and 10 reads, this configuration recovered **68/96 (70.8%)** versus **53/96 (55.2%)** for the previous iterative trace baseline. Calibration runtime was 43.1 s versus the iterative-trace calibration baseline, within the predeclared 4× budget.

## Executed held-out result

The selected held-out 10-read subset contains 960 physical reads. Its descriptive raw-read profile measured about 2.23% substitutions, 1.65% insertions and 1.92% deletions per aligned reference base.

| Reads / strand | Direct | Medoid | Graph/alignment | Iterative trace | **Multi-start trace** | External BBS |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1/96 (1.0%) | 1/96 (1.0%) | 1/96 (1.0%) | 1/96 (1.0%) | 1/96 (1.0%) | 1–2/96 across repeats |
| 5 | 1/96 (1.0%) | 6/96 (6.3%) | 33/96 (34.4%) | 44/96 (45.8%) | **53/96 (55.2%)** | **72–73/96 (75.0–76.0%)**; first repeat 73/96 |
| 10 | 1/96 (1.0%) | 10/96 (10.4%) | 65/96 (67.7%) | 63/96 (65.6%) | **80/96 (83.3%)** | **93/96 (96.9%)** in all five repeats |

The multi-start method significantly improved exact recovery over the previous OligoArk methods on the same held-out clusters. At 5 reads it had 11 exact-only wins versus 2 iterative-only wins (two-sided exact McNemar p = 0.0225), and 20 exact-only wins versus 0 graph/alignment-only wins (p = 1.91e-6). At 10 reads it had 17 exact-only wins versus 0 iterative-only wins (p = 1.53e-5), and 15 exact-only wins versus 0 graph/alignment-only wins (p = 6.10e-5).

BBS remains stronger on exact recovery. At 5 reads, the frozen multi-start method had 4 exact-only wins versus 24 BBS-only wins (p = 1.80e-4). At 10 reads it had 1 exact-only win versus 14 BBS-only wins (p = 0.000977). The exact-recovery gap therefore narrowed substantially but did not close.

Mean edit distance at 5 reads was **0.740** for multi-start, 1.010 for iterative trace, 1.938 for graph/alignment, and 0.625 for the first BBS repeat. At 10 reads it was **0.219** for multi-start, 0.469 for iterative trace, 0.490 for graph/alignment, and 0.354 for BBS. The lower 10-read mean edit distance for multi-start does not make it superior overall: BBS still has substantially higher exact-reference recovery, while its few failures contain more edits.

On the GitHub-hosted runner, multi-start took about **23.8 s** for the 96-cluster 5-read evaluation and **61.9 s** at 10 reads, with about **24 MB peak RSS**. The complete workflow—including dataset checks, pinned BBS build, calibration grid, all held-out methods, five BBS repetitions per coverage, validation, and artifact upload—completed in **6 min 58 sec**, within the 10–15 minute target.

## Metrics and artifacts

The workflow reports exact-reference recovery with Wilson 95% intervals, edit and normalized-edit distance, wall time, peak RSS, deterministic subset hashes, paired exact McNemar tests versus BBS and prior OligoArk methods, calibration candidates, raw per-cluster CSV/JSON, BBS repeat variability, timing records, and provenance metadata. Failures are retained.

## Claim boundary

This is a **physical-read reference-reconstruction benchmark**, not an end-to-end OligoArk archive decode. CNR does not contain a file encoded with OligoArk framing, ECC, fountain/hybrid redundancy, or the combined optimizer, so those archive-level results and original-file SHA-256 recovery are not applicable.

The defensible conclusion is: **the lightweight multi-start refinement materially improves OligoArk on unseen physical nanopore reads and narrows the gap to BBS, especially at 10-read coverage, but OligoArk still underperforms BBS on exact reconstruction and is not state of the art on this benchmark.**

## Targeted one-edit repair: one-shot held-out result

After the earlier multi-start result, parameter development moved to the disjoint 48-cluster calibration split only. The frozen targeted configuration was selected before held-out evaluation from calibration workflow run `37411358445`, where it reached **47/48 (97.9%)** exact recovery at 10 reads and **34/48 (70.8%)** at 5 reads. The original 96 held-out cluster IDs and deterministic read ordering were unchanged.

A single held-out validation was then executed in workflow run `37412327238` (artifact `11390090747`). The bounded targeted repair achieved:

| Reads / strand | Multi-start | **Targeted one-edit repair** | External BBS |
| ---: | ---: | ---: | ---: |
| 5 | 53/96 (55.2%) | **68/96 (70.8%)** | 73/96 (76.0%) |
| 10 | 80/96 (83.3%) | **84/96 (87.5%)** | 93/96 (96.9%) |

At 10 reads, targeted repair reduced mean edit distance from **0.219** to **0.167** and left **8 one-edit failures, 4 two-edit failures, and no 3+ edit failures**. Against multi-start it had 10 exact-only wins and 6 regressions (two-sided exact McNemar p = 0.4545), so the +4 exact recoveries are a real observed improvement but are not statistically significant at this sample size. Against iterative trace and graph/alignment, the targeted method improved exact recovery significantly (p = 0.000324 and p = 0.000878 respectively).

At 5 reads, targeted repair improved from 53/96 to 68/96, with 18 exact-only wins versus 3 regressions relative to multi-start (McNemar p = 0.00149). Against BBS, the targeted method was 68/96 versus 73/96 at 5 reads and the paired exact difference was not statistically significant (p = 0.424). At 10 reads, BBS remained stronger at 93/96 versus 84/96 (p = 0.0117).

The targeted method took about **17.1 s** for the 96-cluster 5-read evaluation and **41.2 s** at 10 reads, with about **24.4 MB peak RSS**. The complete one-shot GitHub Actions workflow, including dataset verification, pinned BBS build, all reconstruction baselines, five BBS repeats per coverage, acceptance checks and artifact upload, finished in roughly **5 minutes**, below the 10-minute budget.

The acceptance rule was satisfied because held-out 10-read exact recovery was strictly greater than the prior 80/96 result and 5-read recovery did not regress. The ideal 86/96 target was not reached, and BBS remains the stronger exact-recovery method on this physical dataset. The defensible conclusion is therefore that the targeted low-compute repair **materially improves OligoArk and narrows the external gap, but does not match or exceed BBS at 10 reads**.


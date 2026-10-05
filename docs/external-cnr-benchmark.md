# Fast external physical-read benchmark

This benchmark evaluates OligoArk v0.6 reconstruction methods on one public physical DNA-storage dataset and compares them with one recent external trace-reconstruction baseline that runs on a standard GitHub-hosted CPU runner.

## Dataset

The benchmark uses Microsoft's **Clustered Nanopore Reads (CNR)** dataset released with *Trellis BMA: coded trace reconstruction on IDS channels for DNA storage* (ISIT 2021). It contains 10,000 reference strands of length 110 and 269,709 clustered Oxford Nanopore MinION reads. `Centers.txt` is the explicit ground truth and the clusters in `Clusters.txt` are in the same order as those references.

The workflow pins dataset commit `6938f44796185902a08381943c2895782886c5c3` and verifies the Git blob checksums for both data files before running. The benchmark also records SHA-256 hashes of the checked-out files in its metadata.

The CNR maintainers added an important limitation in 2024: the generated centers contain unintended long-range dependencies, and some recovered clusters may therefore be malformed. OligoArk preserves this limitation in every benchmark artifact and does not generalize this dataset to all physical DNA-storage channels.

## External baseline

The external baseline is **Bidirectional Beam Search (BBS)** from Gu et al., *Efficient trace reconstruction in DNA storage systems using bidirectional beam search* (iScience, 2025). BBS directly supports the Microsoft CNR cluster format, has an open-source Rust implementation, and is practical on CPU-only CI. The workflow builds the official source at commit `3e4ab46871929819e4f3e34a831c57cac88bb456` with its locked dependencies and runs it with its default algorithm settings on one CPU thread.

HEDGES is not included in this same-read benchmark because HEDGES is an encoding plus decoding code: the CNR physical strands were not synthesized with HEDGES. Retrofitting HEDGES would change the encoded data rather than compare reconstruction on identical reads. The same restriction applies to OligoArk fountain/hybrid redundancy and the combined archive optimizer, so those are explicitly reported as not applicable rather than assigned invented scores.

## Fast deterministic design

To keep the workflow suitable for ordinary GitHub-hosted runners, it selects 96 clusters from those having at least 10 physical reads. Selection is deterministic: eligible cluster indexes are ranked by SHA-256 with seed `20261005`; reconstruction quality and reference content are never used for selection. Within each chosen cluster, reads are independently SHA-256 ranked before taking nested coverages of 1, 5 and 10 reads. Every reconstruction method therefore receives the exact same reads and reference strand at a given coverage.

The evaluated methods are:

- OligoArk direct: the first deterministically selected physical read;
- OligoArk medoid consensus;
- OligoArk graph/alignment using the v0.6 default graph settings and the largest graph component without reference-guided candidate selection;
- OligoArk iterative trace consensus with the frozen v0.6 three-round setting;
- external BBS from the pinned official implementation.

No method is calibrated on the evaluation references. The reference is used only after reconstruction to score accuracy.

## Metrics

For each method and coverage, the workflow reports exact-reference recovery with Wilson 95% confidence intervals, edit distance, normalized edit distance, CPU-runner wall time, per-cluster time, process peak RSS, and deterministic subset SHA-256 comparison. Pairwise exact-recovery differences against BBS use a two-sided exact McNemar test on the same clusters. Raw per-cluster CSV/JSON, summaries, paired comparisons, timing records, deterministic subset files and environment/provenance metadata are uploaded as a 90-day GitHub Actions artifact.

The raw-read profile additionally measures substitution, insertion and deletion counts against the known references on the selected 10-read subset. The full-dataset empty-cluster rate is reported as the observable strand/dropout analogue. These are descriptive measurements of this dataset, not universal sequencing error rates.

## Claim boundary

This is a **physical-read reference-reconstruction benchmark**, not an end-to-end OligoArk archive decode. CNR does not contain a file encoded with OligoArk framing, ECC or redundancy, so file-level SHA-256 recovery, OligoArk fountain/hybrid recovery and the combined archive optimizer are not applicable. The benchmark instead hashes the ordered reference subset and reconstructed subset to provide an exact aggregate integrity check.


## Executed v0.6 result

The fast held-out physical-read benchmark uses 96 deterministically selected clusters from 8,978 CNR clusters having at least 10 reads and evaluates nested coverages of 1, 5 and 10 reads. The selected 10-read subset contained 960 physical reads. Against the known references, the descriptive raw-read profile measured 2.23% substitutions, 1.65% insertions and 1.92% deletions per aligned reference base.

| Reads / strand | Direct | Medoid | OligoArk graph/alignment | OligoArk iterative trace | External BBS |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1/96 (1.0%) | 1/96 (1.0%) | 1/96 (1.0%) | 1/96 (1.0%) | 1/96 (1.0%) in the exact-head run |
| 5 | 1/96 (1.0%) | 6/96 (6.3%) | 33/96 (34.4%) | 44/96 (45.8%) | **73/96 (76.0%)** |
| 10 | 1/96 (1.0%) | 10/96 (10.4%) | 65/96 (67.7%) | 63/96 (65.6%) | **93/96 (96.9%)** |

At one read, exact-recovery differences versus BBS were not statistically significant. An earlier identical official-code BBS execution produced 2/96 instead of 1/96 because the upstream implementation does not define deterministic tie-breaking for equal-score candidates stored in randomized Rust `HashMap`s. The workflow therefore records five unmodified BBS repetitions per coverage so this implementation-level variability remains visible. At five reads, BBS exceeded graph/alignment by 40 paired successes versus one graph-only success (exact two-sided McNemar p = 1.96e-11) and exceeded iterative trace by 29 net paired successes (31 BBS-only versus two trace-only; p = 1.31e-7). At ten reads, BBS exceeded graph/alignment on 28 paired clusters with no graph-only wins (p = 7.45e-9) and exceeded iterative trace on 30 paired clusters with no trace-only wins (p = 1.86e-9).

Mean edit distance at five reads was 0.625 for BBS, 1.010 for iterative trace and 1.938 for graph/alignment. At ten reads it was 0.354 for BBS, 0.469 for iterative trace and 0.490 for graph/alignment. Therefore, on this single public physical CNR benchmark, the defensible result is that **OligoArk v0.6 reconstruction underperforms the external BBS baseline at useful multi-read coverage**. The result does not invalidate OligoArk's internal simulation gains; it shows that those gains do not establish external state-of-the-art performance.

The benchmark intentionally reports OligoArk fountain/hybrid redundancy, the combined archive optimizer and full-file SHA-256 recovery as not applicable because the CNR strands were not encoded with OligoArk. No result is inferred for those methods.

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

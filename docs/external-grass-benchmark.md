# Independent Grass physical-read benchmark

This benchmark tests the **unchanged merged OligoArk confidence-fusion reconstruction method** on an independent physical DNA-storage dataset derived from the Grass et al. experiment. The purpose is cross-dataset generalization, not dataset-specific tuning.

## Dataset and provenance

The benchmark uses the prepared binned Grass dataset published on Zenodo as record `14296588` (DOI `10.5281/zenodo.14296588`) and derived from the physical DNA-storage experiment reported in *Angewandte Chemie International Edition* (DOI `10.1002/anie.201411378`).

The pinned benchmark file is `Grass.txt`:

- published MD5: `b076770ad26e955ff60b94e6344e3a05`
- measured SHA-256: `e1e7257fe023590983c1af7506cdbc34b0cda1d2eb396490f54dc3552247b8aa`
- prepared benchmark license: CC BY
- parsed clusters: **4,982**
- clusters with at least 10 associated reads: **4,911**
- encoded reference length in the prepared file: **117 nt**
- physical experiment: CustomArray electrochemical oligo synthesis
- sequencing: Illumina MiSeq 2×150 TruSeq

Each prepared cluster stores the encoded reference strand followed by its associated physical reads. This makes reference-level reconstruction measurable without inventing a new read-clustering algorithm.

## Frozen OligoArk configuration

The benchmark uses the current merged confidence-fusion configuration unchanged:

`fusion-fast-t2-c8-r0-q025-g005`

Only the mechanically required target length changes from the CNR value of 110 to the Grass reference length of 117. No search, scoring, candidate-generation, or threshold parameter was selected on the Grass held-out references.

A deterministic 48-cluster diagnostic development split uses seed `20261011`. The final held-out split contains 96 disjoint clusters selected with seed `20261012`. Reads within each selected cluster are deterministically ranked with SHA-256 before taking nested 1-, 5-, and 10-read subsets.

The development diagnostic reached **48/48 exact recovery at both 5 and 10 reads**. No configuration change was made from this result.

## Held-out result

The authoritative successful validation is GitHub Actions workflow run `37511642177`, artifact `11434354623`.

| Reads / strand | **OligoArk confidence fusion** | Pinned BBS |
| ---: | ---: | ---: |
| 1 | 35/96 (36.5%) | 35/96 (36.5%) |
| 5 | **94/96 (97.9%)** | 90/96 (93.8%) |
| 10 | **96/96 (100%)** | 95/96 (99.0%) |

Wilson 95% confidence intervals for OligoArk were:

- 1 read: **27.5%–46.4%**
- 5 reads: **92.7%–99.4%**
- 10 reads: **96.2%–100%**

At 5 reads, OligoArk mean edit distance was **0.0208** versus **0.0625** for BBS. OligoArk had two remaining one-edit substitution failures.

At 10 reads, OligoArk reconstructed **all 96/96 held-out strands exactly**, with mean and median edit distance **0** and no one-, two-, or 3+-edit failures. BBS reconstructed 95/96 and had one one-edit substitution failure.

BBS was repeated five times per coverage. Exact recovery was stable at **35/96**, **90/96**, and **95/96** for 1, 5, and 10 reads respectively.

## Paired comparison

Using the first pinned BBS repeat on identical held-out clusters:

| Reads / strand | OligoArk-only exact | BBS-only exact | Exact McNemar p |
| ---: | ---: | ---: | ---: |
| 1 | 0 | 0 | 1.0 |
| 5 | 4 | 0 | 0.125 |
| 10 | 1 | 0 | 1.0 |

The numerical advantage for OligoArk at 5 and 10 reads is **not statistically significant** on this 96-cluster subset. The defensible conclusion is cross-dataset generalization, not superiority.

## Physical-read error profile

Across the 960 held-out reads used in the 10-read condition:

- substitution rate per aligned reference base: **0.491%**
- insertion rate: **0.048%**
- deletion rate: **0.576%**

This Illumina dataset has a different and substantially lower raw-error regime than the Microsoft CNR nanopore benchmark, providing an independent physical sequencing condition for the frozen method.

## Runtime and memory

On the GitHub-hosted CPU runner:

- OligoArk: **37.33 s** at 5 reads and **86.53 s** at 10 reads, about **24 MB peak RSS**
- BBS: about **0.16 s** at 5 reads and **0.25 s** at 10 reads, about **3.5–3.8 MB peak RSS**
- complete workflow: about **6 min 44 s**, within the 10-minute budget

BBS is therefore substantially faster despite the numerical exact-recovery difference on this subset.

## Claim boundary

This is a reference-level physical-read reconstruction benchmark, not end-to-end decoding of the original Grass archive through OligoArk framing/ECC. The result demonstrates that the **same confidence-fusion settings developed on CNR generalize strongly to an independent Illumina physical DNA-storage dataset**. It does not establish a general state-of-the-art claim.

# DNAformer Pilot Illumina physical-read benchmark

This benchmark evaluates the **unchanged merged OligoArk confidence-fusion reconstruction method** on the DNAformer Pilot Illumina physical DNA-storage dataset. It is an additional cross-dataset generalization test beyond Microsoft CNR, Grass, and LCRC; no reconstruction parameter was tuned on the DNAformer held-out references.

## Dataset and provenance

The benchmark uses `BinnedPilotIllumina.txt` from Zenodo record `13896773` (DOI `10.5281/zenodo.13896773`), associated with Bar-Lev et al., *Scalable and robust DNA-based storage via coding theory and deep learning* (Nature Machine Intelligence, DOI `10.1038/s42256-025-01003-z`).

The Zenodo record describes the Pilot Illumina data as a physical dataset synthesized by **Twist Bioscience** and sequenced with **Illumina MiSeq paired-end reads**, with paired reads stitched using **PEAR**. The binned format explicitly stores each encoded reference strand followed by its associated physical reads.

Pinned file provenance:

- file: `BinnedPilotIllumina.txt`
- published MD5: `02b16f63adf5ac85e302c151eaa8ad48`
- measured SHA-256: `c6a88a1a31b2ba59745d127118b309bca41f013da3fe52ec38a76e62915e9102`
- parsed clusters: **1,000**
- clusters with at least 10 reads: **1,000**
- encoded reference length: **140 nt**
- Zenodo states the dataset uses licensing similar to the associated code repository; the associated code repository is MIT licensed.

## Frozen OligoArk configuration

The benchmark uses the already-merged confidence-fusion configuration unchanged:

`fusion-fast-t2-c8-r0-q025-g005`

Only the known target strand length is inferred mechanically from the dataset. The search size, q-gram weight, score-gain threshold, trimming rule, anchor count, and rounds are unchanged from the CNR-developed configuration.

A deterministic 48-cluster diagnostic split uses seed `20261013`. The final held-out split contains 96 disjoint clusters selected with seed `20261014`. Reads within selected clusters are deterministically ranked with SHA-256 before taking nested 1-, 5-, and 10-read subsets.

The diagnostic split reached **48/48 exact recovery** at both 5 and 10 reads. No algorithm or parameter change was made from this diagnostic result.

## Held-out result

The authoritative validation is GitHub Actions workflow run `37566948402`, artifact `11459885755`, artifact digest `sha256:f6743a04ca49a7f3479432811305f842dad4a54d9e59428fb1707f1167409450`.

| Reads / strand | **OligoArk confidence fusion** | Pinned BBS |
| ---: | ---: | ---: |
| 1 | 83/96 (86.5%) | 83/96 (86.5%) |
| 5 | **96/96 (100%)** | **96/96 (100%)** |
| 10 | **96/96 (100%)** | **96/96 (100%)** |

For OligoArk, Wilson 95% confidence intervals were:

- 1 read: **78.2%–91.9%**
- 5 reads: **96.2%–100%**
- 10 reads: **96.2%–100%**

At 1 read, OligoArk mean edit distance was **0.1458**. At 5 and 10 reads, mean and median edit distance were **0**, with no one-, two-, or 3+-edit failures.

BBS was repeated five times at every coverage and was stable at **83/96, 96/96, and 96/96** for 1, 5, and 10 reads respectively.

## Paired comparison

The first pinned BBS repeat and OligoArk produced identical exact-recovery outcomes at every tested coverage:

| Reads / strand | OligoArk-only exact | BBS-only exact | Exact McNemar p |
| ---: | ---: | ---: | ---: |
| 1 | 0 | 0 | 1.0 |
| 5 | 0 | 0 | 1.0 |
| 10 | 0 | 0 | 1.0 |

There is therefore no exact-recovery difference between the methods on this held-out subset.

## Physical-read error profile

Across the 960 held-out reads used for the 10-read condition:

- substitution rate per aligned reference base: **0.0551%**
- insertion rate: **0.0060%**
- deletion rate: **0.0067%**

This is a substantially lower-error physical sequencing regime than CNR and provides a separate experiment and synthesis/readout combination from the earlier benchmark datasets.

## Runtime and memory

On the GitHub-hosted CPU runner:

- OligoArk: **34.02 s** at 5 reads and **73.76 s** at 10 reads, about **24 MB peak RSS**
- BBS: about **0.05 s** at 5 reads and **0.08 s** at 10 reads, about **3.6–3.7 MB peak RSS**
- complete workflow: about **4 min 55 s**, within the 10-minute budget

BBS remains substantially faster.

## Claim boundary

This is a reference-level physical-read reconstruction benchmark using already-binned read clusters, not end-to-end decoding of the original DNAformer archive. The result provides additional evidence that the **same frozen confidence-fusion configuration generalizes across multiple independent physical datasets**, but it does not establish a general state-of-the-art claim.

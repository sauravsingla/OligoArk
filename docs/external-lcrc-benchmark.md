# Third independent LCRC physical-read benchmark

This benchmark evaluates the current merged OligoArk **confidence-fusion** reconstruction method on a third independent physical DNA-storage experiment: the 2026 LCRC HFS-Pool-11.7K high-fidelity oligo pool from Zhang et al., *Science Advances* (DOI `10.1126/sciadv.aec1469`). The confidence-fusion search and scoring parameters are unchanged from the Microsoft CNR and Grass benchmarks.

## Dataset and provenance

The benchmark pins the public `dna-storage-lab/DNAStorage_LCRC` repository at commit:

`ce35bb2720c49ea6d1e6bc17903655b45f9c6c51`

The source experiment is also deposited in SRA under accession `PRJNA1371011`. The repository is MIT-licensed; the sequencing data are publicly available through the repository/SRA, while the repository does not state a separate sequencing-data license.

Files used:

- reference design: `LCRC_SmallScale/reference/DNA_oligoPool/oligoPool_11.7K.fa`
  - Git blob: `3676efd8535e0b5554408154234ed4e2dd1d18ab`
  - measured SHA-256: `316f3a957e37433764666c9c6e44295b931171ff7a85b693418877fcbe2561eb`
- real physical reads: `LCRC_SmallScale/fastq/HFS_Pool_11.7K_NGS_PEmerged.fastq.gz`
  - Git blob: `c9926ac7e348e0aba8906ccf6dfd02133ad6b433`
  - measured SHA-256: `86fc3a7269761abf0c0a7b0d9205f198119b995687618346bb1766640b5d1b28`

The experiment contains **11,745 designed oligos of 200 nt**. HFS-Pool-11.7K was synthesized by **Twist Bioscience** using the study's high-fidelity synthesis condition and sequenced with **150-nt paired-end NGS/Illumina sequencing**; the repository supplies the paired-end merged FASTQ used here.

## Read-to-reference association

Unlike CNR and the prepared Grass benchmark, the raw LCRC FASTQ is not already binned by reference. The benchmark therefore performs a deterministic association step before reconstruction.

Only the **interior** of each public 200-nt design is used to build a unique 13-mer anchor index, excluding the shared 20-nt primer regions at each end. Anchors are sampled every 8 nt. Each observed read and its reverse complement are scanned for unique anchors. A read is assigned only when:

- at least 3 unique anchor hits support the same reference;
- the top reference exceeds the second-best anchor count by at least 1; and
- a final bounded edit-distance sanity check to that reference is at most 15 edits.

This association step uses the design library solely to determine which physical reads belong to which strand. **The reference is never passed to OligoArk reconstruction or candidate scoring.** OligoArk receives only the associated observed reads plus the mechanically known 200-nt target length.

Association results:

- total FASTQ reads: **46,749**
- mapped reads: **46,392 (99.24%)**
- rejected/ambiguous reads: **357**
- populated reference clusters: **11,030**
- clusters with at least 10 mapped reads: **545**
- unique interior anchor k-mers: **222,332**
- median anchor support for accepted reads: **19**
- median assignment edit distance: **0**
- association runtime: **11.0 s**

This mapping approach is intentionally simple, deterministic, codec-agnostic, and independent of the authors' LCRC decoder.

## Frozen OligoArk configuration

The benchmark uses the same confidence-fusion configuration developed before the Grass/LCRC validations:

`fusion-fast-t2-c8-r0-q025-g005`

The only dataset-specific reconstruction input is the published strand length, **200 nt**. No scoring weight, candidate limit, q-gram width, confidence threshold, or search parameter was selected on LCRC references.

A deterministic 48-cluster diagnostic development split uses seed `20261013`. The untouched 96-cluster held-out split uses seed `20261014`, after excluding all development cluster IDs. Both splits are selected only from clusters with at least 10 accepted physical reads.

The diagnostic development split achieved **48/48 exact reconstruction at both 5 and 10 reads**. No OligoArk configuration change was made.

## Held-out result

The authoritative run is GitHub Actions workflow `37564250037`, artifact `11458730863`, artifact digest:

`sha256:078676b9947a9833ff2cee939e412c0cbc97689d0e3a378ae39ac650bcc13127`

| Reads / strand | **OligoArk confidence fusion** | Pinned BBS |
| ---: | ---: | ---: |
| 1 | 69/96 (71.9%) | 69/96 (71.9%) |
| 5 | **96/96 (100%)** | **96/96 (100%)** |
| 10 | **96/96 (100%)** | **96/96 (100%)** |

Wilson 95% confidence intervals for OligoArk:

- 1 read: **62.2%–79.9%**
- 5 reads: **96.2%–100%**
- 10 reads: **96.2%–100%**

At 5 and 10 reads, OligoArk had mean and median edit distance **0**, with **no one-, two-, or 3+-edit failures**.

BBS was repeated five times at each coverage. Exact recovery was stable at **96/96 for both 5 and 10 reads** in all repeats. At one read, exact recovery was also 69/96 in each repeat, although equal-score output sequences can differ between repetitions.

## Paired comparison

Using the first pinned BBS repeat on identical held-out clusters:

| Reads / strand | OligoArk-only exact | BBS-only exact | Exact McNemar p |
| ---: | ---: | ---: | ---: |
| 1 | 0 | 0 | 1.0 |
| 5 | 0 | 0 | 1.0 |
| 10 | 0 | 0 | 1.0 |

There is therefore no exact-recovery difference between the two methods on this held-out subset at any tested coverage.

## Physical-read error profile

Across the 960 held-out reads used in the 10-read condition:

- substitutions: **244**
- insertions: **13**
- deletions: **114**
- substitution rate per aligned reference base: **0.127%**
- insertion rate: **0.0068%**
- deletion rate: **0.059%**

This is a substantially lower-error physical regime than Microsoft CNR and is consistent with the high-fidelity HFS pool and Illumina/NGS readout. It provides independent laboratory, synthesis, indexing, strand-length, and protocol validation, but it is not a harder benchmark than CNR.

## Runtime and memory

On the GitHub-hosted CPU runner:

- read association: **11.0 s**
- OligoArk: **49.64 s** at 5 reads and **121.81 s** at 10 reads
- OligoArk peak RSS: about **25 MB**
- BBS: about **0.07 s** at 5 reads and **0.14 s** at 10 reads
- BBS peak RSS: about **4 MB**
- complete workflow: about **5 min 42 s**

The workflow therefore remains within the 10-minute CPU-only budget. BBS remains much faster.

## Limitations and claim boundary

The held-out subset is selected from the **545 references with at least 10 confidently associated reads**, not from all 11,745 designed oligos. The result therefore measures reconstruction quality conditional on sufficient read coverage and successful association; it is not an estimate of whole-pool retrieval yield.

The reference library is used in the association layer, as required to recover read-to-strand membership from an unbinned FASTQ. It is not used by the reconstruction algorithm itself.

This benchmark is reference-level physical-read reconstruction, **not end-to-end decoding through the authors' LCRC/LDPC archive codec** and not an OligoArk archive decode.

Together with CNR and Grass, the result provides evidence that the **same frozen confidence-fusion settings generalize across three independent physical DNA-storage experiments** with different laboratories, designs, strand lengths, synthesis protocols, sequencing conditions, and error profiles. It does **not** establish a general state-of-the-art claim.

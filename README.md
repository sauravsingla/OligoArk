# OligoArk 🧬

**DNA archival storage for large files with exact recovery.**

OligoArk is a research project for storing digital data as DNA-like sequences and recovering
the original file exactly.

The goal is simple: **store large files, use little RAM, survive damaged or missing DNA
strands, and recover the exact original data.**

OligoArk checks success using **SHA-256**. If the recovered file does not match exactly, the
test is counted as a failure.

## Why DNA storage?

DNA storage is **not mainly about compression**. A 1 GB file does not automatically become a
few MB.

The advantage of DNA is physical storage: very large amounts of data could be stored in a very
small amount of DNA, kept for a long time, and need little or no power while sitting in
storage.

OligoArk focuses on making this process scalable and reliable.

## What works today?

OligoArk has four strong results:

- exact recovery of a **1 GiB mixed dataset** with about **44 MiB peak RAM**;
- exact recovery with **1% and 5% strand loss**;
- much better recovery from **insert/delete and mixed DNA errors** on 248-nt strands;
- about **51% better storage density** on the improved 248-nt profile.

## 1 GiB storage result

OligoArk successfully recovered a **1 GiB dataset** containing different types of data.

| Test | Result | Peak RAM |
| --- | ---: | ---: |
| Clean | ✅ PASS | 43.85 MiB |
| 1% strand loss | ✅ PASS | 44.08 MiB |
| 5% strand loss | ✅ PASS | 43.63 MiB |

At 1% loss, OligoArk recovered **45,338 / 45,338** lost strands.

At 5% loss, it recovered **226,705 / 226,705** lost strands.

The 1 GiB test used **1.646 bits per nucleotide** and 12.5% redundancy.

Exact recovery also passed at:

**1 KiB → 64 KiB → 1 MiB → 10 MiB → 100 MiB → 1 GiB**

## Better recovery on 248-nt strands

The improved compact 248-nt profile now handles insertion/deletion errors much better.

| Protection | Clean | 1% loss | 5% loss | Substitution | Insert/Delete | Mixed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| None | ✅ | ❌ | ❌ | ❌ | ✅ | ❌ |
| XOR | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Fountain | ✅ | ❌ | ❌ | ❌ | ✅ | ❌ |
| Hybrid | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |

The strongest result is that **XOR and Hybrid pass every tested condition**, including
insert/delete errors and mixed noise.

In the indel test, XOR repaired **71 damaged records** and Hybrid repaired **84 damaged
records**, while still recovering the exact original file.

## Better storage density

The improved 248-nt profile stores more useful data in the same strand length.

| Method | Before | Now | Improvement |
| --- | ---: | ---: | ---: |
| XOR | 1.003 bits/nt | **1.519 bits/nt** | about **51%** |
| Hybrid | 0.820 bits/nt | **1.242 bits/nt** | about **51%** |

Higher bits/nt means more useful data can be stored in the same number of DNA bases.

## Comparison with other DNA-storage codecs

The largest fully completed fair comparison is currently **1 MiB**.

| Method | Density | Clean | 1% loss | 5% loss | Substitution |
| --- | ---: | ---: | ---: | ---: | ---: |
| OligoArk Fountain | 0.463 bits/nt | 5/5 | 0/5 | 0/5 | **5/5** |
| DNA Fountain | **1.347 bits/nt** | 5/5 | **5/5** | **5/5** | 0/5 |
| Goldman-style + XOR | 0.515 bits/nt | 5/5 | 0/5 | 0/5 | 0/5 |

In simple terms, **DNA Fountain is still better in this test for density and strand-loss
recovery**, while OligoArk performed better in the tested substitution condition.

OligoArk also attempted **10 MiB and 100 MiB** matched comparisons with 10 trials per
condition, but the workers exceeded the 1,200-second limit. These are reported as
**timeouts, not successful results**.

## Other physical-read datasets

These tests use published DNA sequencing datasets created by other projects. They test
OligoArk's ability to reconstruct reference strands.

| Dataset | Reads/strand | OligoArk | Pinned BBS |
| --- | ---: | ---: | ---: |
| Microsoft CNR | 5 | **73/96** | 72–74/96 |
| Microsoft CNR | 10 | **93/96** | 93/96 |
| Grass et al. | 5 | **94/96** | 90/96 |
| Grass et al. | 10 | **96/96** | 95/96 |
| LCRC HFS-11.7K | 5 | **96/96** | 96/96 |
| LCRC HFS-11.7K | 10 | **96/96** | 96/96 |
| DNAformer Pilot | 1 | **83/96** | 83/96 |
| DNAformer Pilot | 5 | **96/96** | 96/96 |
| DNAformer Pilot | 10 | **96/96** | 96/96 |

These results are useful reconstruction evidence, but they are **not a physical OligoArk
storage experiment** because those DNA strands were created by other systems.

## Wet-lab status

OligoArk can now prepare a real lab experiment. The repo includes:

- synthesis-ready FASTA and CSV files;
- experiment metadata and manifest format;
- sequencing input requirements;
- reconstruction scripts;
- read-depth plans;
- exact SHA-256 verification.

**Current status: ready for a wet-lab experiment, but not yet physically tested end to end.**

A full physical proof still needs:

**OligoArk encode → DNA synthesis → storage → DNA sequencing → reconstruction → exact SHA-256
match**

## Current limitations

The main things still to improve are the large 10–100 MiB matched codec comparison and a true
wet-lab OligoArk storage experiment.

The software results are strong, but OligoArk should not yet be described as better than every
DNA-storage codec or as a fully proven physical DNA-storage system.

More details: [scalable storage](docs/scalable-storage.md) ·
[codec comparison](docs/dna-fountain-baseline.md) ·
[wet-lab protocol](docs/wet-lab-validation.md)

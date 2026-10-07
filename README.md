# OligoArk 🧬

**DNA archival storage for large files with exact recovery.**

## Problem statement

DNA can store data for a very long time, but a practical DNA-storage system must do more than
encode small files. It should store large and mixed types of data with low memory use, recover
the original file exactly when strands are lost or damaged, work with realistic short DNA
strands, and report density, speed, memory and redundancy clearly.

OligoArk uses **SHA-256 exact recovery** as the final success criterion.

## Why DNA storage?

DNA storage is not mainly about making a digital file smaller. Its advantage is the possibility
of storing very large amounts of information in a tiny amount of physical material for long
periods with little or no power while the archive is at rest. OligoArk focuses on making that
storage process scalable, memory-efficient and recoverable.

## Benchmark results

### 1 GiB scalable storage

OligoArk stored and exactly recovered a **1 GiB heterogeneous dataset** while using about
**44 MiB peak RAM**.

| Test | Result | Lost strands recovered | Density | Encode | Decode | Peak RAM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Clean | ✅ PASS | 0 | **1.646 bits/nt** | 9.78 MiB/s | 23.08 MiB/s | 43.85 MiB |
| 1% strand loss | ✅ PASS | **45,338 / 45,338** | **1.646 bits/nt** | 9.77 MiB/s | 17.91 MiB/s | 44.08 MiB |
| 5% strand loss | ✅ PASS | **226,705 / 226,705** | **1.646 bits/nt** | 9.95 MiB/s | 12.12 MiB/s | 43.63 MiB |

**Redundancy:** 12.5% · **Archive overhead:** 1.234× · **Total strands:** 5,096,877

Exact SHA-256 recovery also passed at:

**1 KiB → 64 KiB → 1 MiB → 10 MiB → 100 MiB → 1 GiB**

### Improved 248-nt short-strand recovery

The new compact 248-nt profile adds bounded single-indel realignment and lower framing
overhead. The table below uses a 64 KiB heterogeneous payload.

| Protection | Clean | 1% loss | 5% loss | Substitution | Indel | Mixed | Density |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| None | ✅ | ❌ | ❌ | ❌ | ✅ | ❌ | **1.710 bits/nt** |
| XOR | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | **1.519 bits/nt** |
| Fountain | ✅ | ❌ | ❌ | ❌ | ✅ | ❌ | **1.367 bits/nt** |
| Hybrid | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | **1.242 bits/nt** |

The strongest result is **XOR and Hybrid passing all six tested conditions**, including indel
and mixed noise. In the indel test, XOR repaired **71** damaged records and Hybrid repaired
**84** while still reaching the exact source SHA-256.

Compared with the previous 248-nt profile, density improved by about **51%** for both XOR
(**1.003 → 1.519 bits/nt**) and Hybrid (**0.820 → 1.242 bits/nt**).

The compact profile trades some inner Reed-Solomon protection for density, so unprotected
None/Fountain modes no longer pass the tested substitution condition. Negative results are
kept rather than hidden.

### Matched codec comparison

The established common matched comparison remains **1 MiB**, using the same 152-nt limit,
25% nominal redundancy, fault rates, seeds and exact SHA-256 gate.

| Method | Density | Clean | 1% loss | 5% loss | Substitution |
| --- | ---: | ---: | ---: | ---: | ---: |
| OligoArk Fountain (previous matched run) | 0.463 bits/nt | 5/5 | 0/5 | 0/5 | **5/5** |
| DNA Fountain | **1.347 bits/nt** | 5/5 | **5/5** | **5/5** | 0/5 |
| Goldman-style + XOR | 0.515 bits/nt | 5/5 | 0/5 | 0/5 | 0/5 |

PR #39 also attempted **10 MiB and 100 MiB** matched runs with **10 trials per condition**.
All three methods exceeded the current **1,200-second worker limit**, so those larger sizes
are reported as **timeouts, not successful comparisons**. The largest statistically matched
common completed size therefore remains 1 MiB.

### Other dataset benchmarks

These tests use published physical-read datasets and measure **reference-strand reconstruction**.

| Dataset | Reads/strand | OligoArk | Pinned BBS |
| --- | ---: | ---: | ---: |
| Microsoft CNR (Nanopore) | 5 | **73/96 (76.0%)** | 72–74/96 |
| Microsoft CNR (Nanopore) | 10 | **93/96 (96.9%)** | **93/96 (96.9%)** |
| Grass et al. (Illumina) | 5 | **94/96 (97.9%)** | 90/96 (93.8%) |
| Grass et al. (Illumina) | 10 | **96/96 (100%)** | 95/96 (99.0%) |
| LCRC HFS-11.7K | 5 | **96/96 (100%)** | **96/96 (100%)** |
| LCRC HFS-11.7K | 10 | **96/96 (100%)** | **96/96 (100%)** |
| DNAformer Pilot | 1 | **83/96 (86.5%)** | **83/96 (86.5%)** |
| DNAformer Pilot | 5 | **96/96 (100%)** | **96/96 (100%)** |
| DNAformer Pilot | 10 | **96/96 (100%)** | **96/96 (100%)** |

### Wet-lab readiness

OligoArk now includes a lab handoff workflow that generates synthesis-ready FASTA/CSV files,
experiment metadata, a manifest schema, sequencing input requirements, reconstruction scripts,
read-depth plans and an exact SHA-256 verification gate.

**Status: prepared, not physically executed.**

> **Claim boundary:** The 1 GiB result and 248-nt fault matrix are software benchmarks.
> CNR/Grass/LCRC/DNAformer are reference-strand reconstruction benchmarks. OligoArk has not
> yet demonstrated an end-to-end synthesized → sequenced → SHA-256-verified physical archive.

More details: [scalable storage](docs/scalable-storage.md) ·
[codec comparison](docs/dna-fountain-baseline.md) ·
[wet-lab protocol](docs/wet-lab-validation.md)

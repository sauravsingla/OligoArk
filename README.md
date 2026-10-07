# OligoArk 🧬

**Scalable DNA archival-storage and reconstruction research.**

## Problem statement

DNA archival storage needs to scale to large heterogeneous files, recover the original bytes
exactly under strand loss/noise, and do so with measurable density, redundancy, throughput and
memory under realistic strand constraints. OligoArk uses bounded-memory streaming archives,
configurable 150–250 nt physical-design profiles, redundancy/error simulation, reconstruction,
and **SHA-256 exact recovery as the final success criterion**.

## Storage benchmark

### 1 GiB bounded-memory archive

Validated on a deterministic 1 GiB heterogeneous payload using the `scale-1024` systems
profile with XOR redundancy.

| Condition | SHA-256 | Lost / recovered | Density | Encode | Decode | Peak RSS | Redundancy | Archive overhead |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Clean | ✅ PASS | 0 / 0 | **1.645833 bits/nt** | 9.78 MiB/s | 23.08 MiB/s | 43.85 MiB | 12.5% | 1.234× |
| 1% controlled dropout | ✅ PASS | **45,338 / 45,338** | **1.645833** | 9.77 MiB/s | 17.91 MiB/s | 44.08 MiB | 12.5% | 1.234× |
| 5% controlled dropout | ✅ PASS | **226,705 / 226,705** | **1.645833** | 9.95 MiB/s | 12.12 MiB/s | 43.63 MiB | 12.5% | 1.234× |

**1 GiB archive:** 4,530,557 data strands + 566,320 parity strands = **5,096,877 strands**,
**5,219,201,308 encoded nt**, source SHA-256
`cc6286341b8650694bdb4f565a75372a981fde3935d1577eb415031d9fb12b5e`.

**Scaling:** 1 KiB → 64 KiB → 1 MiB → 10 MiB → 100 MiB → **1 GiB** all passed exact
SHA-256 recovery. Peak-memory log-log slope = **0.0434 (bounded)**; runtime slope =
**0.7935 (linear-or-better)**.

### RS-enabled 248-nt fault matrix

64 KiB heterogeneous payload, `oligoark-248`, 8 RS symbols.

| Scheme | Clean | 1% dropout | 5% dropout | Substitution | Indel | Mixed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| None | ✅ | ❌ | ❌ | ✅ | ❌ | ❌ |
| XOR | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| Fountain | ✅ | ❌ | ❌ | ✅ | ❌ | ❌ |
| Hybrid | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |

The current remaining robustness bottleneck is **insertion/deletion synchronization**; the
tested RS-enabled profiles recover substitutions but not the tested indel or mixed regimes.

### Matched 152-nt codec comparison

1 MiB payload, 25% nominal redundancy, identical fault rates/seeds, five trials per condition,
and the same SHA-256 exact-recovery gate.

| Method | Density | Clean | 1% dropout | 5% dropout | Substitution | Indel | Mixed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **OligoArk Fountain** | 0.463 bits/nt | 5/5 | 0/5 | 0/5 | **5/5** | 0/5 | 0/5 |
| **DNA Fountain clean-room** | **1.347 bits/nt** | 5/5 | **5/5** | **5/5** | 0/5 | 0/5 | 0/5 |
| **Goldman-style rotating + XOR** | 0.515 bits/nt | 5/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 |

### External physical-read reconstruction

| Dataset | Reads/strand | **OligoArk** | Pinned BBS |
| --- | ---: | ---: | ---: |
| Microsoft CNR (Nanopore) | 5 | **73/96 (76.0%)** | 72–74/96 |
| Microsoft CNR (Nanopore) | 10 | **93/96 (96.9%)** | **93/96 (96.9%)** |
| Grass et al. (Illumina) | 5 | **94/96 (97.9%)** | 90/96 (93.8%) |
| Grass et al. (Illumina) | 10 | **96/96 (100%)** | 95/96 (99.0%) |
| LCRC HFS-11.7K | 5 | **96/96 (100%)** | **96/96 (100%)** |
| DNAformer Pilot | 5 | **96/96 (100%)** | **96/96 (100%)** |

> **Claim boundary:** The 1 GiB result is software archive / controlled-channel evidence using
> the `scale-1024` systems profile. The 248-nt RS matrix is realistic-strand software evidence.
> CNR/Grass/LCRC/DNAformer are reference-strand reconstruction benchmarks. None of these is an
> end-to-end wet-lab OligoArk archive claim.

Details: [scalable storage](docs/scalable-storage.md) · [matched codec baselines](docs/dna-fountain-baseline.md)

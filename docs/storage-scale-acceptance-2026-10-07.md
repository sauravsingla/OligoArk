# Scalable-storage acceptance evidence — 2026-10-07

This report records the validation performed for PR #37 (`feature/scalable-dna-storage`).
It separates **software archive scaling**, **controlled software-channel faults**, and the
repository's existing **external physical-read reconstruction** evidence. It is not a wet-lab
OligoArk archive claim.

## Release milestone

**ACHIEVED for the software archive / controlled-loss evidence class.**

> OligoArk successfully archives and SHA-256 recovers a 100 MB heterogeneous dataset under
> clean and controlled-loss conditions, with measured density, throughput, peak memory and
> redundancy overhead.

The clean, 1% dropout and 5% dropout runs used the same 104,857,600-byte heterogeneous
payload with SHA-256
`6684d4fa7596e35f56768319d34ab5b229fad8cde440d3fe3b429f515a83dc59`.

| 100 MiB XOR case | SHA-256 | Dropped / recovered | bits/nt | Encode MiB/s | Decode MiB/s | Peak RSS MiB | Strand redundancy | Compact archive overhead |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| clean | PASS | 0 / 0 | 1.645832 | 11.78 | 43.88 | 108.42 | 12.5001% | 1.234185x |
| 1% controlled dropout | PASS | 4,421 / 4,421 | 1.645832 | 12.00 | 24.99 | 108.42 | 12.5001% | 1.234185x |
| 5% controlled dropout | PASS | 22,109 / 22,109 | 1.645832 | 11.66 | 16.25 | 108.04 | 12.5001% | 1.234185x |

The archive contains 442,438 data strands + 55,305 XOR parity strands = 497,743 software
strands and 509,688,008 encoded nucleotides. The compact binary container is 129,413,643
bytes for the 100 MiB source.

## Scaling

Clean XOR recovery passed SHA-256 at 1 KiB, 64 KiB, 1 MiB, 10 MiB and 100 MiB. Across those
isolated runs the log-log peak-memory slope was **0.01367**, classified as **bounded**, and
the runtime slope was **0.76661**, classified as **linear-or-better**. In the validation
environment Python itself started around 92 MiB RSS; the 100 MiB workload peaked around
108 MiB, so the process did not buffer a second 100 MiB payload in memory.

The earlier GitHub Actions scalable-storage job (run `37572396584`) also passed and measured
a bounded-memory slope of -0.01183 and linear-or-better runtime slope of 0.80368 on its CI
profile.

## Physical-strand profiles

A 1 KiB sanity archive was encoded and SHA-256 recovered with each physical-length profile.
All sampled strands passed the configured GC/homopolymer constraints.

| Profile | Full strand | Payload / strand | Sampled constraint pass | logical bits/nt | SHA-256 |
| --- | ---: | ---: | ---: | ---: | --- |
| oligoark-152 | 152 nt | 11 B | 100% | 0.509706 | PASS |
| oligoark-200 | 200 nt | 23 B | 100% | 0.806617 | PASS |
| oligoark-248 | 248 nt | 35 B | 100% | 0.983670 | PASS |

The separate 1024-nt software-scale profile intentionally disables biochemical constraints
and inner RS to measure container/streaming scalability. Its sampled physical-constraint pass
rate is therefore not a physical-storage claim.

## 1 MiB redundancy / fault matrix

The requested clean, 1%, 5%, substitution, indel and mixed-noise conditions were exercised
against none/XOR/fountain/hybrid. Results are intentionally retained even when negative.

| Scheme | clean | dropout 1% | dropout 5% | substitution | indel | mixed |
| --- | --- | --- | --- | --- | --- | --- |
| none | PASS | FAIL | FAIL | FAIL | FAIL | FAIL |
| XOR | PASS | **PASS** | **PASS** | FAIL | FAIL | FAIL |
| fountain (25%) | PASS | FAIL | FAIL | FAIL | FAIL | FAIL |
| hybrid | PASS | **PASS** | **PASS** | FAIL | FAIL | FAIL |

For XOR at 1 MiB, 44/44 controlled 1% erasures and 224/224 controlled 5% erasures were
recovered. The fast scale profile uses `rs_nsym=0`; therefore substitution/indel/mixed
failures are expected evidence that erasure redundancy alone is not sufficient for nucleotide
errors. Those regimes should use an RS-enabled physical profile and/or multi-read
reconstruction rather than weakening the SHA-256 gate.

Pure fountain recovery was changed to rescan fountain records from disk instead of retaining
all unresolved payloads in RAM. A regression test exercises this bounded-memory path, and a
second assertion confirms rescans reapply the identical channel corruption rather than
accidentally using pristine records.

## DNA Fountain comparison

The passing GitHub Actions scalable-storage artifact contains a matched 152-nt, 25%-nominal-
redundancy smoke comparison (1 KiB, 3 trials/condition):

| Method | clean | 5% dropout | logical bits/nt | strands |
| --- | ---: | ---: | ---: | ---: |
| OligoArk fountain | 3/3 | 0/3 | 0.457756 | 118 |
| clean-room DNA Fountain-style baseline | 3/3 | 3/3 | **1.347368** | 40 |

This is a small smoke experiment with wide confidence intervals, not a claim that DNA
Fountain universally dominates OligoArk. It does show a real current weakness in OligoArk's
152-nt fountain framing/redundancy and is preserved as a negative result.

## External physical-read comparison

The repository's existing held-out physical benchmarks remain the appropriate external
comparison because the public CNR/Grass/LCRC/DNAformer reads were not DNA-Fountain encoded.
The frozen OligoArk confidence-fusion settings are compared against pinned Bidirectional Beam
Search (BBS): CNR 73/96 vs 72–74/96 at 5 reads and 93/96 vs 93/96 at 10 reads; Grass 94/96
vs 90/96 at 5 reads and 96/96 vs 95/96 at 10 reads; LCRC and DNAformer reach 96/96 for both
methods at 5 and 10 reads. These are **reference-strand reconstruction** results, not
end-to-end OligoArk physical archive recovery.

## Fixes discovered during acceptance validation

1. ZIP members in the heterogeneous fixture now use fixed timestamps, so independent workers
   operate on the exact same deterministic payload.
2. Fountain recovery no longer retains unresolved fountain frames in memory; it rescans the
   compact on-disk archive.
3. Fountain rescans reapply the same deterministic channel faults, preventing clean-data
   leakage.
4. The fountain rescan offset is anchored by record ordinal even if the first fountain record
   is corrupted.
5. Regression tests cover bounded-memory fountain recovery and verify that noisy rescans do
   not see pristine data.

## Environment and claim boundary

Acceptance execution environment: Python 3.13.5, Linux 6.18.44 x86_64, glibc 2.41.
Throughput and RSS are environment-specific software measurements. No synthesis cost,
sequencing accuracy, biochemical durability or wet-lab end-to-end archive performance is
inferred from these results.

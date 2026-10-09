# OligoArk 🧬

**An open-source research framework for DNA data-storage coding, noisy-read reconstruction, and reproducible evaluation.**

OligoArk encodes digital files as DNA-like sequences, models strand loss and sequence errors, reconstructs data, and verifies recovered archives against their original SHA-256 digests. It is **simulation-first**: reference reconstruction from published physical sequencing reads has been evaluated, but an OligoArk-encoded archive has **not** undergone end-to-end synthesis and sequencing.

## Research contributions

- **Confidence-fusion read reconstruction:** multi-start, bidirectional trace consensus followed by bounded edits at uncertain positions, ranked against observed reads. The [reconstruction implementation](src/oligoark/reconstruct.py) and [independent physical-read benchmarks](docs/external-cnr-benchmark.md) document this approach.
- **Measured adaptive codec selection:** deterministic searches across coding, redundancy, and reconstruction configurations; disjoint calibration/evaluation seeds; inspectable trade-offs among recovery, overhead, runtime, and optional user-supplied lifecycle costs. See [research methods and ablations](docs/research.md).
- **Reproducible scaling and comparisons:** 152-nt constrained codec tests, exact-recovery trial accounting, explicit negative results, and a separate bounded-memory streaming-archive evaluation. See [benchmark methods](docs/benchmarking.md) and [scalable storage](docs/scalable-storage.md).

**Novelty boundary:** OligoArk combines and evaluates established ideas; fountain coding, Reed–Solomon protection, graph reconstruction, and sequence-constraint screening are not claimed as inventions. The confidence-fusion combination is a research contribution candidate, **not an established state-of-the-art or first-of-its-kind claim**.

## Quick start

Requires Python 3.10+ and Git:

```bash
git clone https://github.com/sauravsingla/OligoArk.git
cd OligoArk
python -m pip install -e .
python examples/quick_start.py
```

Expected output includes `PASS: original data recovered exactly`, `Bytes: 61`, a 64-character SHA-256 digest, and a note that the channel is simulated. The example uses deterministic seed 42; no lab equipment is required. [Source](examples/quick_start.py) · [End-to-end example](examples/end_to_end.py).

## Research results

**The 10 MiB codec comparison, partial 100 MiB research, 1 GiB streaming-storage test, and external physical-read reconstruction benchmarks are distinct experiments.** They use different configurations and must not be conflated. Exact archive recovery requires the original payload's matching SHA-256 digest; reference-strand reconstruction is assessed separately.

### 10 MiB matched-codec comparison — completed

One deterministic 10 MiB payload, **152-nt maximum strands**, approximately **25% redundancy**, identical simulated channel conditions and trial seeds, **10 trials per condition**:

| Condition | OligoArk compact-v3 hybrid | DNA Fountain-style baseline | Goldman-style + XOR |
| --- | ---: | ---: | ---: |
| Clean | **10/10** | **10/10** | **10/10** |
| 1% strand loss | **10/10** | **10/10** | 0/10 |
| 5% strand loss | **10/10** | **10/10** | 0/10 |
| Substitutions | **9/10** | 0/10 | 0/10 |
| Insertions/deletions | 0/10 (10 timeouts) | 0/10 | 0/10 |
| Mixed errors | **9/10** | 0/10 | 0/10 |

The fourth evaluated method, **OligoArk compact hybrid**, recovered 10/10 clean trials and 0/10 in each noisy condition. The 10 MiB compact-v3 indel-only trials all exceeded the **75-second per-trial deadline**: these are failures under the benchmark, not recoveries.

| 10 MiB measure | OligoArk compact-v3 | DNA Fountain-style baseline | Goldman + XOR |
| --- | ---: | ---: | ---: |
| Logical density (bits/nt) | 1.263 | 1.347 | 0.515 |
| Encoding time | 48.89 s | 271.99 s | Not reported here |
| Peak RAM | 436.33 MiB | 1.89 GiB | Not reported here |

OligoArk compact-v3's measured strand redundancy was **25.0007%**. Times are GitHub-runner observations, **not controlled proof of an algorithmic speedup**. These baselines are independently implemented approximations, not bit-compatible reproductions of the original published encoders. [Codec baseline and provenance](docs/dna-fountain-baseline.md).

### 100 MiB matched-codec research — partial validation

[Workflow run 37766665792](https://github.com/sauravsingla/OligoArk/actions/runs/37766665792) on [commit `d7b0f89`](https://github.com/sauravsingla/OligoArk/commit/d7b0f89a3e7521f02381a872be92b37342607f71) requested 10 trials in each of six matching conditions. **Three methods completed their trial outcomes; DNA Fountain did not.** This is **not** a complete four-method comparison.

| Condition | OligoArk compact-v3 | OligoArk compact hybrid | Goldman + XOR | DNA Fountain |
| --- | ---: | ---: | ---: | ---: |
| Clean | **10/10** | **10/10** | 0/10 | N/V |
| 1% strand loss | **10/10** | 0/10 | 0/10 | N/V |
| 5% strand loss | **1/10** | 0/10 | 0/10 | N/V |
| Substitutions | 0/10 | 0/10 | 0/10 | N/V |
| Insertions/deletions | 0/10 | 0/10 | 0/10 | N/V |
| Mixed errors | 0/10 | 0/10 | 0/10 | N/V |
| **Exact recoveries / 60 requested** | **21/60** | **10/60** | **0/60** | **N/V** |

**N/V = no verified trial outcome, not 0/10 attempted recoveries.** DNA Fountain exceeded its **19,800-second (5.5-hour) worker budget** and left 60 explicitly unverified placeholders. All **60 Goldman trials timed out** (including clean); compact-v3 had 10 indel and 10 mixed-error timeouts, and compact hybrid had 30 substitution/indel/mixed timeouts. Timed-out trials are recorded negative results within the specified deadline; they do not establish intrinsic codec impossibility.

| 100 MiB measure | Compact-v3 | Compact hybrid | Goldman + XOR |
| --- | ---: | ---: | ---: |
| Encode time | 602.58 s | 871.49 s | 594.03 s |
| Peak RAM | 3,810.70 MiB | 3,130.94 MiB | 3,213.53 MiB |
| Density (bits/nt) | 1.263157 | 1.221052 | 0.515436 |
| Maximum strand | 152 nt | 152 nt | 149 nt |
| Measured strand redundancy | 25.0001% | 25.00% | 25.00% |

The experiment used **104,857,600 source bytes**, a **152-nt ceiling**, **25% nominal redundancy**, seeds **20260000–20260009**, and a **300-second per-trial deadline**. Runtime/memory figures are single-run measurements; DNA Fountain has no completed measurements at this size. [Per-method artifacts and failure logs](https://github.com/sauravsingla/OligoArk/actions/runs/37766665792).

### Separate 1 GiB streaming-storage experiment — completed

The bounded-memory **software archive** test achieved exact SHA-256 recovery using **1.646 bits/nt** and **12.5% redundancy**:

| Controlled channel | SHA-256 | Peak RAM |
| --- | --- | ---: |
| Clean | PASS | 43.85 MiB |
| 1% data-strand loss | PASS | 44.08 MiB |
| 5% data-strand loss | PASS | 43.63 MiB |

These are **streaming archive / controlled-erasure results**, not the same 152-nt four-codec experiment. [Scale acceptance evidence](docs/storage-scale-acceptance-2026-10-07.md).

### Published physical-read reconstruction — reference level only

OligoArk's **frozen confidence-fusion** configuration was evaluated on held-out clusters from independent published sequencing datasets, against a pinned **Bidirectional Beam Search (BBS)** implementation:

| Dataset | 5 reads: OligoArk / BBS | 10 reads: OligoArk / BBS |
| --- | ---: | ---: |
| Microsoft CNR | 73/96 / 72–74/96 | **93/96 / 93/96** |
| Grass et al. | **94/96 / 90/96** | **96/96 / 95/96** |
| LCRC HFS-11.7K | **96/96 / 96/96** | **96/96 / 96/96** |
| DNAformer Pilot | **96/96 / 96/96** | **96/96 / 96/96** |

This shows cross-dataset generalisation on selected **reference-level** read clusters. The small exact-recovery differences are **not statistically significant** on these subsets, and **BBS is substantially faster**. Published reads were **not encoded with OligoArk's archive format**, so these are not end-to-end OligoArk storage recoveries.

[Microsoft CNR](docs/external-cnr-benchmark.md) · [Grass](docs/external-grass-benchmark.md) · [LCRC](docs/external-lcrc-benchmark.md) · [DNAformer](docs/external-dnaformer-pilot-benchmark.md).

## Limitations and physical validation

- **152-nt indel recovery at 10 MiB remains 0/10** under the specified timeout; substitution and mixed errors also limit the 100 MiB results.
- The **100 MiB four-method comparison is incomplete** because the DNA Fountain worker timed out before producing verified trial results.
- The repository has synthesis-ready FASTA/CSV and a [wet-lab validation protocol](docs/wet-lab-validation.md), but **physical OligoArk archive status remains: prepared, not physically executed**. A valid claim requires OligoArk encoding → DNA synthesis → sequencing → archive reconstruction → exact SHA-256 match.

## Further reading

[Research scope and prior work](docs/research.md) · [Architecture](docs/architecture.md) · [Benchmark methodology](docs/benchmarking.md) · [Configuration](docs/configuration.md) · [Contributing](CONTRIBUTING.md).

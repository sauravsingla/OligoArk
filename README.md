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

## Benchmark highlight

**Featured: 10 MiB matched-codec comparison — completed.** One deterministic 10 MiB payload, **152-nt maximum strands**, approximately **25% redundancy**, the same simulated channel conditions and trial seeds, and **10 trials per condition**.

| Simulated condition | OligoArk compact-v3 | DNA Fountain-style | Goldman-style + XOR |
| --- | ---: | ---: | ---: |
| Clean | **10/10** | **10/10** | **10/10** |
| 1% strand loss | **10/10** | **10/10** | 0/10 |
| 5% strand loss | **10/10** | **10/10** | 0/10 |
| Substitutions | **9/10** | 0/10 | 0/10 |
| Insertions/deletions | **0/10** (10 timeouts) | 0/10 | 0/10 |
| Mixed errors | **9/10** | 0/10 | 0/10 |

**Measured compact-v3 profile:** 1.263 logical bits/nt, 48.89 s encoding, and 436.33 MiB peak RAM in this run. These are **simulated, workload-specific results**, not physical DNA storage validation or proof of a general speedup. The baselines are independently implemented approximations, not bit-compatible reproductions of the original published encoders. All compact-v3 indel trials exceeded the **75-second per-trial deadline**. The additional compact-hybrid method recovered 10/10 clean trials but 0/10 in each noisy condition.

### All benchmark reports and evidence

- **10 MiB matched-codec comparison:** [Baseline implementation, methodology, and provenance](docs/dna-fountain-baseline.md) · [General benchmark methodology](docs/benchmarking.md).
- **100 MiB matched-codec research (partial):** [Workflow, per-method artifacts, and failure logs](https://github.com/sauravsingla/OligoArk/actions/runs/37766665792) · [Evaluated commit](https://github.com/sauravsingla/OligoArk/commit/d7b0f89a3e7521f02381a872be92b37342607f71). Compact-v3 verified **21/60** exact recoveries; DNA Fountain produced **no verified outcomes** before its worker budget expired. **Not a complete four-method comparison.**
- **1 GiB streaming-storage experiment (completed):** [Scale acceptance evidence](docs/storage-scale-acceptance-2026-10-07.md) · [Streaming architecture and scope](docs/scalable-storage.md). Exact SHA-256 recovery under clean, 1%, and 5% controlled strand-loss conditions; **separate from the 152-nt matched-codec comparison**.
- **Published physical-read reconstruction (reference level only):** [Microsoft CNR](docs/external-cnr-benchmark.md) · [Grass et al.](docs/external-grass-benchmark.md) · [LCRC HFS-11.7K](docs/external-lcrc-benchmark.md) · [DNAformer Pilot](docs/external-dnaformer-pilot-benchmark.md). Published reads were **not encoded with OligoArk**; the small differences versus BBS are not statistically significant on the selected subsets, and BBS is substantially faster.

**Interpretation:** Codec-recovery, streaming-archive, and published-read reconstruction results use **different protocols** and must not be combined into one success rate. Exact archive recovery is checked against the original payload's SHA-256 digest; reference-strand reconstruction is evaluated separately.

## Limitations and physical validation

- **152-nt indel recovery at 10 MiB remains 0/10** under the specified timeout; substitution and mixed errors also limit the 100 MiB results.
- The **100 MiB four-method comparison is incomplete** because the DNA Fountain worker timed out before producing verified trial results.
- The repository has synthesis-ready FASTA/CSV and a [wet-lab validation protocol](docs/wet-lab-validation.md), but **physical OligoArk archive status remains: prepared, not physically executed**. A valid claim requires OligoArk encoding → DNA synthesis → sequencing → archive reconstruction → exact SHA-256 match.

## Further reading

[Research scope and prior work](docs/research.md) · [Architecture](docs/architecture.md) · [Benchmark methodology](docs/benchmarking.md) · [Configuration](docs/configuration.md) · [Contributing](CONTRIBUTING.md).

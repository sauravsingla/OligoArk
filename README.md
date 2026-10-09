# 🧬 OligoArk

**DNA data-storage research you can run, stress-test, and reproduce.**

An open-source Python framework to encode files into **DNA-like strands**, simulate noisy reads, reconstruct data, and verify recovery with **SHA-256**.

[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](pyproject.toml)
[![MIT License](https://img.shields.io/badge/License-MIT-2ea44f)](LICENSE)
[![Research software](https://img.shields.io/badge/Status-Research%20software-orange)](docs/research.md)

![OligoArk animated software workflow: file, DNA encoding, simulated errors, and verified recovery](assets/oligoark-workflow.svg)

*Encode → Simulate → Reconstruct → Verify · Software simulation, not physical DNA synthesis.*

[**Quick start**](#quick-start) · [**Benchmark**](#featured-benchmark-10-mib) · [**Results**](#benchmark-reports) · [**Contribute**](CONTRIBUTING.md)

## Quick start

Requires Python 3.10+ and Git; no laboratory equipment needed.

```bash
git clone https://github.com/sauravsingla/OligoArk.git
cd OligoArk
python -m pip install -e .
python examples/quick_start.py
```

Expected: `PASS: original data recovered exactly` (deterministic seed 42, simulated channel). [Demo source](examples/quick_start.py) · [End-to-end example](examples/end_to_end.py).

## What you can explore

- **Noisy-read reconstruction:** confidence-fusion consensus and bounded edits. [Code](src/oligoark/reconstruct.py) · [Physical-read evaluation](docs/external-cnr-benchmark.md).
- **Adaptive coding:** compare redundancy and reconstruction strategies on held-out simulation seeds. [Research](docs/research.md).
- **Reproducible benchmarks:** inspect exact recoveries, failures, overhead, runtime and separate streaming-archive results. [Methods](docs/benchmarking.md) · [Streaming storage](docs/scalable-storage.md).

## Featured benchmark: 10 MiB

**Matched simulated codec comparison:** one 10 MiB payload, 152-nt maximum strands, ~25% redundancy, identical channel settings/seeds, **10 trials per condition**. Results are complete, SHA-256-verified archive recoveries.

| Simulated condition | OligoArk compact-v3 | DNA Fountain-style | Goldman + XOR |
| --- | ---: | ---: | ---: |
| Clean | **10/10** | 10/10 | 10/10 |
| 1% / 5% strand loss | **10/10 each** | 10/10 each | 0/10 each |
| Substitutions | **9/10** | 0/10 | 0/10 |
| Insertions/deletions | **0/10 (timeouts)** | 0/10 | 0/10 |
| Mixed errors | **9/10** | 0/10 | 0/10 |

Compact-v3's **10 indel trials exceeded the 75-second deadline**. The fourth tested method (OligoArk compact hybrid) recovered 10/10 clean and 0/10 noisy trials. Baselines are **independent approximations**, not bit-compatible published implementations; this software result does **not** establish general superiority. [Baseline provenance](docs/dna-fountain-baseline.md).

### Benchmark reports

- **10 MiB codec comparison:** [Baseline](docs/dna-fountain-baseline.md) · [Methodology](docs/benchmarking.md).
- **100 MiB codec research (partial):** [Workflow and failure logs](https://github.com/sauravsingla/OligoArk/actions/runs/37766665792) · [Commit](https://github.com/sauravsingla/OligoArk/commit/d7b0f89a3e7521f02381a872be92b37342607f71). Compact-v3: **21/60**; DNA Fountain: **no verified results** before timeout.
- **1 GiB streaming archive (separate experiment):** [Scale evidence](docs/storage-scale-acceptance-2026-10-07.md) · [Design](docs/scalable-storage.md). Exact recovery in clean, 1%, and 5% controlled-loss tests.
- **Published physical reads (reference-level reconstruction only):** [Microsoft CNR](docs/external-cnr-benchmark.md) · [Grass](docs/external-grass-benchmark.md) · [LCRC](docs/external-lcrc-benchmark.md) · [DNAformer](docs/external-dnaformer-pilot-benchmark.md). Small differences from BBS were not statistically significant on selected subsets; BBS was faster.
- **Held-out v0.6 ablations:** [1,680-trial study](docs/research.md) · [Reproduction](docs/benchmarking.md).

**These experiments use different protocols and cannot be combined into one success rate.**

## Research limitations

OligoArk is **simulation-first**. Although published physical reads have been tested for reference reconstruction, **no OligoArk-encoded archive has completed synthesis → sequencing → exact recovery**. [Wet-lab protocol](docs/wet-lab-validation.md) · [Physical pilot plan (not yet executed)](docs/physical-pilot-v1.md).

Indel recovery and 100 MiB noisy-channel results remain limited; the four-method 100 MiB comparison is unfinished. OligoArk combines existing techniques rather than claiming to invent fountain coding, Reed–Solomon protection, graph reconstruction or sequence screening.

## Explore and contribute

[Architecture](docs/architecture.md) · [Configuration](docs/configuration.md) · [Research](docs/research.md) · [Benchmarking](docs/benchmarking.md) · [Contributing](CONTRIBUTING.md) · [Open an issue](https://github.com/sauravsingla/OligoArk/issues)

**License:** [MIT](LICENSE).

# OligoArk 🧬

[![CI](https://github.com/sauravsingla/OligoArk/actions/workflows/ci.yml/badge.svg)](https://github.com/sauravsingla/OligoArk/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/oligoark.svg)](https://pypi.org/project/oligoark/)

**Research software for adaptive DNA-data encoding, error simulation, and reliable reconstruction.**

OligoArk converts digital data into DNA-like sequences, simulates errors such as substitutions, insertions, deletions and missing reads, reconstructs the data with multiple strategies, and verifies exact recovery using SHA-256.

> OligoArk is a **software/simulation research framework**. It does not physically synthesize or store DNA and makes no wet-lab performance claim.

## v0.6 physical benchmark results

| Physical dataset | Reads / strand | **OligoArk confidence fusion** | Pinned BBS |
| --- | ---: | ---: | ---: |
| Microsoft CNR (Nanopore) | 5 | **73/96 (76.0%)** | 72–74/96 |
| Microsoft CNR (Nanopore) | 10 | **93/96 (96.9%)** | **93/96 (96.9%)** |
| Grass et al. (Illumina) | 5 | **94/96 (97.9%)** | 90/96 (93.8%) |
| Grass et al. (Illumina) | 10 | **96/96 (100%)** | 95/96 (99.0%) |
| LCRC HFS-11.7K (Illumina PE150) | 5 | **96/96 (100%)** | **96/96 (100%)** |
| LCRC HFS-11.7K (Illumina PE150) | 10 | **96/96 (100%)** | **96/96 (100%)** |

The **same frozen confidence-fusion settings** were used across all three physical benchmarks; only the known strand length changed mechanically (110 nt CNR, 117 nt Grass, 200 nt LCRC). On the third independent LCRC experiment, both OligoArk and BBS reconstructed all 96 held-out strands exactly at 5 and 10 reads. BBS remains substantially faster. These are reference-level physical-read reconstruction results, not end-to-end external archive decoding or a general state-of-the-art claim.

See [CNR benchmark](docs/external-cnr-benchmark.md), [Grass benchmark](docs/external-grass-benchmark.md), and [LCRC benchmark](docs/external-lcrc-benchmark.md).

## Install

```bash
pip install oligoark
```

## Quick start

```bash
printf 'OligoArk demo data\n' > demo.txt
oligoark archive demo.txt --output demo.oligoark.json
oligoark recover demo.oligoark.json --output recovered.txt
cmp demo.txt recovered.txt
```

## Core capabilities

- Adaptive encoding and redundancy selection
- Substitution, insertion/deletion, dropout and duplication simulation
- Direct, medoid, graph/alignment, iterative-trace, multi-start, targeted and confidence-fusion reconstruction
- SHA-256-verified recovery
- Reproducible held-out benchmarking and optimizer evaluation

## More details

See [Research](docs/research.md), [Benchmarking](docs/benchmarking.md), [Architecture](docs/architecture.md), and [Configuration](docs/configuration.md).

## License

MIT — see [LICENSE](LICENSE).

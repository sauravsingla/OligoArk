# OligoArk 🧬

[![CI](https://github.com/sauravsingla/OligoArk/actions/workflows/ci.yml/badge.svg)](https://github.com/sauravsingla/OligoArk/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/oligoark.svg)](https://pypi.org/project/oligoark/)

**Research software for adaptive DNA-data encoding, error simulation, and reliable reconstruction.**

OligoArk converts digital data into DNA-like sequences, simulates errors such as substitutions, insertions, deletions and missing reads, reconstructs the data with multiple strategies, and verifies exact recovery using SHA-256.

> OligoArk is a **software/simulation research framework**. It does not physically synthesize or store DNA and makes no wet-lab performance claim.

## v0.6 key result

- **1,680** untouched held-out simulation trials.
- Graph/alignment and iterative-trace reconstruction: **223/240 = 92.9%** overall recovery.
- Moderate-indel recovery: **30/30** at both tested read-coverage levels.
- Graph/trace reconstruction rescued **38/38** paired direct-method failures in the two moderate-indel regimes, with **0 regressions**.
- Combined robust optimizer: **221/240 = 92.1%** overall recovery.

The iterative-trace method is substantially slower, learned policies did not outperform the simpler heuristic, and a physical-read CNR benchmark found the external BBS method stronger than OligoArk reconstruction at 5- and 10-read coverage. End-to-end decoding of an external DNA archive has not been validated.

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
- Direct, medoid, graph/alignment and iterative-trace reconstruction
- SHA-256-verified recovery
- Reproducible held-out benchmarking and optimizer evaluation

## More details

See [Research](docs/research.md), [Benchmarking](docs/benchmarking.md), [Architecture](docs/architecture.md), and [Configuration](docs/configuration.md).

## License

MIT — see [LICENSE](LICENSE).

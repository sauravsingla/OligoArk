# OligoArk 🧬

[![CI](https://github.com/sauravsingla/OligoArk/actions/workflows/ci.yml/badge.svg)](https://github.com/sauravsingla/OligoArk/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/oligoark.svg)](https://pypi.org/project/oligoark/)

### Encode files into DNA-like strands, simulate storage errors, reconstruct them, and verify exact recovery.

**OligoArk** is an open-source research framework for DNA archival storage. It combines encoding, redundancy, noisy-channel simulation, multi-read reconstruction and SHA-256 verification in one reproducible pipeline.

## Try it in 30 seconds

```bash
pip install oligoark
printf 'OligoArk demo data\n' > demo.txt
oligoark archive demo.txt --output demo.oligoark.json
oligoark recover demo.oligoark.json --output recovered.txt
cmp demo.txt recovered.txt
```

A successful `cmp` confirms byte-for-byte recovery.

## Large files and scale validation

For large files, use the compact bounded-memory streaming container:

```bash
oligoark archive-stream large.bin --output large.oligoark.bin --profile scale-1024
oligoark recover-stream large.oligoark.bin --output large.recovered.bin
cmp large.bin large.recovered.bin
```

The 100 MiB software acceptance milestone has been achieved for clean, 1% controlled-dropout,
and 5% controlled-dropout recovery with exact SHA-256 verification. The validated clean path
uses bounded-memory streaming; detailed density, throughput, peak-RSS, redundancy, negative
noise results, and claim boundaries are recorded in
[the acceptance evidence](docs/storage-scale-acceptance-2026-10-07.md). This is software
archive evidence, not a wet-lab end-to-end storage claim.

## Why it exists

DNA storage research often evaluates encoding, channel errors, reconstruction and integrity separately. OligoArk makes those stages runnable together so researchers can ask practical questions such as:

- How much redundancy is needed under a given dropout/error profile?
- Which reconstruction strategy works best with only a few noisy reads?
- Does a recovered archive exactly match the original bytes?
- Do reconstruction settings generalize across independent physical-read datasets?

Its current research focus is **low-compute confidence-fusion reconstruction**: deterministic reconstruction signals are combined with bounded confidence-guided repair, without requiring neural models, GPUs or unrestricted search.

> **Scope:** OligoArk is software research, not a wet-lab DNA-storage platform. Physical-read benchmarks use published sequencing data and do not by themselves establish end-to-end wet-lab archival performance or a general state-of-the-art claim.

**Best way to help:** try a new public DNA-storage dataset, reproduce a benchmark, contribute a channel/reconstruction method, or report a failure case.

## What is distinctive

OligoArk is not presented as a new biochemical storage medium or a replacement for established DNA-storage codecs. Its main research contribution is the **integration and validation** of several ideas in one reproducible system:

- **Confidence-fusion reconstruction** combines multi-start, alignment, trace, observed-read distance, q-gram consistency, and local ambiguity evidence in a tightly bounded deterministic candidate search.
- **Leakage-controlled evaluation** separates calibration, development, and untouched held-out subsets so reconstruction parameters are frozen before final evaluation.
- **Cross-dataset generalization** tests the same frozen confidence-fusion settings on multiple independent physical DNA-storage datasets, changing only mechanically required strand length.
- **Measured adaptive codec selection** evaluates candidate redundancy/reconstruction configurations through the real encode → corrupt → recover → SHA-256 verification path instead of relying only on proxy scores.
- **Explicit negative results and claim boundaries** are retained so improvements, regressions, runtime trade-offs, and failed learning approaches remain visible.

> OligoArk is **software research**, not a wet-lab DNA-storage platform. Physical benchmark results below measure reference-level reconstruction from published sequencing reads; they are not end-to-end external archive decodes and do not establish a general state-of-the-art claim.

## Physical benchmark results

| Physical dataset | Reads / strand | **OligoArk confidence fusion** | Pinned BBS |
| --- | ---: | ---: | ---: |
| Microsoft CNR (Nanopore) | 5 | **73/96 (76.0%)** | 72–74/96 |
| Microsoft CNR (Nanopore) | 10 | **93/96 (96.9%)** | **93/96 (96.9%)** |
| Grass et al. (Illumina) | 5 | **94/96 (97.9%)** | 90/96 (93.8%) |
| Grass et al. (Illumina) | 10 | **96/96 (100%)** | 95/96 (99.0%) |
| LCRC HFS-11.7K (Illumina PE150) | 5 | **96/96 (100%)** | **96/96 (100%)** |
| LCRC HFS-11.7K (Illumina PE150) | 10 | **96/96 (100%)** | **96/96 (100%)** |
| DNAformer Pilot (Illumina MiSeq) | 5 | **96/96 (100%)** | **96/96 (100%)** |
| DNAformer Pilot (Illumina MiSeq) | 10 | **96/96 (100%)** | **96/96 (100%)** |

The **same frozen confidence-fusion settings** were used across all four physical benchmarks; only known strand length changed mechanically: 110 nt for CNR, 117 nt for Grass, 200 nt for LCRC, and 140 nt for DNAformer Pilot. On the DNAformer Pilot held-out split, both OligoArk and BBS reconstructed all 96 strands exactly at 5 and 10 reads; at 1 read both reached 83/96. BBS remains substantially faster.

See the detailed benchmark reports for [CNR](docs/external-cnr-benchmark.md), [Grass](docs/external-grass-benchmark.md), [LCRC](docs/external-lcrc-benchmark.md), and [DNAformer Pilot](docs/external-dnaformer-pilot-benchmark.md).

## Research capabilities

- Adaptive encoding and redundancy selection
- XOR, LT-style fountain, and hybrid redundancy strategies
- Substitution, insertion/deletion, dropout, and duplication simulation
- Direct, medoid, graph/alignment, iterative-trace, multi-start, targeted, and confidence-fusion reconstruction
- Constraint-aware DNA-like sequence generation
- Deterministic experiment selection and held-out evaluation
- SHA-256-verified recovery and reproducible benchmark artifacts
- Policy-learning experiments with retained negative results
- Physical-read adapters and external baseline comparison
- Bounded-memory compact binary archives for large-file scale validation
- 152/200/248-nt physical strand profiles and clean-room DNA Fountain comparison

## Research philosophy

OligoArk separates **simulation evidence**, **physical-read reconstruction evidence**, and **external published evidence**. It avoids treating software channel simulations as wet-lab validation, does not use hidden references during reconstruction candidate selection, and preserves negative experiments when a proposed method fails to improve untouched held-out results.

For the detailed novelty, validation design, ablations, prior-work positioning, and limitations, see [Research](docs/research.md).

## Documentation

- [Research scope and novelty](docs/research.md)
- [Benchmarking](docs/benchmarking.md)
- [Architecture](docs/architecture.md)
- [Configuration](docs/configuration.md)
- [Scalable archival storage](docs/scalable-storage.md)
- [DNA Fountain comparison](docs/dna-fountain-baseline.md)
- [100 MiB acceptance evidence](docs/storage-scale-acceptance-2026-10-07.md)

## License

MIT — see [LICENSE](LICENSE).

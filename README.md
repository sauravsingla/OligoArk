# OligoArk 🧬

[![CI](https://github.com/sauravsingla/OligoArk/actions/workflows/ci.yml/badge.svg)](https://github.com/sauravsingla/OligoArk/actions/workflows/ci.yml) [![PyPI](https://img.shields.io/pypi/v/oligoark.svg)](https://pypi.org/project/oligoark/)

**AI-Native DNA Archival Storage** — an open-source research framework for adaptive DNA encoding, software channel simulation, graph-assisted reconstruction, empirical policy learning, and explainable heterogeneous storage tiering.

**What the name means:** “Oligo” refers to oligonucleotides—short DNA strands—and “Ark” represents safeguarding information for the future. Together, **OligoArk** expresses the idea of preserving digital information in DNA for long-term archival storage.

> **Research status:** alpha. OligoArk produces DNA-like sequence encodings and software simulations. It does **not** claim wet-lab validation, present-day commercial DNA-storage economics, or physical-media performance.

## Why OligoArk

DNA data storage research has demonstrated high-density archival concepts, random access, error correction, and long-horizon preservation. A practical software research stack also needs to decide **when** a future DNA tier might be appropriate, **how** encoding/redundancy should adapt to channel conditions, and **how** noisy/duplicated reads should be reconstructed. OligoArk treats these questions as one reproducible system rather than only mapping bits to `A/C/G/T`.

## Research contributions

OligoArk v0.6 exposes seven testable research layers through one **AI-Native Archival Intelligence Layer**:

1. **Explainable heterogeneous storage tiering** — evaluates SSD, object archive, tape, and an explicitly experimental `dna_future` tier from retention, access, mutability, durability, retrieval urgency, redundancy, energy, and user-supplied normalized economic assumptions.
2. **Search-based adaptive codec optimisation** — actually encodes, simulates, reconstructs, SHA-256-verifies, measures and ranks candidate configurations across chunk size, Reed–Solomon strength, XOR/fountain/hybrid redundancy, sequence constraints and reconstruction mode.
3. **Explicit graph + multi-trace reconstruction** — sequencing reads are graph nodes, qualifying similarities are weighted edges, connected components form clusters, and medoid, single-pass alignment, plus iterative trace consensus provide deterministic indel-aware ablations.
4. **Policy learning** — instance-based, deterministic ridge-regression, and deterministic RBF-kernel utility baselines learn only from reproducible observations; no proprietary model is required.
5. **Lifecycle-aware archival intelligence** — every tier returns a decomposable score, and caller-supplied lifecycle cost/energy/latency can influence both tier selection and codec optimisation.
6. **Robust held-out validation** — calibration seeds are disjoint from evaluation seeds, robust search penalizes cross-seed instability, calibration rotates across several payload contents, and learning is tested on unseen seeds and channel regimes.
7. **Physical-read evaluation pathway** — FASTA/FASTQ reads can be compared against an explicit reference-oligo FASTA with provenance metadata; OligoArk never infers an external archive format or labels such reconstruction as OligoArk end-to-end decoding.

`plan_archive()` remains a lightweight explainable heuristic. `optimize_archive_plan()` performs a real measured search, while `evaluate_optimized_archive_plan()` freezes the calibrated winner and measures it on disjoint held-out seeds. The CLI and API return the selected tier, redundancy scheme, sequence constraints, reconstruction strategy, search metadata, and per-candidate objective breakdowns. Reconstruction remains extensible through typed `ReadReconstructor` and `EdgeScorer` interfaces so future PyTorch/PyTorch-Geometric models can plug in without changing the archive format.

The independently implemented **fountain-style seeded XOR/peeling baseline is now integrated into the real archive/recovery pipeline** and can be selected as `fountain` or combined with XOR as `hybrid`. It is not represented as the published DNA Fountain implementation.

## Architecture

```mermaid
flowchart LR
  A[File / bytes] --> B[Chunker]
  B --> C[Reed-Solomon ECC]
  C --> D[Hard GC/homopolymer constraint search]
  D --> E[DNA strand framing + XOR/Fountain/Hybrid]
  E --> F[Software archive]
  F --> G[Channel simulator]
  G --> H[Direct verified decode]
  H -->|insufficient| I[Multi-threshold graph + alignment/trace consensus]
  I --> J[Decode + XOR erasure recovery]
  H -->|sufficient| J
  J --> K[SHA-256 verification]

  L[Channel profile] --> M[Heuristic codec policy]
  N[Candidate search space] --> O[Balanced robust/full-grid measured optimizer]
  L --> O
  T[Calibration seeds] --> O
  O --> U[Frozen codec plan]
  V[Disjoint evaluation seeds] --> W[Held-out validation]
  U --> W
  U --> B

  P[Workload profile] --> Q[AI-Native Archival Intelligence]
  R[User economic/lifecycle assumptions] --> Q
  L --> Q
  Q --> S[Decomposed tier + codec plan]
```

See [`docs/architecture.md`](docs/architecture.md) for module-level details.

## Install

For normal use, install the published package from PyPI:

```bash
pip install oligoark
```

Optional API support:

```bash
pip install "oligoark[api]"
```

For contributors and research development:

```bash
git clone https://github.com/sauravsingla/OligoArk.git
cd OligoArk
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -U pip
pip install -e ".[dev,api,bench]"
```

The core package has **no mandatory third-party runtime dependency**. FastAPI, plotting, and development tools are optional extras.

## Quick start

```bash
printf 'OligoArk demo data\n' > demo.txt
oligoark archive demo.txt --output demo.oligoark.json
oligoark inspect demo.oligoark.json
oligoark recover demo.oligoark.json --output recovered.txt
cmp demo.txt recovered.txt
```

`archive` reports measured software facts such as strand count, encoded nucleotide count, GC statistics, homopolymer maximum, logical bits/nucleotide, and SHA-256. These are **software encoding measurements**, not physical synthesis claims.

### Simulate errors and reconstruct reads

```bash
oligoark simulate demo.oligoark.json \
  --substitution 0.001 \
  --dropout 0.01 \
  --duplicate 0.20 \
  --seed 7 \
  --output reads.txt

oligoark recover-reads \
  demo.oligoark.json reads.txt \
  --similarity-threshold 0.90 \
  --output recovered-from-reads.txt
```

The recovery report states whether direct decoding succeeded or graph reconstruction was needed. Recovery is accepted only after the original SHA-256 is reproduced.

### Adaptive codec recommendation

```bash
oligoark policy \
  --substitution 0.01 \
  --deletion 0.002 \
  --dropout 0.08 \
  --durability-priority 0.9 \
  --storage-overhead-priority 0.1
```

### Heuristic archival plan

```bash
oligoark plan \
  --retention-years 100 \
  --accesses-per-year 0.1 \
  --durability-priority 1.0 \
  --energy-priority 0.8 \
  --substitution 0.01 \
  --dropout 0.05
```

This command is intentionally labeled heuristic. For a real measured search, use:

```bash
oligoark optimize-plan demo.txt \
  --retention-years 100 \
  --substitution 0.01 \
  --dropout 0.05 \
  --calibration-seeds 9401,9402,9403,9404,9405,9406 \
  --evaluation-seeds 31001,31002,31003,31004 \
  --search-method balanced_robust \
  --search-seed 6060 \
  --reconstruction-modes direct,graph,trace \
  --max-candidates 36
```

`optimize-plan` evaluates actual archive/simulation/recovery candidates and accepts success only after SHA-256 verification. `--seeds` is retained as a backward-compatible alias for calibration seeds; when `--evaluation-seeds` is supplied, the calibrated winner is frozen before unseen-seed evaluation.

### Storage-tier recommendation

```bash
oligoark recommend \
  --retention-years 100 \
  --accesses-per-year 0.1 \
  --mutability 0.0 \
  --retrieval-urgency 0.1 \
  --durability-priority 1.0 \
  --energy-priority 0.8 \
  --redundancy-priority 0.8 \
  --cost-priority 0.5
```

The default economics are deliberately **neutral across tiers**; OligoArk does not invent vendor pricing or future DNA costs. Real normalized inputs can be supplied through the Python API, REST API, or `--economics-json`. For example:

```json
{
  "storage_cost_index": {"ssd": 0.8, "object_archive": 0.3, "tape": 0.2, "dna_future": 0.5},
  "retrieval_cost_index": {"ssd": 0.1, "object_archive": 0.4, "tape": 0.7, "dna_future": 0.9}
}
```

Lower normalized values mean lower assumed cost. These values are user inputs, not prices asserted by OligoArk. For explicit lifecycle analysis, `--lifecycle-json` accepts caller-supplied per-tier storage cost, retrieval cost, idle/retrieval energy and retrieval latency values; OligoArk provides no fabricated lifecycle defaults.

## Python SDK

```python
from oligoark import ArchiveConfig, archive_bytes, archive_statistics, recover_bytes

payload = b"long-lived research artifact"
archive = archive_bytes(payload, ArchiveConfig(rs_nsym=12))
print(archive_statistics(archive).to_dict())
recovered = recover_bytes(archive)
assert recovered == payload
```

### Empirical policy-learning baseline

```python
from oligoark import EmpiricalPolicyModel, PolicyObservation

model = EmpiricalPolicyModel().fit(observations)
recommendation = model.recommend(target_channel)
print(recommendation.policy)
```

`observations` must come from explicitly labeled simulated or measured experiments. The model is transparent instance-based learning, not a pretrained black box.

### AI-Native Archival Intelligence Layer

```python
from oligoark import ChannelProfile, WorkloadProfile, plan_archive

workload = WorkloadProfile(
    retention_years=100,
    accesses_per_year=0.1,
    mutability=0.0,
    retrieval_urgency=0.1,
    durability_priority=1.0,
    energy_priority=0.8,
)
channel = ChannelProfile(substitution_rate=0.01, dropout_rate=0.05)

plan = plan_archive(workload, channel)
print(plan.to_dict())
```

### Reconstruction extension interface

Custom reconstruction methods implement `ReadReconstructor`; custom graph edge models implement `EdgeScorer`. The default `GraphConsensusReconstructor` uses normalized Levenshtein similarity. A future PyTorch/PyTorch-Geometric scorer can be injected with `use_qgram_prefilter=False` without changing encoding, archive metadata, or recovery verification.

## REST API

```bash
uvicorn oligoark.api:app --host 127.0.0.1 --port 8000
```

Open `/docs` for the generated OpenAPI UI. Endpoints include:

- `GET /health`
- `POST /encode`
- `POST /recover`
- `POST /recover-reads`
- `POST /reconstruction-diagnostics`
- `POST /simulate`
- `POST /recommend`
- `POST /policy`
- `POST /plan`
- `POST /optimize-plan`

Runtime settings can be supplied via `OLIGOARK_LOG_LEVEL`, `OLIGOARK_RECONSTRUCTION_THRESHOLD`, and `OLIGOARK_MAX_API_PAYLOAD_BYTES`. See [`docs/configuration.md`](docs/configuration.md).

## Reproducible benchmark

```bash
python benchmarks/run_benchmark.py
python benchmarks/run_experiments.py --profile smoke
python benchmarks/run_verified_demo.py
bandit -r src/oligoark -q
pip-audit --local --skip-editable
```

The benchmark writes:

- `benchmark-results/results.json`
- `benchmark-results/results.csv`
- `benchmark-results/metadata.json`
- `benchmark-results/recovery_by_regime.png` when Matplotlib is installed
- `benchmark-results/overhead_by_regime.png` when Matplotlib is installed

Metadata records the OligoArk version, Python version, platform, deterministic seed, payload size, payload SHA-256, and claim scope.

### Experiment and ablation framework

The continuity benchmark is retained, while v0.6 adds untouched leakage-controlled multi-trace ablations:

```bash
python benchmarks/run_experiments.py --profile smoke
python benchmarks/run_experiments.py --profile publication
```

The experiment framework compares **fixed**, **adaptive**, **adaptive + fountain/hybrid redundancy**, **adaptive + medoid graph consensus**, **adaptive + graph/alignment reconstruction**, **adaptive + iterative trace reconstruction**, and the **combined robust measured optimizer**. Calibration seeds are disjoint from evaluation seeds. It records raw CSV/JSON, aggregated CSV/JSON, Wilson 95% recovery intervals, paired differences versus fixed on identical realizations, runtime, encoded overhead, strand/read counts, graph rescue, candidate-selection metadata, and reproducible plots. The publication workflow deterministically separates combined-optimizer calibration from the non-optimizer baselines, shards the baselines by payload size, channel regime, and two disjoint evaluation-seed partitions, then aggregates every failure and success without reducing the scientific design.

The same artifacts feed a held-out learning comparison between the deterministic heuristic, empirical instance-based model, ridge-regression model, RBF-kernel model, adaptive+fountain baseline, and measured-search result:

```bash
python benchmarks/run_learning_evaluation.py --experiment-dir publication-results
python benchmarks/run_graph_rescue.py
```

`benchmarks/run_verified_demo.py` separately exercises adaptive optimisation followed by substitution, indel, dropout and mixed corruption with graph/alignment reconstruction and SHA-256 verification.

See [`docs/benchmarking.md`](docs/benchmarking.md) for interpretation rules.

## Research Validation

The v0.5 publication workflow was executed on commit `180618c9f5bdc1260d00c0b60a09dd6442c1a569` with Python 3.13.15 on Linux. It evaluated **960 held-out trials**: 3 payload sizes (512, 2048 and 8192 bytes) × 8 channel regimes × 5 strategies × 8 evaluation seeds. The combined optimizer used calibration seeds `9201..9204`, completely disjoint from evaluation seeds `2026..2033`, and evaluated 24 deterministically balanced candidates from a 360-candidate logical search space per payload/regime. The artifact also preserves 576 calibration-candidate evaluations and their objective breakdowns.

Across all 192 held-out trials per strategy, SHA-256-verified recovery was: **adaptive+fountain/hybrid 66.7%** (95% Wilson CI 59.7–73.0%), **adaptive 59.9%** (52.8–66.6%), **adaptive+graph 59.9%** (52.8–66.6%), **combined measured optimizer 57.3%** (50.2–64.1%), and **fixed 45.8%** (38.9–52.9%). Mean encoded-overhead ratios were 2.114, 1.557, 1.557, 1.895 and 1.492 respectively; mean software runtimes were 0.959 s, 0.728 s, 3.379 s, 1.841 s and 0.036 s in that run environment. These runtimes are not cross-machine performance claims.

The result is deliberately reported even where it is negative. The combined optimizer improved over fixed overall, but **did not beat adaptive+fountain/hybrid** and showed calibration-to-held-out degradation in dropout, indel and mixed regimes. For example, at 10% configured dropout it selected a low-overhead XOR configuration after 100% calibration recovery, but held-out recovery fell with payload size; this is an observed generalization limitation of the current calibration budget rather than a hidden success.

Controlled reconstruction produced the required verified rescue chain in **2/2 constructed cases**. In both cases direct recovery failed, an explicit 7-node/21-edge graph formed one component, consensus was generated, and normal archive recovery passed SHA-256. In the substitution case both medoid and alignment consensus rescued recovery; in the insertion/deletion case **medoid failed while alignment-aware consensus succeeded**. In the broader 960-trial publication sweep, however, graph reconstruction did not produce additional direct-failure rescues, so the controlled rescue result is evidence of capability rather than evidence of broad superiority.

Held-out policy learning trained on seeds 2026–2029 and the first four channel regimes, then tested on seeds 2030–2033 and four unseen channel regimes. On 48 held-out groups, the heuristic and ridge model each recovered **41.7%** with mean regret 0.0520; the empirical model recovered **37.5%** with regret 0.0740; and the measured-search result recovered **39.6%** with regret 0.1630. The learned ridge baseline therefore did **not** outperform the heuristic in this experiment.

All numbers above are software/simulation results from the preserved GitHub Actions publication artifact, not wet-lab DNA-storage measurements or physical-media claims.

### v0.6 validation design

v0.6 deliberately does **not** reuse the v0.5 publication test seeds. The new publication profile uses calibration seeds `9401–9406`, untouched evaluation seeds `31001–31010`, three payload sizes, eight channel/coverage regimes, seven strategies, three independent calibration payload contents, and `balanced_robust` candidate scoring. Moderate-indel channels are evaluated at both five-copy and eight-copy per-strand coverage so reconstruction is tested with genuine multi-trace evidence rather than only one optional duplicate.

The DNA-Aeon provenance manifest records public sequencing accessions, but the repository does not bundle a verified external read-to-reference oligo mapping. `oligoark physical-evaluate` therefore requires an explicit reference FASTA and reports only reference-reconstruction accuracy; it does not claim end-to-end decoding of an external DNA-storage system.

The v0.5 results below remain the historical baseline until the v0.6 untouched publication workflow is executed and its artifact is inspected.

## Tests and quality gates

```bash
pytest --cov=oligoark --cov-report=term-missing
ruff check .
mypy src/oligoark
python -m build
python -m compileall -q src tests examples benchmarks
python examples/end_to_end.py
python examples/noisy_recovery.py
python examples/policy_learning.py
python examples/archival_plan.py
python examples/custom_reconstructor.py
python benchmarks/run_benchmark.py
```

The test suite includes unit, integration, API, CLI, reconstruction plug-in, archival-intelligence, ECC, policy-learning, archive-validation, version-synchronization, and Hypothesis property-based tests. CI runs supported Python versions independently, enforces at least 80% coverage, builds the package and Docker image, runs all examples, executes the deterministic benchmark, verifies CSV/JSON/metadata outputs, generates both benchmark plots, and uploads benchmark artifacts. Tag pushes matching `v*` build installable release artifacts in a separate release workflow.

## Scientific assumptions and limitations

- The base codec is a reversible 2-bit mapping wrapped in protected frames; it is not claimed to be capacity-optimal.
- GC range and maximum homopolymer length are configurable **hard software constraints** enforced by deterministic mask search. They are research constraints, not claims about a particular synthesis platform.
- Reed–Solomon protects frame bytes; archive-level redundancy can be `none`, `xor`, `fountain`, or `hybrid`. The fountain implementation is an LT-style research baseline, not DNA Fountain.
- The simulator uses independent substitution/insertion/deletion/dropout probabilities plus explicit read coverage and optional extra duplication. Real channels can be correlated and platform-specific.
- Reconstruction baselines use explicit weighted graphs, connected components, medoid/alignment consensus and a deterministic iterative trace-consensus extension. These remain research baselines and are not claimed equivalent to HEDGES or state-of-the-art trace reconstruction.
- The `dna_future` storage tier is scenario analysis only. It is not a statement that DNA is presently cheaper, faster, or operationally superior to SSD/object/tape.
- Economic defaults are neutral. Production decisions require externally sourced and time-appropriate cost, energy, durability, and retrieval assumptions.
- The reference API is a research service, not a hardened multi-tenant production system.

## Research roadmap

1. Evaluate learned `EdgeScorer`/`ReadReconstructor` plug-ins against the new explicit graph/alignment baseline.
2. Validate the search space and learned policy models on suitably licensed published physical DNA-storage read/error datasets rather than simulation alone.
3. Add platform-specific synthesis/sequencing channel models and correlation structure.
4. Compare the LT-style fountain baseline with published rateless/fountain implementations under common constraints.
5. Run lifecycle sensitivity analyses from dated externally sourced cost/energy assumptions across 10/50/100/500-year horizons.
6. Rust acceleration for codec, edit distance, q-gram indexing, and large-read clustering hot paths.
7. Optional real synthesis/sequencing adapters isolated behind interfaces so simulation results cannot be confused with wet-lab evidence.

## Prior work and attribution

OligoArk is independently implemented and does not vendor or copy another DNA-storage repository. Relevant foundations include:

- Church, Gao & Kosuri (2012), *Next-generation digital information storage in DNA*, **Science**. DOI: `10.1126/science.1226355`.
- Goldman et al. (2013), *Towards practical, high-capacity, low-maintenance information storage in synthesized DNA*, **Nature**. DOI: `10.1038/nature11875`.
- Grass et al. (2015), *Robust Chemical Preservation of Digital Information on DNA in Silica with Error-Correcting Codes*, **Angewandte Chemie International Edition**. DOI: `10.1002/anie.201411378`.
- Erlich & Zielinski (2017), *DNA Fountain enables a robust and efficient storage architecture*, **Science**. DOI: `10.1126/science.aaj2038`.
- Yazdi, Gabrys & Milenkovic (2017), *Portable and Error-Free DNA-Based Data Storage*, **Scientific Reports**. DOI: `10.1038/s41598-017-05188-1`.
- Organick et al. (2018), *Random access in large-scale DNA data storage*, **Nature Biotechnology**. DOI: `10.1038/nbt.4079`.

See [`docs/research.md`](docs/research.md) for research framing and claim boundaries, and [`docs/benchmarking.md`](docs/benchmarking.md) for the exact v0.5 held-out design.

## Security

OligoArk parses untrusted archive-like input and performs potentially expensive reconstruction. Do not expose the reference API directly to untrusted networks without authentication, rate limiting, request limits, and deployment hardening. See [`SECURITY.md`](SECURITY.md).

## Contributing

Contributions are welcome. Read [`CONTRIBUTING.md`](CONTRIBUTING.md), especially the requirement to label results as **measured software result**, **simulation result**, **external published result**, or **hypothesis**.

## License

MIT — see [`LICENSE`](LICENSE).

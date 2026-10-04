# OligoArk 🧬

[![CI](https://github.com/sauravsingla/OligoArk/actions/workflows/ci.yml/badge.svg)](https://github.com/sauravsingla/OligoArk/actions/workflows/ci.yml) [![PyPI](https://img.shields.io/pypi/v/oligoark.svg)](https://pypi.org/project/oligoark/)

**AI-Native DNA Archival Storage** — an open-source research framework for adaptive DNA encoding, software channel simulation, graph-assisted reconstruction, empirical policy learning, and explainable heterogeneous storage tiering.

**What the name means:** “Oligo” refers to oligonucleotides—short DNA strands—and “Ark” represents safeguarding information for the future. Together, **OligoArk** expresses the idea of preserving digital information in DNA for long-term archival storage.

> **Research status:** alpha. OligoArk produces DNA-like sequence encodings and software simulations. It does **not** claim wet-lab validation, present-day commercial DNA-storage economics, or physical-media performance.

## Why OligoArk

DNA data storage research has demonstrated high-density archival concepts, random access, error correction, and long-horizon preservation. A practical software research stack also needs to decide **when** a future DNA tier might be appropriate, **how** encoding/redundancy should adapt to channel conditions, and **how** noisy/duplicated reads should be reconstructed. OligoArk treats these questions as one reproducible system rather than only mapping bits to `A/C/G/T`.

## Research contributions

OligoArk v0.4 exposes five testable research layers through one **AI-Native Archival Intelligence Layer**:

1. **Explainable heterogeneous storage tiering** — evaluates SSD, object archive, tape, and an explicitly experimental `dna_future` tier from retention, access, mutability, durability, retrieval urgency, redundancy, energy, and user-supplied normalized economic assumptions.
2. **Search-based adaptive codec optimisation** — actually encodes, simulates, reconstructs, SHA-256-verifies, measures and ranks candidate configurations across chunk size, Reed–Solomon strength, XOR/fountain/hybrid redundancy, sequence constraints and reconstruction mode.
3. **Explicit graph + alignment reconstruction** — sequencing reads are graph nodes, qualifying similarities are weighted edges, connected components form clusters, and a medoid-anchored global-alignment consensus provides a deterministic indel-aware baseline.
4. **Policy learning** — both an instance-based baseline and a deterministic ridge-regression utility model learn only from caller-supplied reproducible observations; no proprietary model is required.
5. **Lifecycle-aware archival intelligence** — storage tiering can include caller-supplied cost, energy and retrieval-latency assumptions in addition to retention, access, mutability, durability and redundancy priorities.

`plan_archive()` remains a lightweight explainable heuristic. `optimize_archive_plan()` and the `oligoark optimize-plan` command perform a real measured search over candidate configurations and return the selected tier, redundancy scheme, sequence constraints and reconstruction strategy. Reconstruction remains extensible through typed `ReadReconstructor` and `EdgeScorer` interfaces so future PyTorch/PyTorch-Geometric models can plug in without changing the archive format.

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
  H -->|insufficient| I[Explicit weighted graph + alignment consensus]
  I --> J[Decode + XOR erasure recovery]
  H -->|sufficient| J
  J --> K[SHA-256 verification]

  L[Channel profile] --> M[Heuristic codec policy]
  N[Candidate search space] --> O[Measured codec optimizer]
  L --> O
  O --> B

  P[Workload profile] --> Q[AI-Native Archival Intelligence]
  R[User economic assumptions] --> Q
  L --> Q
  Q --> S[Explainable tier + codec plan]
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
  --seeds 2026,2027 \
  --max-candidates 24
```

`optimize-plan` evaluates actual archive/simulation/recovery candidates and accepts success only after SHA-256 verification.

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

The single-seed benchmark is retained for continuity, while v0.4 adds multi-seed ablations:

```bash
python benchmarks/run_experiments.py --profile smoke
python benchmarks/run_experiments.py --profile publication
```

The experiment framework compares **fixed**, **adaptive**, **adaptive + fountain/hybrid redundancy**, **adaptive + graph/alignment reconstruction**, and the **combined measured optimizer**. It records raw CSV/JSON, aggregated CSV/JSON, Wilson 95% recovery intervals, runtime, encoded overhead, strand/read counts and graph usage, plus reproducible plots. CI runs the smoke profile; the heavier publication profile is a manual GitHub Actions workflow.

`benchmarks/run_verified_demo.py` separately exercises adaptive optimisation followed by substitution, indel, dropout and mixed corruption with graph/alignment reconstruction and SHA-256 verification.

See [`docs/benchmarking.md`](docs/benchmarking.md) for interpretation rules.

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
- The simulator uses independent substitution/insertion/deletion/dropout/duplication probabilities. Real channels can be correlated and platform-specific.
- The default reconstruction baseline builds an explicit weighted similarity graph, uses connected components, and applies medoid-anchored global-alignment consensus. It improves the algorithmic treatment of indels, but it is not claimed to match state-of-the-art reconstruction or HEDGES-style indel coding.
- The `dna_future` storage tier is scenario analysis only. It is not a statement that DNA is presently cheaper, faster, or operationally superior to SSD/object/tape.
- Economic defaults are neutral. Production decisions require externally sourced and time-appropriate cost, energy, durability, and retrieval assumptions.
- The reference API is a research service, not a hardened multi-tenant production system.

## Research roadmap

1. Evaluate learned `EdgeScorer`/`ReadReconstructor` plug-ins against the new explicit graph/alignment baseline.
2. Calibrate the search space and learned policy models on published physical DNA-storage datasets rather than simulation alone.
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

See [`docs/research.md`](docs/research.md) for research framing and claim boundaries.

## Security

OligoArk parses untrusted archive-like input and performs potentially expensive reconstruction. Do not expose the reference API directly to untrusted networks without authentication, rate limiting, request limits, and deployment hardening. See [`SECURITY.md`](SECURITY.md).

## Contributing

Contributions are welcome. Read [`CONTRIBUTING.md`](CONTRIBUTING.md), especially the requirement to label results as **measured software result**, **simulation result**, **external published result**, or **hypothesis**.

## License

MIT — see [`LICENSE`](LICENSE).

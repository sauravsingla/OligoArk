# Changelog

## 0.5.0 - 2026-10-04

### Added
- Disjoint calibration and held-out evaluation for the measured optimizer.
- Canonical order-independent balanced candidate sampling plus full-grid search.
- Explicit optimizer terms for redundancy and optional caller-supplied lifecycle cost, energy, retrieval cost, and latency.
- Per-candidate objective contribution breakdowns and search metadata.
- Fully decomposable per-tier storage score contributions.
- Controlled graph-rescue validation comparing direct, medoid-graph, and alignment-graph recovery.
- Graph diagnostics for nodes, candidate pairs, edges, connected components, cluster sizes, consensus lengths, and runtime.
- Experiment-record to policy-learning dataset conversion and held-out learning evaluation.
- Deterministic ridge-model serialization and restoration.
- Paired strategy effects on identical simulated channel realizations.
- Publication plots for error/recovery, overhead/recovery, runtime/recovery, graph rescue, ablation, and optimizer generalization.
- Payload-sharded publication workflow with merged 90-day artifacts.
- CLI/API controls for calibration/evaluation seeds, search method/seed, objective weights, and reconstruction diagnostics.

### Validation
- Executed the publication workflow on commit `180618c9f5bdc1260d00c0b60a09dd6442c1a569`: 960 held-out trials, 576 retained calibration-candidate evaluations, 3 payload sizes, 8 channel regimes, 5 strategies, disjoint calibration/evaluation seeds, paired effects, Wilson intervals, graph-rescue diagnostics, and held-out policy-learning outputs.
- Best overall held-out recovery in this software experiment was adaptive+fountain/hybrid at 66.7% (128/192; Wilson 95% CI 59.7–73.0%). Combined measured search recovered 57.3% (110/192) and did not outperform the simpler redundancy baseline.
- Controlled alignment-graph reconstruction rescued 2/2 direct-failure cases; medoid consensus also rescued the substitution case but failed the insertion/deletion case that alignment consensus recovered.
- Ridge policy learning matched the heuristic on the unseen seed/channel split (41.7% recovery each) rather than outperforming it; negative results are retained.

### Changed
- Version advanced to 0.5.0.
- Publication evaluation seeds are explicitly disjoint from optimizer calibration seeds.
- Lifecycle terms can influence the codec objective when caller-supplied physical assumptions are provided.
- Publication validation distinguishes calibration score from unseen-seed recovery.

## 0.4.0 - 2026-10-04

### Added
- Search-based codec optimisation that actually encodes, simulates, recovers, SHA-256 verifies, measures, and ranks candidate configurations.
- Adaptive redundancy selection across XOR, LT-style fountain, and hybrid archive strategies.
- Fountain symbols integrated into the main archive/recovery path.
- Hard GC-content and homopolymer constraints with deterministic mask search and explicit failure when constraints cannot be satisfied.
- Explicit weighted similarity graphs, connected-component clustering, and medoid-anchored global-alignment consensus.
- Caller-supplied lifecycle storage/retrieval cost, energy, and latency estimates.
- Dependency-free ridge-regression utility policy model alongside the existing empirical baseline.
- Multi-seed/multi-payload ablation experiments with Wilson 95% recovery intervals.
- CI smoke ablations, a verified adaptive corruption/recovery demo, and a manual publication experiment workflow.
- Bandit, pip-audit, and CodeQL security gates.
- CLI/API support for redundancy, hard constraints, lifecycle inputs, and measured optimisation.

### Changed
- Version advanced to 0.4.0.
- `plan_archive()` is explicitly documented as heuristic; `optimize_archive_plan()` is the measured optimizer.
- Graph reconstruction now builds real node/edge graphs instead of representative-only clusters.
- The continuity benchmark's fixed baseline explicitly disables hard sequence constraints for a fair historical comparison.

## 0.3.0 - 2026-10-04

### Added
- Formal `EdgeScorer` and `ReadReconstructor` protocols for future learned/GNN reconstruction.
- `GraphConsensusReconstructor` with pluggable edge scoring and optional q-gram prefilter bypass.
- Pluggable `recover_from_reads()` reconstruction fallback with strategy diagnostics.
- `ArchivalIntelligencePlan`, `objective_from_workload()`, and `plan_archive()` as a unified AI-native planning layer.
- Integrated `oligoark plan` CLI command and `POST /plan` FastAPI endpoint.
- Public examples for archival planning and custom reconstruction scoring.
- Tests for plug-in reconstruction, intelligence planning, API/CLI integration, and package-version synchronization.
- CI execution of the reproducible benchmark with mandatory CSV/JSON/metadata outputs and plot generation/upload.

### Changed
- Package version advanced to 0.3.0.
- API health response now includes the package version.
- Documentation now describes typed future PyTorch/PyTorch-Geometric integration boundaries explicitly.

## 0.2.0 - 2026-10-04

### Added
- Graph-consensus reconstruction is now wired into the noisy-read recovery path.
- Q-gram prefiltering before Levenshtein graph edges to reduce unnecessary edit-distance work.
- Archive inspection statistics including GC range, homopolymer maximum, logical density, and strand counts.
- Neutral/user-configurable economic assumptions and redundancy/cost priorities in storage tiering.
- Channel/objective-aware codec policy controls.
- Dependency-free empirical policy-learning baseline trained from prior observations.
- Runtime JSON/environment configuration and logging helpers.
- API endpoints for simulation, noisy-read recovery, and codec policy recommendation.
- CLI `inspect` command and reconstruction diagnostics.
- API/CLI/integration/property-based tests and >=80% coverage gate.
- Expanded deterministic benchmarks for substitution, insertion, deletion, dropout, and mixed channels.
- Reproducibility metadata and separate recovery/overhead plots.
- Configuration and benchmarking documentation.

### Fixed
- GitHub CI lint failures from the initial v0.1 release.
- Direct `pytest` support for the `src/` package layout.
- Self-referential optional dependency metadata.

## 0.1.0 - 2026-10-04

### Added
- Reversible binary-to-DNA framing with protected metadata, CRC and SHA-256 verification.
- Pure-Python Reed-Solomon correction and XOR single-erasure parity groups.
- Fountain-style seeded XOR redundancy and peeling-decoder research baseline.
- Deterministic simulator for substitutions, insertions, deletions, dropout, and duplication.
- Explainable storage-tier and adaptive codec policy engines.
- Graph-style read clustering/reconstruction baseline.
- CLI, Python SDK, FastAPI surface, benchmarks, tests, CI, Docker, and research documentation.

# Changelog

## Unreleased

### Added
- Lightweight deterministic `multistart_trace_consensus()` using multiple observed-read anchors, bidirectional refinement, candidate pruning and optional known-length scoring without GPU/deep-learning dependencies.
- Leakage-controlled CNR calibration: 48 deterministic calibration clusters are excluded from the unchanged 96-cluster held-out physical-read benchmark before configuration selection.
- Paired held-out comparisons against prior OligoArk graph/alignment and iterative-trace methods in addition to the pinned external BBS baseline.

### Validation
- Calibration selected the 3-anchor, 1-round bidirectional configuration at 68/96 exact reconstructions across 5- and 10-read calibration trials versus 53/96 for the prior iterative trace baseline, within the predeclared 4× runtime budget.
- On the unchanged 96-cluster physical CNR held-out set, multi-start trace improved exact recovery from 44/96 to 53/96 versus iterative trace at five reads, and from 63/96 to 80/96 at ten reads; it also exceeded graph/alignment's 33/96 and 65/96.
- BBS remained stronger at 72–73/96 across five five-read repetitions and 93/96 in all five ten-read repetitions, so no state-of-the-art claim is made.
- Multi-start peak RSS was about 24 MB and measured held-out runtimes were 23.8 s at five reads and 61.9 s at ten reads. The final complete CPU-only GitHub workflow finished in 6 min 58 sec.

## 0.6.0 - 2026-10-05

### Added
- Iterative medoid-refined trace consensus for insertion/deletion-heavy read clusters.
- Multi-threshold explicit-graph reconstruction that emits multiple independently derived candidates while retaining normal CRC/ECC/SHA-256 acceptance gates.
- Explicit deterministic multi-read channel coverage through `copies_per_strand`, preserving legacy probabilistic duplication behavior.
- `balanced_robust` optimizer scoring with cross-seed fold-instability penalties.
- Multi-payload-content optimizer calibration and larger calibration payloads to reduce the v0.5 small-payload overfitting failure.
- Direct/graph/trace reconstruction as an optimizer policy dimension.
- Entirely new v0.6 calibration and held-out seed sets not used by the v0.5 publication study.
- Dependency-free deterministic RBF-kernel policy-learning baseline.
- Held-out learning comparison against heuristic, empirical, ridge, kernel, adaptive+fountain and measured-search strategies.
- FASTA/FASTQ/gzip physical-read reconstruction adapter with an explicit reference-oligo input.
- Verified DNA-Aeon provenance manifest for DOI `10.1038/s41467-023-36297-3`, BioProject `PRJNA855029`, and the paper-listed SRA runs.
- Strategy-partitioned publication execution: non-optimizer baselines are sharded by payload/channel/evaluation seeds, while the combined optimizer calibrates once per payload/channel before all ten untouched held-out seeds are evaluated.

### Validation
- Executed the v0.6 publication study with 1,680 untouched held-out trials, 24 optimizer calibrations and 864 calibration-candidate evaluations across 3 payload sizes, 8 regimes, 7 strategies and 10 evaluation seeds.
- Adaptive graph/alignment and iterative trace each recovered 223/240 trials (92.9%, Wilson 95% CI 89.0–95.5%); combined robust search recovered 221/240 (92.1%, 88.0–94.9%); adaptive+fountain recovered 200/240 (83.3%, 78.1–87.5%).
- In the two moderate-indel regimes, graph and trace each rescued all 38/38 direct-adaptive failures with zero paired regressions; at 8192 bytes direct recovered 0/10 and 1/10 at the two coverage levels while graph/trace recovered 10/10 in both.
- The robust optimizer was 6/6 on calibration in every payload×scenario cell, but only 19/24 cells stayed perfect on untouched seeds; the dropout-10% and substitution gaps are retained as negative generalization evidence.
- Held-out empirical, ridge and RBF-kernel policy learners did not beat the heuristic on recovery; the measured-search comparator reached 100% recovery at much higher software runtime/utility cost.
- Public DNA-Aeon provenance was verified, but no physical end-to-end decode is claimed because the repository does not bundle an explicit reproducible read-to-reference oligo mapping.

### Changed
- Version advanced to 0.6.0.
- Publication optimizer calibration uses six seeds, three independent calibration payload contents, up to 512 bytes, 36 candidates and instability-aware scoring; wall-clock runtime/retrieval weights are disabled for publication selection so seed shards choose the same deterministic winner.
- Publication experiments now include explicit medoid, alignment-graph, and iterative-trace reconstruction ablations plus two moderate-indel coverage levels.
- Physical sequencing data remains claim-limited: public reads may be evaluated only with an explicit reproducible reference mapping; no external archive format is inferred.

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

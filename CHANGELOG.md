# Changelog

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

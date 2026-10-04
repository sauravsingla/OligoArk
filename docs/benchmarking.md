# Benchmarking and reproducibility

OligoArk provides three complementary experiment paths. All results are software measurements or simulations unless explicitly linked to an external physical dataset.

## 1. Continuity benchmark

```bash
python benchmarks/run_benchmark.py
```

This retains the original fixed-vs-adaptive deterministic comparison and writes CSV/JSON, environment metadata and recovery/overhead plots.

## 2. Ablation experiment framework

Quick CI-sized sweep:

```bash
python benchmarks/run_experiments.py --profile smoke
```

Publication-oriented sweep:

```bash
python benchmarks/run_experiments.py --profile publication --output-dir publication-results
```

The publication profile varies five deterministic seeds, three payload sizes and eight channel regimes. It compares:

- `fixed`
- `adaptive`
- `adaptive_fountain`
- `adaptive_graph`
- `combined` measured optimisation

Raw and aggregated results are stored as JSON and CSV. Aggregation includes Wilson 95% intervals for verified recovery rates, mean encoded overhead, mean runtime and graph-use rate. The script also generates recovery and overhead plots.

CI runs only the smoke profile. The GitHub `Publication experiments` workflow runs the larger profile manually and uploads artifacts for 90 days.

## 3. Verified optimisation/reconstruction demonstration

```bash
python benchmarks/run_verified_demo.py
```

The demo uses separate calibration and evaluation seeds. For substitution, indel, dropout and mixed scenarios it performs:

```text
payload
  -> candidate search / adaptive optimisation
  -> hard-constrained DNA encoding
  -> selected XOR / fountain / hybrid redundancy
  -> simulated corruption
  -> explicit graph clustering
  -> alignment-aware consensus
  -> frame + erasure recovery
  -> SHA-256 verification
```

It writes `verified-demo-results.json` and intentionally reports both successes and failures rather than treating failure as a test harness error.

## Metrics

- `recovered`: exact payload equality after SHA-256-verified recovery.
- `recovery_rate`: successes / trials for the aggregate group.
- `recovery_ci95_low/high`: Wilson 95% interval for the binomial recovery proportion.
- `encoded_nucleotides`: measured length of generated software DNA strings.
- `overhead_ratio`: encoded nucleotide count divided by the ideal raw 2-bit mapping length.
- `strand_count` and `read_count`: encoded strands and simulated reads.
- `runtime_seconds`: wall-clock software runtime in the current environment.
- `graph_reconstruction_used`: whether graph/alignment fallback was needed.
- codec fields such as redundancy scheme, RS symbols and chunk size.

## Reproducibility rules

- Report the exact git commit/release, Python version, platform, experiment profile, seeds and payload sizes.
- Do not compare runtime across machines without reporting environment metadata.
- Do not convert nucleotide counts into synthesis prices unless a dated external pricing model is explicitly supplied.
- Do not describe simulated rates as observed synthesis/sequencing rates.
- Do not discard failed recoveries from aggregate statistics.
- When tuning and evaluation use different seeds, report both.

# Benchmarking and reproducibility

OligoArk provides continuity, smoke, controlled-reconstruction, held-out learning, and publication-scale experiment paths. Every result is a measured software result or a software-channel simulation unless a physical dataset is explicitly named.

## Continuity benchmark

```bash
python benchmarks/run_benchmark.py
```

This retains the historical fixed-vs-adaptive comparison and writes CSV/JSON, environment metadata, and recovery/overhead plots.

## CI smoke validation

```bash
python benchmarks/run_experiments.py --profile smoke
python benchmarks/run_learning_evaluation.py
python benchmarks/run_graph_rescue.py
python benchmarks/run_verified_demo.py
```

The smoke profile is deliberately small and is a regression gate, **not** a basis for statistical significance claims. It verifies that held-out calibration, aggregation, learning evaluation, graph-rescue evidence, plotting, and the normal recovery integrity gates remain executable.

## Publication experiment

The publication profile is defined in `oligoark.experiments.publication_profile()`. v0.5 uses:

- calibration seeds `9201,9202,9203,9204`;
- disjoint held-out evaluation seeds `2026..2033`;
- payload sizes `512, 2048, 8192` bytes;
- clean, two substitution, two indel, two dropout, and one mixed regime;
- fixed, heuristic adaptive, adaptive+hybrid redundancy, adaptive+graph/alignment, and combined measured-search strategies.

The combined optimizer is calibrated **only** on calibration seeds, frozen, and then evaluated on unseen evaluation seeds. Its budgeted candidate search is deterministic and order-independent. `balanced` search distributes coverage across redundancy/reconstruction groups; `full_grid` evaluates the complete valid grid when practical.

A local publication run is:

```bash
python benchmarks/run_experiments.py \
  --profile publication \
  --output-dir publication-results
```

GitHub Actions uses deterministic payload-size shards. Each shard runs the same seed/scenario/strategy design for one payload size, uploads its raw artifact, and the aggregation job combines all trials without dropping failures:

```text
512 B shard  ─┐
2048 B shard ├─> aggregate -> held-out learning -> graph rescue -> 90-day artifact
8192 B shard ┘
```

The aggregation artifact contains:

- `raw.json` / `raw.csv`: every held-out strategy trial;
- `summary.json` / `summary.csv`: Wilson 95% recovery intervals, mean overhead/runtime, graph-rescue rate;
- `paired-effects.json` / `paired-effects.csv`: paired differences versus the fixed strategy on identical seed/scenario/payload realizations;
- `calibration.json`: selected optimizer plans and the complete nested candidate evaluations for each payload/scenario calibration;
- `calibration-candidates.csv`: flattened candidate configurations, objective contributions, rejected-candidate reasons, and selected-winner flags;
- `metadata.json`: commit, Python/platform information, calibration/evaluation seeds, calibration payload limit, search method/seed, payload sizes and claim scope;
- plots for recovery vs configured error rate, overhead vs recovery, runtime vs recovery, strategy ablation, graph rescue, and calibration-to-held-out optimizer generalization.

## Controlled graph-rescue evidence

```bash
python benchmarks/run_graph_rescue.py
```

This intentionally starts from multiple noisy observations of a real OligoArk strand. It compares:

1. direct normal archive recovery;
2. explicit graph clustering with medoid/non-alignment consensus;
3. explicit graph clustering with alignment-aware consensus.

A rescue is counted only when direct recovery fails, reconstruction creates consensus candidate(s), and ordinary frame/ECC/CRC/SHA-256 recovery succeeds. Diagnostics include node count, candidate-pair count, retained edges, connected components, cluster sizes, consensus lengths, reconstruction runtime, and whether reconstruction changed the final result.

No expected strand is inserted into the recovery path.

## Held-out policy learning

```bash
python benchmarks/run_learning_evaluation.py \
  --experiment-dir publication-results \
  --output-dir learning-results
```

The learning pipeline converts reproducible experiment records into `PolicyObservation` data, splits evaluation seeds into disjoint training/test sets **and** holds out a disjoint subset of channel scenarios, fits both the instance-based empirical model and the deterministic ridge-regression utility model, and compares them with the deterministic heuristic and the measured combined-search result.

Reported metrics include:

- SHA-256-verified recovery rate;
- mean encoded-overhead ratio;
- mean runtime;
- policy-selection accuracy relative to the best evaluated codec candidate;
- mean utility regret;
- serialized ridge coefficient/model state for reproducibility.

If a learned model performs worse than the heuristic, the result is retained and reported.

## Metrics

- **recovered** — exact payload recovery accepted by the archive SHA-256 integrity gate.
- **recovery_rate** — successes divided by held-out trials for a group.
- **recovery_ci95_low/high** — Wilson 95% interval for the binomial recovery proportion.
- **encoded_nucleotides** — measured length of generated software DNA strings.
- **overhead_ratio** — encoded nucleotide count divided by the ideal raw 2-bit mapping length.
- **redundancy_ratio** — explicit optimizer reliability/redundancy term derived from selected XOR/fountain configuration.
- **strand_count / read_count** — generated strands and simulated reads.
- **runtime_seconds** — wall-clock software runtime in the reported environment.
- **graph_reconstruction_used** — whether reconstruction rescued a direct failure in the evaluated trial.
- **paired recovery/overhead/runtime difference** — within-realization difference from the fixed baseline.

## Optimizer objective

`OptimizationWeights` exposes recovery, overhead, redundancy, runtime, retrieval, durability, lifecycle storage cost, lifecycle retrieval cost, energy, and latency terms. Each candidate returns an `ObjectiveBreakdown` whose signed contributions reproduce the final score.

Normalized software quantities and physical/user-supplied lifecycle quantities remain separate. Lifecycle physical terms are omitted entirely when the caller does not provide lifecycle assumptions.

## Reproducibility rules

- Report the exact commit/release, Python version, platform, search method/seed, calibration seeds, evaluation seeds, and payload sizes.
- Calibration/training seeds must not overlap final evaluation/test seeds.
- Do not compare runtime across machines without environment metadata.
- Do not convert nucleotide counts into synthesis prices unless a dated external model is explicitly supplied.
- Do not describe configured software error probabilities as measured synthesis/sequencing rates.
- Do not discard failed recoveries.
- Do not call smoke-test differences statistically significant.
- Preserve negative learned-policy, optimizer, and graph-reconstruction results.

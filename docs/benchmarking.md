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

The publication profile is defined in `oligoark.experiments.publication_profile()`. v0.6 deliberately uses seed sets untouched by the v0.5 study:

- calibration seeds `9401..9406`;
- disjoint untouched evaluation seeds `31001..31010`;
- three independent calibration payload contents, up to 512 bytes;
- payload sizes `512, 2048, 8192` bytes;
- clean, 1% substitution, low indel, moderate indel at two coverage levels, two dropout regimes, and mixed noise;
- fixed, heuristic adaptive, adaptive+hybrid redundancy, adaptive+medoid graph consensus, adaptive+graph/alignment, adaptive+iterative-trace, and combined robust-search strategies;
- explicit per-strand read coverage from 1 to 8 traces depending on regime.

The combined optimizer is calibrated **only** on calibration seeds, frozen, and then evaluated on untouched evaluation seeds. Its budgeted candidate search is deterministic and order-independent. `balanced_robust` distributes coverage across redundancy/reconstruction groups and subtracts a cross-seed fold-instability penalty; `full_grid` evaluates the complete valid grid when practical. The publication profile sets runtime and retrieval objective weights to zero during candidate selection so repeated calibration on different GitHub runners cannot change the frozen winner through wall-clock jitter; runtime is still measured and reported as an outcome.

A local publication run is:

```bash
python benchmarks/run_experiments.py \
  --profile publication \
  --output-dir publication-results
```

GitHub Actions separates the combined optimizer from the non-optimizer baselines. Baselines use deterministic payload-size × scenario × evaluation-seed shards, with five untouched seeds per shard and no repeated optimizer calibration. A separate combined job calibrates once per payload/scenario and evaluates the frozen winner on all ten untouched seeds. Aggregation recombines every held-out trial without dropping failures:

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

## Executed v0.5 results

The publication workflow completed successfully on commit `180618c9f5bdc1260d00c0b60a09dd6442c1a569` and uploaded a 90-day validation artifact. The merged design contains **960 held-out trials** and **576 calibration-candidate evaluations**.

| Strategy | Successes / trials | Recovery | 95% Wilson CI | Mean overhead | Mean runtime |
| --- | ---: | ---: | ---: | ---: | ---: |
| fixed | 88 / 192 | 45.8% | 38.9–52.9% | 1.492× | 0.036 s |
| heuristic adaptive | 115 / 192 | 59.9% | 52.8–66.6% | 1.557× | 0.728 s |
| adaptive + fountain/hybrid | 128 / 192 | 66.7% | 59.7–73.0% | 2.114× | 0.959 s |
| adaptive + graph/alignment | 115 / 192 | 59.9% | 52.8–66.6% | 1.557× | 3.379 s |
| combined measured optimizer | 110 / 192 | 57.3% | 50.2–64.1% | 1.895× | 1.841 s |

Paired recovery differences versus fixed on identical held-out realizations were +14.1 percentage points for heuristic adaptive, +20.8 points for adaptive+fountain/hybrid, +14.1 points for adaptive+graph/alignment, and +11.5 points for combined search. The combined search had two paired regressions versus fixed.

Recovery by regime, aggregated over all three payload sizes and eight held-out seeds per size:

| Regime | Fixed | Adaptive | Adaptive + fountain | Adaptive + graph | Combined |
| --- | ---: | ---: | ---: | ---: | ---: |
| clean | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| substitution 0.1% | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| substitution 1% | 4.2% | 100.0% | 100.0% | 100.0% | 79.2% |
| indel low | 12.5% | 12.5% | 25.0% | 12.5% | 29.2% |
| indel moderate | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| dropout 2% | 91.7% | 91.7% | 100.0% | 91.7% | 95.8% |
| dropout 10% | 45.8% | 62.5% | 87.5% | 62.5% | 41.7% |
| mixed | 12.5% | 12.5% | 20.8% | 12.5% | 12.5% |

The combined optimizer's calibration recovery was frequently 100% while held-out recovery degraded with larger payloads in dropout, substitution and low-indel regimes. The current optimizer is therefore demonstrated as a functioning measured search, **not** as a uniformly superior policy. The four-seed, at-most-256-byte calibration budget is a documented generalization limitation.

Runtimes above are wall-clock measurements from Python 3.13.15 on the recorded GitHub-hosted Linux environment and should only be compared within this run.

## Controlled graph-rescue evidence

```bash
python benchmarks/run_graph_rescue.py
```

This intentionally starts from multiple noisy observations of a real OligoArk strand. It compares:

1. direct normal archive recovery;
2. explicit graph clustering with medoid/non-alignment consensus;
3. explicit graph clustering with alignment-aware consensus;
4. multi-threshold graph clustering with iterative trace consensus.

A rescue is counted only when direct recovery fails, reconstruction creates consensus candidate(s), and ordinary frame/ECC/CRC/SHA-256 recovery succeeds. Diagnostics include node count, candidate-pair count, retained edges, connected components, cluster sizes, consensus lengths, reconstruction runtime, and whether reconstruction changed the final result.

No expected strand is inserted into the recovery path.

## Held-out policy learning

```bash
python benchmarks/run_learning_evaluation.py \
  --experiment-dir publication-results \
  --output-dir learning-results
```

The learning pipeline converts reproducible experiment records into `PolicyObservation` data, splits evaluation seeds into disjoint training/test sets **and** holds out a disjoint subset of channel scenarios, fits instance-based empirical, deterministic ridge-regression and deterministic RBF-kernel utility models, and compares them with the heuristic, adaptive+fountain baseline and measured combined-search result.

Reported metrics include:

- SHA-256-verified recovery rate;
- mean encoded-overhead ratio;
- mean runtime;
- policy-selection accuracy relative to the best evaluated codec candidate;
- mean utility regret;
- serialized ridge coefficient/model state for reproducibility.

The executed held-out learning result used training seeds 2026–2029 on four channel regimes and test seeds 2030–2033 on four unseen regimes. Across 48 test groups: heuristic recovery was 41.7% with mean regret 0.0520; empirical recovery was 37.5% with regret 0.0740; ridge recovery was 41.7% with regret 0.0520; measured-search recovery was 39.6% with regret 0.1630. The ridge model therefore matched but did not outperform the heuristic, and the empirical/search baselines were worse on this split. These negative results are retained rather than hidden.

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


## Physical-read adapter

For an external physical dataset with an explicit reference oligo FASTA:

```bash
python benchmarks/run_physical_dataset.py \
  --manifest datasets/dna_aeon.json \
  --reads reads.fastq.gz \
  --references references.fasta \
  --output physical-reconstruction-results.json
```

The checked DNA-Aeon manifest records public SRA provenance only. The evaluator does not infer another project's archive format and does not turn reference-reconstruction accuracy into an OligoArk end-to-end decoding claim.

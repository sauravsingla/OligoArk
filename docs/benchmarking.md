# Benchmarking and reproducibility

All OligoArk benchmark results are **software measurements or seeded simulations** unless an external physical dataset is explicitly identified.

## Continuity benchmark

```bash
python benchmarks/run_benchmark.py
```

This retains the historical fixed-vs-adaptive comparison and produces JSON/CSV plus recovery and overhead plots.

## Held-out ablation experiments

CI-sized validation:

```bash
python benchmarks/run_experiments.py --profile smoke
python benchmarks/run_learning_evaluation.py
python benchmarks/run_graph_rescue.py
```

Publication profile:

```bash
python benchmarks/run_experiments.py --profile publication --output-dir publication-results
```

The preferred full publication execution is the **Publication experiments** GitHub Actions workflow. It shards by payload size, runs the same publication profile in each shard, downloads all raw shards, merges them, recomputes aggregate statistics, evaluates learned policies, runs graph-rescue validation, and uploads a 90-day validation artifact.

### Seed discipline

v0.5 separates:

- **optimizer calibration seeds** — used to choose the combined-system configuration;
- **evaluation seeds** — never seen during optimizer selection;
- **learning train/test seeds** — a further disjoint split of experiment evaluation records.

The metadata file records all seed sets, the git commit, platform, Python version, search method, search seed, candidate budget, strategies, scenarios, and payload sizes.

## Publication strategies

- `fixed`
- `adaptive` — deterministic heuristic policy
- `adaptive_fountain` — heuristic policy with hybrid XOR/fountain redundancy
- `adaptive_graph` — heuristic policy with graph/alignment fallback
- `combined` — measured candidate search frozen on calibration seeds and tested on unseen seeds

The learned empirical and ridge-regression policies are evaluated from the experiment output in a separate held-out learning analysis so training/test provenance stays explicit.

## Statistics and artifacts

Raw and aggregate outputs include:

- per-trial JSON/CSV;
- aggregate JSON/CSV;
- Wilson 95% recovery intervals;
- paired recovery/overhead/runtime differences versus the fixed baseline on identical seed/scenario/payload trials;
- encoded nucleotide overhead;
- runtime;
- graph rescue rate;
- selected codec/redundancy fields;
- calibration search score and held-out recovery;
- learned-policy recovery, overhead, runtime, selection accuracy, and regret.

Plots include recovery by strategy, overhead by strategy, recovery vs configured error rate, overhead vs recovery, runtime vs recovery, strategy ablation, graph rescue rate, and optimizer calibration-vs-held-out behavior.

## Controlled graph rescue

```bash
python benchmarks/run_graph_rescue.py
```

This benchmark constructs noisy observations from ordinary OligoArk strands and compares:

```text
direct decode
graph + medoid consensus
graph + alignment-aware consensus
```

A rescue is counted only when direct decoding fails and a reconstructed path succeeds through the normal frame/ECC/CRC/SHA-256 verification pipeline. The original strand is not injected into recovery.

## Learning evaluation

```bash
python benchmarks/run_learning_evaluation.py \
  --experiment-dir experiment-results \
  --output-dir learning-results
```

Experiment records become explicit `PolicyObservation` training rows. The evaluation reports heuristic, empirical, ridge-regression, and measured-search results on disjoint held-out seeds. The ridge model coefficients and normalization state are saved to JSON.

## Interpretation rules

- Never describe configured simulation error rates as observed sequencing/synthesis rates.
- Never drop failed recovery trials.
- Report calibration and evaluation seed sets separately.
- Do not call a smoke-test difference statistically significant.
- Runtime comparisons require the recorded machine/software environment.
- Do not translate nucleotide counts to synthesis cost without a dated, externally sourced cost model.
- Physical lifecycle cost, energy, or latency terms are used only when supplied by the caller.
- Publication conclusions should come from the completed held-out publication workflow, not CI smoke data.

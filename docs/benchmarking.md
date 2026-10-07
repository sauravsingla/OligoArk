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

## Executed v0.6 results

The v0.6 publication validation completed with **1,680 held-out trials** (3 payload sizes × 8 regimes × 7 strategies × 10 untouched evaluation seeds), 24 optimizer calibration records and **864 calibration-candidate evaluations**. The aggregate artifact records publication commit `5c1330de9e1a3bfde0ed05cd31aa30033cba0065`, Python 3.13.15 and Linux. Execution-only recovery sharding was used to finish the longest 8192-byte moderate-indel baselines without changing any scientific settings.

| Strategy | Successes / trials | Recovery | 95% Wilson CI | Mean overhead | Mean runtime |
| --- | ---: | ---: | ---: | ---: | ---: |
| fixed | 148 / 240 | 61.7% | 55.4–67.6% | 1.492× | 0.058 s |
| heuristic adaptive | 177 / 240 | 73.8% | 67.8–78.9% | 1.557× | 0.943 s |
| adaptive + fountain/hybrid | 200 / 240 | 83.3% | 78.1–87.5% | 2.114× | 1.251 s |
| adaptive + medoid | 177 / 240 | 73.8% | 67.8–78.9% | 1.557× | 29.880 s |
| adaptive + graph/alignment | 223 / 240 | **92.9%** | 89.0–95.5% | 1.557× | 36.080 s |
| adaptive + iterative trace | 223 / 240 | **92.9%** | 89.0–95.5% | 1.557× | 763.326 s |
| combined robust optimizer | 221 / 240 | **92.1%** | 88.0–94.9% | 1.710× | 33.645 s |

The two moderate-indel regimes provide the key reconstruction ablation:

| Regime | Direct adaptive | Medoid | Adaptive + fountain | Graph/alignment | Trace | Combined |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.15% insertion + 0.15% deletion, 5 copies/strand | 7/30 | 7/30 | 11/30 | **30/30** | **30/30** | **30/30** |
| same indel rate, 8 copies/strand | 15/30 | 15/30 | 23/30 | **30/30** | **30/30** | **30/30** |

Graph/alignment and iterative trace each rescued **38/38** direct-adaptive failures across these two moderate-indel regimes with zero paired regressions. Across all eight regimes, each rescued **46/63** direct failures (5 low-indel, 23 moderate-indel, 15 high-coverage moderate-indel and 3 mixed-noise cases), again with no paired regression. `graph_reconstruction_used` counts a verified change from direct failure to integrity-checked success, not merely invocation of a reconstruction strategy.

The combined robust optimizer selected winners that were 6/6 on calibration in all 24 payload×scenario cells. On untouched held-out seeds, **19/24 cells stayed at 100%** and aggregate recovery was 221/240. Generalization was strong in indel regimes but imperfect elsewhere: 10% dropout recovered 14/30 versus 21/30 for adaptive+fountain, and 1% substitution recovered 28/30 versus 30/30. Calibration success should therefore not be reported as an unbiased estimate of held-out recovery.

Runtime is a software measurement from the recorded GitHub-hosted environment. The very high iterative-trace mean runtime is a real computational-cost result, especially at 8192-byte high coverage; it must not be hidden when comparing reconstruction methods.

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

The executed v0.6 held-out learning result used training seeds `31001–31004`, validation seeds `31005–31006`, and test seeds `31007–31010`. Training scenarios were clean, 10% dropout, 2% dropout and low-indel; test scenarios were moderate-indel, high-coverage moderate-indel, mixed and 1% substitution. Across 48 test groups: heuristic recovery was **68.8%** with mean regret 0.0378; empirical recovery was **47.9%** with regret 0.2037; linear-ridge and RBF-kernel recovery were each **47.9%** with regret 0.1923; adaptive+fountain recovered **77.1%** with regret 0.1163; measured search recovered **100%** but with 66.6 s mean software runtime and mean utility regret 1.4279. Linear and kernel selection accuracy was 79.2%, but neither learned model improved recovery over the heuristic. These negative learned-policy results and the measured-search cost trade-off are retained rather than hidden.

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

### CNR smoke example

`datasets/cnr.json` records provenance for Microsoft's Clustered Nanopore Reads dataset (MIT license): the pinned upstream commit and the SHA-256 and Git blob hashes of `Centers.txt` and `Clusters.txt`. The data is not bundled. Fetch the pinned commit, then prepare a subset:

```bash
git clone https://github.com/microsoft/clustered-nanopore-reads-dataset.git cnr
git -C cnr checkout 6938f44796185902a08381943c2895782886c5c3
python benchmarks/convert_cnr_to_physical.py \
  --centers cnr/Centers.txt \
  --clusters cnr/Clusters.txt \
  --limit 20 \
  --max-reads-per-cluster 10 \
  --reads-out cnr-reads.fasta \
  --references-out cnr-references.fasta \
  --supplied-output cnr-supplied-results.json
python benchmarks/run_physical_dataset.py \
  --manifest datasets/cnr.json \
  --reads cnr-reads.fasta \
  --references cnr-references.fasta \
  --output cnr-nearest-reference-results.json
```

The converter refuses inputs whose SHA-256 differs from the manifest. It takes the first `--limit` non-empty clusters in file order and skips empty clusters, so every selected reference has reads.

The two outputs answer different questions:

- `cnr-supplied-results.json` (`"assignment": "supplied_clusters"`) keeps CNR's own association: the reads of cluster *i* are scored only against center *i*. This is the faithful reading of the dataset.
- `cnr-nearest-reference-results.json` (`"assignment": "nearest_reference"`) exercises the generic physical-read adapter, which reassigns every read to its most similar reference by edit distance. It compares every read with every reference, so keep the subset small and raise `--max-reads` above 5,000 if needed.

On the default subset of the first 20 non-empty clusters (194 reads, under 20 seconds each), both paths recovered 17 of those 20 selected references exactly with trace consensus and left no read unassigned. Empty clusters are skipped by the converter, so this is recovery on the selected non-empty references, not a strand-level rate that counts dropout or empty clusters. Neither output is the authoritative CNR result; that is the held-out BBS comparison in [external-cnr-benchmark.md](external-cnr-benchmark.md).

Upstream limitation: the CNR README notes (8/12/2024) that the 10,000 source sequences contain long-range dependencies instead of being uniformly random, due to an error in their generation. The clustering algorithm may therefore behave unexpectedly and some recovered clusters may be malformed, which makes trace reconstruction harder.

Four reproducible external physical-read benchmarks now use the same frozen confidence-fusion reconstruction settings. On Microsoft's 110-base Clustered Nanopore Reads (CNR) dataset, OligoArk reaches **73/96 at five reads** and **93/96 at ten reads**, versus BBS **72–74/96** and **93/96**. On the independent Grass et al. Illumina dataset, the unchanged settings—apart from the mechanical 117-base target length—reach **94/96 at five reads** and **96/96 at ten reads**, versus BBS **90/96** and **95/96**. On the independent 2026 LCRC HFS-Pool-11.7K Illumina PE150 experiment, the same settings with the published 200-base target length reach **96/96 at both five and ten reads**, exactly tying BBS. On the DNAformer Pilot Illumina dataset, the same frozen settings with the mechanical 140-base target length reach **96/96 at both five and ten reads**, again exactly tying BBS; at one read both methods reach **83/96**. The LCRC association layer maps reads to the public design library before reconstruction and is kept separate from candidate scoring, while the CNR, Grass, and DNAformer benchmarks use explicit published bins. All four physical workflows are CPU-only and complete within the 10-minute benchmark budget; BBS remains substantially faster. See [external-cnr-benchmark.md](external-cnr-benchmark.md), [external-grass-benchmark.md](external-grass-benchmark.md), [external-lcrc-benchmark.md](external-lcrc-benchmark.md), and [external-dnaformer-pilot-benchmark.md](external-dnaformer-pilot-benchmark.md) for provenance, split isolation, confidence intervals, paired tests, runtime/memory, association rules, and claim boundaries.

# Research scope, novelty, validation, and claim boundaries

OligoArk is a software research framework, not a wet-lab DNA-storage system. v0.6 focuses on **untouched held-out validation** of multi-trace indel reconstruction, robust codec optimisation, policy learning, and explicit-reference physical-read evaluation while preserving the v0.5 leakage controls.

## OligoArk system contributions

1. **Measured adaptive codec optimisation.** Candidate configurations are actually encoded, corrupted by a seeded software channel, recovered through the normal decoder, SHA-256 verified, measured, and ranked.
2. **Leakage-controlled evaluation.** Optimizer calibration seeds are disjoint from evaluation seeds; the frozen winner is tested on unseen stochastic channel realizations.
3. **Balanced deterministic search.** Candidate enumeration is canonical and order-independent. Budgeted search uses deterministic balanced coverage across redundancy/reconstruction groups; full-grid mode remains available.
4. **Adaptive redundancy selection.** XOR, an independently implemented LT-style fountain baseline, and hybrid XOR+fountain strategies share the main archive/recovery path.
5. **Hard constrained encoding.** Configured GC bounds and homopolymer limits are acceptance constraints. Deterministic mask search either finds a valid sequence or fails explicitly.
6. **Explicit graph reconstruction.** Reads are nodes, qualifying similarities are weighted edges, and connected components define reconstruction clusters.
7. **Layered indel consensus.** Medoid and single-pass medoid-anchored alignment remain as ablations; v0.6 adds iterative trace consensus that repeatedly realigns to the prior consensus.
8. **Multi-threshold trace reconstruction.** The trace reconstructor builds explicit graphs at several thresholds and emits multiple independent candidate consensuses; normal frame CRC, ECC, and archive SHA-256 remain the acceptance gate.
9. **Transparent policy learning.** OligoArk includes inverse-distance empirical selection, deterministic ridge-regression utility learning, and a dependency-free Gaussian RBF-kernel utility baseline, all trained only from reproducible experiment records.
10. **Lifecycle-aware archival intelligence.** Caller-supplied cost, energy, and retrieval-latency inputs are decomposed by tier and can contribute to the measured codec objective. OligoArk supplies no fabricated physical price/energy defaults.
11. **Paired ablation experiments.** Fixed, heuristic-adaptive, redundancy-enabled, graph-enabled, and combined systems are compared on identical held-out channel realizations with raw trial preservation and Wilson recovery intervals.

## Validation design

The v0.5 publication profile uses separate seed sets:

- **Calibration seeds:** used only to select the combined optimizer configuration.
- **Evaluation seeds:** unseen by the optimizer and used to measure final recovery/overhead/runtime behavior.
- **Learning split:** publication evaluation records are further split into disjoint training and held-out test subsets for the empirical and ridge-regression policy models.

Publication artifacts preserve the exact git commit, Python/platform metadata, calibration seeds, evaluation seeds, payload sizes, search method/seed, raw trials, aggregate summaries, paired differences, policy-model state, and graph-rescue diagnostics.

## Measured v0.5 publication results

The publication workflow completed successfully on commit `180618c9f5bdc1260d00c0b60a09dd6442c1a569`. It produced 960 held-out trials, 576 retained calibration-candidate evaluations, paired effects, Wilson confidence intervals, graph-rescue diagnostics, learned-policy outputs and plots.

Overall SHA-256-verified recovery across 192 held-out trials per strategy was:

| Strategy | Recovery | 95% Wilson CI | Mean encoded overhead | Mean runtime |
| --- | ---: | ---: | ---: | ---: |
| adaptive + fountain/hybrid | 66.7% | 59.7–73.0% | 2.114× | 0.959 s |
| heuristic adaptive | 59.9% | 52.8–66.6% | 1.557× | 0.728 s |
| adaptive + graph/alignment | 59.9% | 52.8–66.6% | 1.557× | 3.379 s |
| combined measured optimizer | 57.3% | 50.2–64.1% | 1.895× | 1.841 s |
| fixed | 45.8% | 38.9–52.9% | 1.492× | 0.036 s |

The paired recovery-rate differences versus fixed, evaluated on identical simulated channel realizations, were +20.8 percentage points for adaptive+fountain/hybrid, +14.1 points for heuristic adaptive, +14.1 points for adaptive+graph/alignment, and +11.5 points for the combined optimizer. The combined optimizer had two paired regressions versus fixed and did not outperform the simpler adaptive+fountain/hybrid strategy overall.

Regime-level results are also mixed and therefore informative. All strategies recovered 100% in clean and 0.1% substitution regimes. At 1% substitution, fixed recovered 4.2%, adaptive/adaptive+fountain/adaptive+graph each recovered 100%, and the combined optimizer recovered 79.2%. At 10% dropout, adaptive+fountain recovered 87.5%, adaptive/adaptive+graph 62.5%, fixed 45.8%, and combined 41.7%. No strategy recovered the moderate-indel regime. In the low-indel regime, combined recovered 29.2%, adaptive+fountain 25.0%, and fixed/adaptive/adaptive+graph 12.5%.

The measured search therefore demonstrates a real, inspectable optimization mechanism, but its current four-seed/256-byte calibration budget can overfit stochastic conditions. This is reported as a **negative generalization result**, not tuned away after observing the held-out test set. A future optimizer study should increase calibration diversity or use sequential/uncertainty-aware search with a new untouched evaluation set.

The controlled graph-rescue experiment independently establishes capability. Both constructed cases started with direct decode failure and formed an explicit 7-node, 21-edge, one-component graph. In the substitution case, medoid and alignment consensus both produced verified recovery. In the insertion/deletion case, medoid failed while alignment-aware consensus produced a consensus of length 280 and recovered through the normal decoder with SHA-256 verification. In the larger publication sweep, graph reconstruction produced no additional direct-failure rescues, so OligoArk does not claim broad graph superiority from this release.

The policy-learning evaluation used disjoint seed and channel splits: training seeds 2026–2029 on clean, 0.1% substitution, 1% substitution and low-indel regimes; test seeds 2030–2033 on moderate-indel, 2% dropout, 10% dropout and mixed regimes. Across 48 held-out groups, the heuristic and ridge model each recovered 41.7% with mean regret 0.0520; empirical recovered 37.5% with regret 0.0740; measured search recovered 39.6% with regret 0.1630. The learned ridge baseline did not beat the heuristic.

All of these are software/simulation results. Runtime values are specific to the recorded GitHub Actions environment and must not be interpreted as physical DNA-system latency.

## v0.6 untouched evaluation design

The v0.6 simulation study uses seed sets not used in the v0.5 publication study. Optimizer calibration uses seeds `9401–9406`; final evaluation uses `31001–31010`. Calibration rotates across three independently generated payload contents and uses up to 512 bytes rather than the earlier 256-byte limit. `balanced_robust` retains deterministic, order-independent sampling but penalizes disagreement between alternating calibration-seed folds. For the sharded publication study, wall-clock runtime/retrieval terms are excluded from candidate selection so runner-speed noise cannot change the frozen winner; runtime remains a reported metric.

The channel model now separates **coverage** from extra stochastic duplication. `copies_per_strand` produces a declared number of independently corrupted traces for every surviving strand, while `duplicate_rate` preserves the legacy probability of one additional trace. This is still a software model, but it allows graph/trace reconstruction to be tested on genuine multi-read clusters rather than one or two traces.

The v0.6 publication profile contains three payload sizes, eight channel/coverage regimes, seven strategies and ten untouched evaluation seeds. It includes both moderate-indel coverage levels so the effect of trace count can be measured rather than assumed.

## Physical-data pathway

The DNA-Aeon paper (DOI `10.1038/s41467-023-36297-3`) states that sequencing data are publicly deposited under BioProject `PRJNA855029` and lists SRA runs `SRR19954693` through `SRR19954697` (with the paper's listed ordering); its Data Availability statement says no data restrictions apply. OligoArk records those identifiers in `datasets/dna_aeon.json`.

v0.6 does **not** claim an external physical archive decode. The public read accessions alone do not establish a machine-verifiable read-to-reference oligo mapping inside this repository. The physical adapter therefore requires the user to supply an explicit reference FASTA and reports medoid/alignment/iterative-trace **reference reconstruction** only. Any later end-to-end external dataset claim must document the exact reference source, preprocessing, read assignment and archive-format interpretation.

## Claim boundaries

Every reported statement should be classified as one of:

- **Measured software result** — observed directly from OligoArk execution, such as encoded nucleotide count, runtime, graph edge count, or selected policy.
- **Simulation result** — recovery under a fully specified seeded software channel.
- **External published evidence** — a result attributed to cited literature or a documented dataset.
- **Hypothesis** — a proposed effect not yet established by OligoArk evidence.

OligoArk does **not** claim that its simulator reproduces a specific sequencing platform, that its mask search is a biochemical synthesis model, that its LT-style fountain baseline is DNA Fountain, that its graph/alignment consensus is HEDGES, that its normalized storage traits are measured physical properties, or that software simulations establish wet-lab performance.

## Relationship to prior work

The following works are comparison points and scientific context; OligoArk is independently implemented and does not vendor their code.

- Church, Gao & Kosuri, *Next-generation digital information storage in DNA*, Science (2012). DOI: `10.1126/science.1226355`.
- Goldman et al., *Towards practical, high-capacity, low-maintenance information storage in synthesized DNA*, Nature (2013). DOI: `10.1038/nature11875`.
- Grass et al., *Robust Chemical Preservation of Digital Information on DNA in Silica with Error-Correcting Codes*, Angewandte Chemie International Edition (2015). DOI: `10.1002/anie.201411378`.
- Erlich & Zielinski, *DNA Fountain enables a robust and efficient storage architecture*, Science (2017). DOI: `10.1126/science.aaj2038`. OligoArk's LT-style fountain code is not an implementation of DNA Fountain.
- Bornholt et al., *A DNA-Based Archival Storage System*, ASPLOS (2016). DOI: `10.1145/2872362.2872397`. This is prior archival-system architecture work, not OligoArk tiering.
- Organick et al., *Random access in large-scale DNA data storage*, Nature Biotechnology (2018). DOI: `10.1038/nbt.4079`.
- Ceze, Nivala & Strauss, *Molecular digital data storage using DNA*, Nature Reviews Genetics (2019). DOI: `10.1038/s41576-019-0125-3`. This review frames DNA as an archival medium and discusses systems challenges.
- Press et al., *HEDGES error-correcting code for DNA storage corrects indels and allows sequence constraints*, PNAS (2020). DOI: `10.1073/pnas.2004821117`. OligoArk does not implement HEDGES.
- Matange, Tuck & Keung, *DNA stability: a central design consideration for DNA data storage systems*, Nature Communications (2021). DOI: `10.1038/s41467-021-21587-5`. OligoArk does not translate this literature into built-in lifetime/energy numbers.
- Welzel et al., *DNA-Aeon provides flexible arithmetic coding for constraint adherence and error correction in DNA storage*, Nature Communications (2023). DOI: `10.1038/s41467-023-36297-3`.
- Sabary et al., *Reconstruction algorithms for DNA-storage systems*, Scientific Reports (2024). DOI: `10.1038/s41598-024-51730-3`. This formalizes reconstruction from multiple traces with insertion, deletion, and substitution errors.
- Schwarz & Freisleben, *Data recovery methods for DNA storage based on fountain codes*, Computational and Structural Biotechnology Journal (2024). DOI: `10.1016/j.csbj.2024.04.048`.
- Schwarz & Freisleben, *Optimizing fountain codes for DNA data storage*, Computational and Structural Biotechnology Journal (2024). DOI: `10.1016/j.csbj.2024.10.038`.
- *Robust multi-read reconstruction from noisy clusters using deep neural network for DNA storage* (RobuSeqNet), Computational and Structural Biotechnology Journal (2024). DOI: `10.1016/j.csbj.2024.02.019`. OligoArk uses no proprietary or pretrained neural model.
- *ReLume: Enhancing DNA storage data reconstruction with flow network and graph partitioning*, Methods (2025). DOI: `10.1016/j.ymeth.2025.03.022`. ReLume motivates graph-based reconstruction comparisons; OligoArk does not claim equivalent methodology or performance.

These sources support the importance of constrained coding, fountain/rateless recovery, indel-aware reconstruction, multi-read consensus, and graph-based reconstruction. They do not establish OligoArk's own performance.

## Current research questions

- Does balanced search-based policy selection generalize to unseen seeds better than fixed or heuristic policies?
- What recovery probability is gained per additional encoded nucleotide?
- When does XOR, fountain, or hybrid redundancy provide the best held-out recovery/overhead trade-off?
- Can alignment graph consensus rescue direct-decoding failures, and when does medoid consensus fail?
- Does the ridge-regression policy model improve held-out selection regret relative to heuristic and empirical baselines?
- How stable are findings across payload sizes and error regimes?
- How sensitive are archival-tier recommendations and codec decisions to caller-supplied lifecycle assumptions?
- Can future learned edge scorers outperform deterministic Levenshtein graph construction while preserving the same integrity gate?

## Remaining physical-validation gap

The release remains simulation-first for end-to-end OligoArk archive recovery. v0.6 adds a reproducible external-read adapter and verified public DNA-Aeon provenance, but the exact external read-to-reference mapping remains an explicit required input. Until that mapping and method-equivalent archive protocol are documented, no OligoArk result should be described as an end-to-end physical DNA-storage decode.

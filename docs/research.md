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

The executed v0.6 publication study uses seed sets untouched by v0.5. Optimizer calibration uses seeds `9401–9406`; final evaluation uses `31001–31010`. Calibration rotates across three independently generated payload contents and uses up to 512 bytes. `balanced_robust` uses deterministic candidate sampling plus a cross-seed instability penalty. Publication selection disables wall-clock runtime/retrieval terms so runner-speed jitter cannot change the frozen winner; runtime remains a reported outcome.

The channel model separates deterministic multi-read coverage from optional extra duplication. `copies_per_strand` emits independently corrupted traces for every surviving strand while `duplicate_rate` retains the legacy probability of one additional trace. The publication profile spans three payload sizes, eight channel/coverage regimes, seven strategies and ten untouched evaluation seeds: **1,680 held-out trials**. The aggregate also retains 24 optimizer calibration records and 864 candidate evaluations. Calibration/evaluation seed sets are disjoint and failed trials are never removed.

## Measured v0.6 publication results

The final v0.6 validation artifact was produced from the unchanged scientific implementation merged in PR #11 and records publication commit `5c1330de9e1a3bfde0ed05cd31aa30033cba0065`, Python 3.13.15 and Linux. The long 8192-byte moderate-indel baseline shards were completed through execution-only recovery sharding; scientific profile, seeds, channel settings and strategies were unchanged.

Overall SHA-256-verified recovery across 240 held-out trials per strategy was:

| Strategy | Recovery | 95% Wilson CI | Mean encoded overhead | Mean runtime |
| --- | ---: | ---: | ---: | ---: |
| adaptive + graph/alignment | **92.9% (223/240)** | 89.0–95.5% | 1.557× | 36.080 s |
| adaptive + iterative trace | **92.9% (223/240)** | 89.0–95.5% | 1.557× | 763.326 s |
| combined robust optimizer | **92.1% (221/240)** | 88.0–94.9% | 1.710× | 33.645 s |
| adaptive + fountain/hybrid | 83.3% (200/240) | 78.1–87.5% | 2.114× | 1.251 s |
| heuristic adaptive | 73.8% (177/240) | 67.8–78.9% | 1.557× | 0.943 s |
| adaptive + medoid | 73.8% (177/240) | 67.8–78.9% | 1.557× | 29.880 s |
| fixed | 61.7% (148/240) | 55.4–67.6% | 1.492× | 0.058 s |

Runtime values are GitHub-hosted software wall-clock measurements, not physical DNA-system latency and not cross-machine performance claims. In particular, iterative trace reconstruction is substantially more expensive than alignment reconstruction in the 8192-byte high-coverage regime.

The main v0.6 scientific result is the moderate-indel reconstruction improvement on untouched held-out trials. Aggregated over all three payload sizes:

| Regime | Direct adaptive | Medoid | Adaptive + fountain | Graph/alignment | Iterative trace | Combined |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| moderate indel, 5 copies/strand | 7/30 | 7/30 | 11/30 | **30/30** | **30/30** | **30/30** |
| moderate indel, 8 copies/strand | 15/30 | 15/30 | 23/30 | **30/30** | **30/30** | **30/30** |

At 8192 bytes specifically, direct adaptive recovered 0/10 in the five-copy regime and 1/10 in the eight-copy regime, while graph/alignment and iterative trace each recovered 10/10 in both. Paired on identical realizations, graph and trace each rescued **38/38 direct failures across the two moderate-indel regimes with zero regressions**. Across all eight publication regimes they each rescued **46/63** direct-adaptive failures: 5 low-indel, 23 moderate-indel, 15 high-coverage moderate-indel and 3 mixed-noise rescues, with no paired regression.

The combined optimizer also improved substantially over v0.5, but its calibration is not perfectly predictive. Every selected winner achieved 6/6 calibration recovery in all 24 payload×scenario cells. On untouched held-out seeds, 19/24 cells remained perfect and aggregate recovery was 221/240. It achieved 30/30 in both moderate-indel regimes, 30/30 in low-indel and mixed noise, but only 14/30 at 10% dropout versus 21/30 for adaptive+fountain, and 28/30 at 1% substitution versus 30/30 for adaptive+fountain. The calibration-to-held-out gap is therefore retained as a negative generalization result.

The controlled graph-rescue experiment independently remains positive. Both constructed cases start from direct failure. Alignment and iterative trace recover both with the normal integrity gate; medoid recovers the substitution case but fails the insertion/deletion case. Unlike v0.5, the v0.6 publication sweep now shows broad verified rescues outside the controlled examples.

The held-out policy-learning evaluation uses training seeds `31001–31004`, validation seeds `31005–31006`, and test seeds `31007–31010`. Training scenarios are clean, 10% dropout, 2% dropout and low-indel; test scenarios are moderate-indel, high-coverage moderate-indel, mixed and 1% substitution. Across 48 test groups: heuristic recovery was 68.8% with mean regret 0.0378; empirical recovery 47.9% with regret 0.2037; linear-ridge recovery 47.9% with regret 0.1923; RBF-kernel recovery 47.9% with regret 0.1923; adaptive+fountain recovery 77.1% with regret 0.1163; and measured search recovery 100% with mean regret 1.4279 because its measured runtime/utility cost was much higher. Linear and kernel selection accuracy was 79.2%, but neither learned model improved held-out recovery over the heuristic. These negative learning results are retained.

## Historical v0.5 baseline

The v0.5 workflow evaluated 960 held-out trials on commit `180618c9f5bdc1260d00c0b60a09dd6442c1a569`. Its best overall strategy was adaptive+fountain/hybrid at 66.7% (128/192; Wilson 95% CI 59.7–73.0%); combined measured search recovered 57.3% (110/192). No v0.5 strategy recovered the moderate-indel regime. The controlled alignment graph rescued 2/2 constructed direct-failure cases, but the broader v0.5 sweep showed no additional direct-failure rescues. v0.6 therefore addresses the specific moderate-indel and broad-rescue gaps while retaining optimizer, learned-policy, runtime and physical-validation limitations.

All results in this section are software/simulation measurements. They do not establish sequencing-platform fidelity, synthesis performance, wet-lab recovery, or physical-media economics.

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

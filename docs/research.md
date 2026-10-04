# Research scope, novelty, validation, and claim boundaries

OligoArk is a software research framework, not a wet-lab DNA-storage system. v0.5 focuses on **held-out validation** of constrained coding, redundancy selection, graph/alignment reconstruction, policy learning, and lifecycle-aware archival planning.

## OligoArk system contributions

1. **Measured adaptive codec optimisation.** Candidate configurations are actually encoded, corrupted by a seeded software channel, recovered through the normal decoder, SHA-256 verified, measured, and ranked.
2. **Leakage-controlled evaluation.** Optimizer calibration seeds are disjoint from evaluation seeds; the frozen winner is tested on unseen stochastic channel realizations.
3. **Balanced deterministic search.** Candidate enumeration is canonical and order-independent. Budgeted search uses deterministic balanced coverage across redundancy/reconstruction groups; full-grid mode remains available.
4. **Adaptive redundancy selection.** XOR, an independently implemented LT-style fountain baseline, and hybrid XOR+fountain strategies share the main archive/recovery path.
5. **Hard constrained encoding.** Configured GC bounds and homopolymer limits are acceptance constraints. Deterministic mask search either finds a valid sequence or fails explicitly.
6. **Explicit graph reconstruction.** Reads are nodes, qualifying similarities are weighted edges, and connected components define reconstruction clusters.
7. **Alignment-aware consensus.** Each cluster uses a medoid-anchored global alignment so insertion/deletion evidence can affect consensus.
8. **Controlled graph-rescue validation.** Dedicated experiments compare direct decoding, medoid graph consensus, and alignment graph consensus and require ordinary frame/ECC/CRC/SHA-256 recovery for success.
9. **Transparent policy learning.** OligoArk includes inverse-distance empirical selection and deterministic ridge-regression utility learning, trained only from reproducible experiment records.
10. **Lifecycle-aware archival intelligence.** Caller-supplied cost, energy, and retrieval-latency inputs are decomposed by tier and can contribute to the measured codec objective. OligoArk supplies no fabricated physical price/energy defaults.
11. **Paired ablation experiments.** Fixed, heuristic-adaptive, redundancy-enabled, graph-enabled, and combined systems are compared on identical held-out channel realizations with raw trial preservation and Wilson recovery intervals.

## Validation design

The v0.5 publication profile uses separate seed sets:

- **Calibration seeds:** used only to select the combined optimizer configuration.
- **Evaluation seeds:** unseen by the optimizer and used to measure final recovery/overhead/runtime behavior.
- **Learning split:** publication evaluation records are further split into disjoint training and held-out test subsets for the empirical and ridge-regression policy models.

Publication artifacts preserve the exact git commit, Python/platform metadata, calibration seeds, evaluation seeds, payload sizes, search method/seed, raw trials, aggregate summaries, paired differences, policy-model state, and graph-rescue diagnostics.

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

The v0.5 release remains simulation-first. Public wet-lab/NGS datasets should be added only when licensing, provenance, preprocessing, and exact evaluation protocol can be documented reproducibly. Until then, no simulation result should be described as physical DNA-storage performance.

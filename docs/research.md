# Research scope, novelty, and claim boundaries

OligoArk is a software research framework, not a wet-lab DNA storage system. v0.4 focuses on **joint, reproducible evaluation** of constrained coding, redundancy selection, graph/alignment reconstruction, policy optimisation and lifecycle-aware archival planning.

## What is implemented as OligoArk research machinery

1. **Measured adaptive codec optimisation.** Candidate configurations are actually encoded, corrupted by a seeded software channel, recovered, SHA-256 verified, measured and ranked. This is distinct from the fast heuristic policy.
2. **Adaptive redundancy selection.** XOR, an independently implemented LT-style fountain baseline, and hybrid XOR+fountain strategies are integrated in the production archive/recovery path.
3. **Hard constrained encoding.** User-configured GC bounds and homopolymer limits are acceptance constraints, not only soft scores. Deterministic mask search either finds a valid strand or raises an explicit constraint failure.
4. **Explicit graph reconstruction.** Reads form nodes; qualifying pairwise similarity scores form weighted edges; connected components define clusters.
5. **Alignment-aware consensus.** Each graph cluster uses a medoid-anchored global alignment so insertion/deletion evidence can affect consensus rather than being discarded by same-length voting.
6. **Two learned-policy baselines.** OligoArk retains inverse-distance empirical selection and adds deterministic ridge-regression utility learning from caller-supplied experiment observations.
7. **Lifecycle-aware archival intelligence.** When the caller supplies per-tier values, the system exposes retention-horizon storage cost, retrieval cost, idle/retrieval energy and retrieval latency estimates. No vendor or future-DNA prices are built in.
8. **Ablation experiments.** Fixed, adaptive, redundancy-enabled, graph-enabled and combined systems are compared under identical seeded scenarios. Publication sweeps vary seed, payload size and error regime and report Wilson recovery intervals.

## Claim boundaries

Every result must be labeled as one of:

- **Measured software result** — runtime, encoded nucleotide count, graph edges, etc. observed from OligoArk code.
- **Simulation result** — recovery under an explicit stochastic software channel and deterministic seed.
- **External published result** — a statement attributed to cited literature or data.
- **Hypothesis** — a proposed effect that has not been validated.

OligoArk does **not** claim that its mask search is a biophysical synthesis model, that its LT-style fountain baseline is DNA Fountain, that its alignment consensus is state of the art, that its normalized tier traits are measured physical properties, or that simulated results establish wet-lab performance.

## Relationship to prior work

- Church, Gao & Kosuri demonstrated early large-scale digital information encoding in synthetic DNA. DOI: `10.1126/science.1226355`.
- Goldman et al. demonstrated a practical DNA-storage encoding architecture with redundancy. DOI: `10.1038/nature11875`.
- Grass et al. combined DNA preservation with error-correcting codes. DOI: `10.1002/anie.201411378`.
- Erlich & Zielinski introduced **DNA Fountain**, using fountain coding with screening of sequence constraints. DOI: `10.1126/science.aaj2038`. OligoArk's fountain module is **not** an implementation of DNA Fountain.
- Organick et al. demonstrated random access in large-scale DNA storage. DOI: `10.1038/nbt.4079`.
- Press et al. introduced **HEDGES**, an indel-capable code that also supports sequence constraints. DOI: `10.1073/pnas.2004821117`. OligoArk does not implement or claim equivalence to HEDGES.
- Welzel et al. presented **DNA-Aeon**, supporting user-defined GC/homopolymer constraints and correction of substitutions, indels and strand loss. DOI: `10.1038/s41467-023-36297-3`.
- A position-limited constrained DNA coding study explicitly addressed GC balance, homopolymer avoidance and multiple error correction. DOI: `10.1093/bib/bbac484`.
- Sabary et al. studied DNA reconstruction from multiple noisy traces with insertion, deletion and substitution errors and dynamic-programming reconstruction algorithms. DOI: `10.1038/s41598-024-51730-3`.
- Schwarz & Freisleben studied optimisation of fountain codes specifically for DNA-storage channel properties. DOI: `10.1016/j.csbj.2024.10.038`.
- RobuSeqNet explored attention/deep-network reconstruction of noisy and contaminated read clusters. DOI: `10.1016/j.csbj.2024.02.019`.
- ReLume explored flow networks and graph partitioning for large-scale DNA-storage reconstruction. DOI: `10.1016/j.ymeth.2025.03.022`.

These citations motivate comparison points; they do not establish that OligoArk matches the published methods.

## Current hypotheses

- Search-based policy selection can improve verified recovery per encoded nucleotide over a fixed policy in some simulated channels.
- Hybrid erasure protection can improve dropout tolerance relative to XOR-only or fountain-only configurations at additional overhead.
- Explicit graph clustering plus alignment consensus can improve recovery from duplicate indel-corrupted reads relative to direct decoding and medoid-only consensus.
- Learned utility models can rank candidate policies across nearby channel conditions better than simple nearest-neighbor selection when trained on adequate experiment grids.
- Joint lifecycle/codec planning may expose trade-offs that are hidden when storage tier and DNA codec are selected independently.

## Evaluation questions

- Under which channel regimes does each adaptive component improve SHA-256-verified recovery?
- How much encoded-nucleotide overhead is required for each gain?
- Does alignment consensus improve insertion/deletion recovery without harming substitution-only clusters?
- Which graph thresholds best balance fragmented clusters against cross-strand contamination?
- Does fountain/hybrid selection improve dropout recovery, and at what redundancy ratio?
- Does the measured optimizer generalize to evaluation seeds that were not used for calibration?
- How stable are conclusions across payload sizes and deterministic seeds?
- How sensitive are tier recommendations to user-supplied lifecycle cost/energy/latency assumptions?
- Can future learned graph scorers outperform the deterministic Levenshtein graph while preserving the same SHA-256 verification gate?

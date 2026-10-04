# Research scope, novelty, and claim boundaries

OligoArk is a software research framework, not a wet-lab DNA storage system. Its intended contribution is the **integration and joint evaluation** of archival intelligence, adaptive coding, and reconstruction layers while preserving deterministic baselines.

## Research hypotheses

1. **Heterogeneous archival intelligence:** a transparent workload model can make future-DNA tier decisions explainable when retention, access, mutability, durability, energy, redundancy, and externally supplied economics are considered together.
2. **Channel-aware policy adaptation:** changing chunk size and redundancy from a channel profile can improve recovery per encoded nucleotide relative to a fixed policy in some regimes.
3. **Graph-assisted reconstruction:** read-similarity graphs and consensus reconstruction can create valid frames that direct decoding cannot recover, especially when multiple noisy observations of the same strand exist. OligoArk exposes typed scorer/reconstructor interfaces so learned graph methods can be evaluated against the same deterministic baseline and checksum gate.
4. **Empirical policy learning:** policy choices learned from prior reproducible experiments can complement deterministic rules without requiring a proprietary pretrained model.
5. **Integrated archival intelligence:** workload priorities and channel conditions can be composed into one explainable tier + codec plan rather than optimized in isolation.

These are testable hypotheses. OligoArk does **not** claim that the current heuristic is globally optimal, that the graph baseline outperforms published decoders, or that DNA is currently cheaper/faster than established storage media.

## Evidence labels

Every result contributed to OligoArk should be identified as one of:

- **Measured software result** — runtime, encoded nucleotide count, memory use, etc., measured from OligoArk code.
- **Simulation result** — recovery under a declared stochastic software channel and deterministic seed.
- **External published result** — a claim attributed to a cited paper or dataset.
- **Hypothesis** — an unvalidated research proposition.

Wet-lab language must never be attached to a software simulation result.

## Evaluation questions

- Under which simulated error regimes does adaptive policy selection improve recovery probability per encoded nucleotide?
- What is the trade-off among redundancy, recovery success, runtime, and graph-reconstruction cost?
- Which q-gram/graph thresholds preserve true duplicate-read neighborhoods without joining unrelated strands?
- Can alignment-aware graph methods or learned `EdgeScorer`/`ReadReconstructor` implementations improve insertion/deletion recovery while preserving checksum-verified correctness?
- Can an empirical policy model generalize across channel profiles better than deterministic rules while remaining explainable?
- At what workload boundaries would a future DNA tier become favorable under user-supplied cost/latency/durability assumptions?
- Does joint workload/channel planning outperform independently chosen tier and codec policies on recovery-per-overhead and lifecycle-cost objectives?

## Selected foundational references

- Church, G. M., Gao, Y., & Kosuri, S. (2012). *Next-generation digital information storage in DNA*. Science 337(6102), 1628. DOI: `10.1126/science.1226355`.
- Goldman, N. et al. (2013). *Towards practical, high-capacity, low-maintenance information storage in synthesized DNA*. Nature 494, 77–80. DOI: `10.1038/nature11875`.
- Grass, R. N. et al. (2015). *Robust Chemical Preservation of Digital Information on DNA in Silica with Error-Correcting Codes*. Angewandte Chemie International Edition 54, 2552–2555. DOI: `10.1002/anie.201411378`.
- Erlich, Y., & Zielinski, D. (2017). *DNA Fountain enables a robust and efficient storage architecture*. Science 355(6328), 950–954. DOI: `10.1126/science.aaj2038`.
- Yazdi, S. M. H. T., Gabrys, R., & Milenkovic, O. (2017). *Portable and Error-Free DNA-Based Data Storage*. Scientific Reports 7, 5011. DOI: `10.1038/s41598-017-05188-1`.
- Organick, L. et al. (2018). *Random access in large-scale DNA data storage*. Nature Biotechnology 36, 242–248. DOI: `10.1038/nbt.4079`.

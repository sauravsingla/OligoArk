# Research scope and novelty

OligoArk is a software research framework, not a wet-lab DNA storage system. Its intended novelty is the **integration and joint evaluation** of layers often studied separately:

1. **Explainable heterogeneous archival tiering** across active and cold-storage classes, with an explicitly experimental `dna_future` tier.
2. **Channel-aware codec policy selection** that changes chunk size, Reed-Solomon strength, parity grouping, and sequence-mask behavior from an error profile.
3. **Graph-assisted reconstruction interfaces** with a deterministic similarity-graph baseline that can later be replaced by a GNN without changing the codec.
4. **Fountain-style redundancy baseline** using seeded overlapping XOR symbols and peeling decode, independently implemented and not represented as the published DNA Fountain algorithm.

These are research hypotheses. OligoArk does not claim that its current policy heuristic is globally optimal, that its graph baseline is state of the art, or that DNA is currently cheaper/faster than established media.

## Evaluation questions

- Under which simulated error regimes does adaptive policy selection improve recovery probability per encoded nucleotide?
- What is the trade-off between redundancy, recovery success, and runtime?
- Can graph-derived clustering or learned edge scoring improve consensus recovery from noisy/duplicate reads?
- At what workload boundaries would a future DNA tier become favorable under externally supplied assumptions?

## References

- Church, G. M., Gao, Y., & Kosuri, S. (2012). *Next-generation digital information storage in DNA*. Science 337(6102), 1628. DOI: `10.1126/science.1226355`.
- Erlich, Y., & Zielinski, D. (2017). *DNA Fountain enables a robust and efficient storage architecture*. Science 355(6328), 950–954. DOI: `10.1126/science.aaj2038`.
- Organick, L. et al. (2018). *Random access in large-scale DNA data storage*. Nature Biotechnology 36, 242–248. DOI: `10.1038/nbt.4079`.

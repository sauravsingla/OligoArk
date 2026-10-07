# OligoArk 🧬

**DNA archival storage for large files with exact SHA-256 recovery.**

OligoArk is a research project for encoding digital data into DNA-like sequences and recovering
the original bytes exactly. A run counts as successful only when the recovered payload matches
the source SHA-256 digest.

## Latest validated results

### 10 MiB matched codec comparison

The latest completed matched comparison uses the same deterministic **10 MiB** payload,
maximum **152 nt** strands, approximately **25% redundancy**, the same fault definitions and
trial seeds, and **10 trials per condition**.

| Method | Density | Clean | 1% loss | 5% loss | Substitution | Insert/Delete | Mixed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **OligoArk compact-v3 hybrid** | **1.263 bits/nt** | **10/10** | **10/10** | **10/10** | **9/10** | timeout 10/10 | **9/10** |
| DNA Fountain clean-room | **1.347 bits/nt** | 10/10 | **10/10** | **10/10** | 0/10 | 0/10 | 0/10 |
| Goldman-style + XOR | 0.515 bits/nt | 10/10 | 0/10 | 0/10 | 0/10 | 0/10 | 0/10 |

OligoArk compact-v3 uses **1.263 bits/nt** at **25.0007% measured strand redundancy** while
staying within the **152-nt** strand limit.

For the 10 MiB encode, OligoArk v3 took **89.34 s**, compared with **399.55 s** for the
clean-room DNA Fountain implementation. Peak RSS observed across v3 conditions was about
**443 MiB**, versus about **1.89 GiB** for DNA Fountain.

The current 152-nt result is:

- **10/10 exact recovery** in clean conditions;
- **10/10 exact recovery** at **1% strand loss**;
- **10/10 exact recovery** at **5% strand loss**;
- **9/10 exact recovery** under the tested substitution condition;
- **9/10 exact recovery** under the tested mixed-fault condition;
- indel-only recovery remains unresolved at this scale because all 10 trials exceeded the
  configured **75-second per-trial deadline**.

A matched **100 MiB** codec comparison has not yet been claimed as complete.

### 1 GiB scalable storage

OligoArk also has a completed **1 GiB** bounded-memory storage result.

| Test | Result | Peak RAM |
| --- | ---: | ---: |
| Clean | ✅ PASS | 43.85 MiB |
| 1% strand loss | ✅ PASS | 44.08 MiB |
| 5% strand loss | ✅ PASS | 43.63 MiB |

At 1% loss, OligoArk recovered **45,338 / 45,338** lost strands.

At 5% loss, it recovered **226,705 / 226,705** lost strands.

The 1 GiB run used **1.646 bits per nucleotide** and **12.5% redundancy**.

## Latest physical-read reconstruction evidence

These results use published sequencing datasets created by other projects. They evaluate
OligoArk reconstruction, not an end-to-end physical OligoArk storage experiment.

| Dataset | Reads/strand | OligoArk | Pinned BBS |
| --- | ---: | ---: | ---: |
| Microsoft CNR | 5 | **73/96** | 72–74/96 |
| Microsoft CNR | 10 | **93/96** | 93/96 |
| Grass et al. | 5 | **94/96** | 90/96 |
| Grass et al. | 10 | **96/96** | 95/96 |
| LCRC HFS-11.7K | 5 | **96/96** | 96/96 |
| LCRC HFS-11.7K | 10 | **96/96** | 96/96 |
| DNAformer Pilot | 1 | **83/96** | 83/96 |
| DNAformer Pilot | 5 | **96/96** | 96/96 |
| DNAformer Pilot | 10 | **96/96** | 96/96 |

## Wet-lab status

The repository includes synthesis-ready FASTA/CSV output, experiment metadata, sequencing
input requirements, reconstruction tooling, read-depth planning, and exact SHA-256
verification.

**Current status: prepared, not physically executed.**

An end-to-end physical result still requires:

**OligoArk encode → DNA synthesis → storage → DNA sequencing → reconstruction → exact SHA-256
match**

## Current limitations

The main open items are:

- improve **152-nt indel recovery** at 10 MiB scale;
- complete the matched **100 MiB** codec comparison;
- execute a true end-to-end OligoArk wet-lab experiment.

OligoArk should not yet be described as universally better than every DNA-storage codec or as
a fully proven physical DNA-storage system.

More details: [scalable storage](docs/scalable-storage.md) ·
[codec comparison](docs/dna-fountain-baseline.md) ·
[wet-lab protocol](docs/wet-lab-validation.md)

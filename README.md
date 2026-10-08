# OligoArk 🧬

**Explore how digital files can be encoded into DNA sequences and recovered accurately.**

OligoArk is an open-source **research project** about DNA data storage. It converts digital information into DNA-like sequences, simulates common data errors, and tests whether the original file can be recovered.

## At a glance

| Question | Answer |
| --- | --- |
| What does it do? | Encode digital data as DNA-like sequences and reconstruct it. |
| How is recovery checked? | A recovered file must match the original **SHA-256** checksum exactly. |
| What has been tested? | A matched **10 MiB** codec comparison and a separate **1 GiB** scalable-storage test. |
| What is still difficult? | Insertion/deletion errors at 10 MiB and full 100 MiB matched validation. |

## How it works

1. **Encode:** Turn digital bytes into DNA-like sequences.
2. **Simulate errors:** Test missing strands and changed, inserted, or deleted bases.
3. **Recover:** Try to reconstruct the original bytes.
4. **Verify:** Count success only when the reconstructed file's SHA-256 checksum matches exactly.

## Quick start: try OligoArk in minutes

Requires **Python 3.10 or newer** and Git. No DNA laboratory equipment or external dataset is needed.

```bash
git clone https://github.com/sauravsingla/OligoArk.git
cd OligoArk
python -m pip install -e .
python examples/quick_start.py
```

**Expected output** (the checksum is printed by the program):

```text
PASS: original data recovered exactly
Bytes: 61
SHA-256: <64-character SHA-256 digest>
Note: simulated DNA channel; no physical synthesis or sequencing.
```

This example uses a fixed input and **seed 42** to make the simulation repeatable. It encodes the message, simulates a DNA read channel, recovers the original bytes, and verifies an exact SHA-256 match. If any check fails, the script exits with an error.

[Read the example source](examples/quick_start.py) · [Explore the longer end-to-end example](examples/end_to_end.py)

## Research results

Results below are from separate experiments with different configurations; **do not compare the 1 GiB storage test directly with the 10 MiB codec benchmark**.

### Visual benchmark comparison

**Exact recovery at 10 MiB** — each block represents **1 successful trial out of 10**. A full bar means all 10 trials recovered the original file with a matching SHA-256 checksum.

| Error condition | OligoArk v3 | DNA Fountain | Goldman + XOR |
| --- | --- | --- | --- |
| No errors | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩 **10/10** | 🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦 **10/10** | 🟪🟪🟪🟪🟪🟪🟪🟪🟪🟪 **10/10** |
| 1% strand loss | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩 **10/10** | 🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦 **10/10** | — **0/10** |
| 5% strand loss | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩 **10/10** | 🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦 **10/10** | — **0/10** |
| Substitutions | 🟩🟩🟩🟩🟩🟩🟩🟩🟩 **9/10** | — **0/10** | — **0/10** |
| Insertions/deletions | — **0/10 (10 timeouts)** | — **0/10** | — **0/10** |
| Mixed errors | 🟩🟩🟩🟩🟩🟩🟩🟩🟩 **9/10** | — **0/10** | — **0/10** |

**Encoding time at 10 MiB** — shorter is better. These are the completed measurements currently reported on `main`.

| Method | Relative time (visual) | Measured time |
| --- | --- | ---: |
| OligoArk v3 | 🟩🟩🟩🟩 | **89.34 seconds** |
| DNA Fountain | 🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦 | **399.55 seconds** |

*Bars are approximate visual guides; use the numbers for exact comparisons. These are computational tests, not a completed OligoArk wet-lab experiment. Insertion/deletion timeouts count as failures. A matched 100 MiB comparison is not yet complete.*

### 10 MiB: benchmark setup and details

The matched comparison uses one deterministic **10 MiB** payload, **152-nt maximum** strands, approximately **25% redundancy**, identical error definitions and trial seeds, and **10 trials per condition**. The visual table above shows the recovery outcomes.

| Method | Density (bits/nt) | Encode time | Peak memory |
| --- | ---: | ---: | ---: |
| OligoArk compact-v3 hybrid | **1.263** | **89.34 s** | **443 MiB** |
| DNA Fountain clean-room | **1.347** | **399.55 s** | **1.89 GiB** |
| Goldman-style + XOR | **0.515** | Not reported here | Not reported here |

OligoArk v3 used **25.0007% measured strand redundancy**. Insertion/deletion-only recovery was **0/10** because all ten trials exceeded the configured **75-second per-trial deadline**; timeouts are failures, not recoveries.

The matched **100 MiB** codec comparison is still incomplete.

### 1 GiB: scalable storage test

OligoArk also has a completed **1 GiB** bounded-memory storage result.

| Test | Result | Peak RAM |
| --- | ---: | ---: |
| Clean | ✅ PASS | 43.85 MiB |
| 1% strand loss | ✅ PASS | 44.08 MiB |
| 5% strand loss | ✅ PASS | 43.63 MiB |

At 1% loss, OligoArk recovered **45,338 / 45,338** lost strands.

At 5% loss, it recovered **226,705 / 226,705** lost strands.

The 1 GiB run used **1.646 bits per nucleotide** and **12.5% redundancy**.

## Reconstruction using published sequencing data

These results evaluate OligoArk reconstruction using published sequencing datasets from other projects, not an end-to-end physical OligoArk experiment.

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

## Physical experiment status

The repository includes synthesis-ready FASTA/CSV output, experiment metadata, sequencing
input requirements, reconstruction tooling, read-depth planning, and exact SHA-256
verification.

**Current status: prepared, not physically executed.**

An end-to-end physical result still requires:

**OligoArk encode → DNA synthesis → storage → DNA sequencing → reconstruction → exact SHA-256
match**

## Known limitations and next steps

- improve **152-nt indel recovery** at 10 MiB scale;
- complete the matched **100 MiB** codec comparison;
- execute a true end-to-end OligoArk wet-lab experiment.

## Learn more

Technical documentation: [scalable storage](docs/scalable-storage.md) ·
[codec comparison](docs/dna-fountain-baseline.md) ·
[wet-lab protocol](docs/wet-lab-validation.md)

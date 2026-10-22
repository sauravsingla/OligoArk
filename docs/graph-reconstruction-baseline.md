# Reproducible graph-reconstruction baseline fixture

This document defines and documents the reproducible graph-reconstruction baseline fixture for
OligoArk (#43, referencing #18).

The baseline provides a small, CPU-only, fully deterministic synthetic-read dataset and evaluation
runner that executes in seconds on standard CI runners without requiring external datasets or
network downloads.

## Purpose and scope

The graph-reconstruction module (`oligoark.reconstruct`) reconstructs candidate consensus strands
from noisy, duplicated DNA reads containing substitution and indel errors.

This baseline establishes a frozen reference point for:
1. **Explicit graph construction and clustering**: Partitioning noisy reads into connected
   components (clusters) using pairwise normalized Levenshtein similarity.
2. **Deterministic trace consensus**: Applying alignment-based multi-trace voting across insertion
   slots and base columns to recover ground-truth strands.
3. **Consensus recovery vs end-to-end archive verification**: Clearly distinguishing
   sequence-level consensus recovery from complete SHA-256-verified payload recovery.

> [!NOTE]
> This baseline is a deterministic software simulation fixture for algorithm benchmarking.
> It does not represent or assert physical synthesis or wet-lab sequencing validation.

## Fixture specification

The fixture is generated purely in memory by `generate_synthetic_graph_baseline_fixture()` in
`tests/test_reconstruct.py` and `benchmarks/run_graph_reconstruction_baseline.py`.

### Ground truth definition
- **Payload**: 69 bytes of deterministic ASCII text (`b"OligoArk reproducible graph-reconstruction baseline fixture (#43, #18)"`).
- **Archive framing**: Encoded using standard `ArchiveConfig(chunk_size=32, rs_nsym=8, redundancy_scheme="none")`, yielding **3 distinct ground-truth strands** with lengths of 236, 236, and 132 nucleotides.
- **Inter-strand separation**: Pairwise Levenshtein similarities between the distinct ground-truth strands range between 0.42 and 0.51, well below the clustering threshold of 0.90, preventing spurious cross-cluster edges.

### Channel simulation parameters
The noisy reads are produced by passing the ground-truth strands to `simulate_channel()` with a
strictly documented `SimulationConfig`:

| Parameter | Value | Description |
| --- | --- | --- |
| `seed` | `42` | Fixed pseudorandom seed ensuring 100% deterministic channel mutations |
| `substitution_rate` | `0.015` (1.5%) | Probability of base substitution |
| `insertion_rate` | `0.005` (0.5%) | Probability of single-nucleotide insertion |
| `deletion_rate` | `0.005` (0.5%) | Probability of single-nucleotide deletion |
| `dropout_rate` | `0.0` (0.0%) | Zero strand-level loss; ensures each cluster retains trace evidence |
| `duplicate_rate` | `0.0` (0.0%) | No extra uncoordinated duplications |
| `copies_per_strand` | `5` | 5 noisy reads generated per ground-truth strand (15 reads total) |

### Ground-truth isolation
The ground truth strands and archive metadata are strictly segregated from the reconstruction
process. Only the noisy `reads: list[str]` are supplied to `GraphConsensusReconstructor.reconstruct()`.
Ground truth is used solely during post-reconstruction evaluation to measure exact recovery.

## Baseline reconstruction behavior

The baseline is evaluated with `GraphConsensusReconstructor(threshold=0.90, scorer=LevenshteinEdgeScorer(), consensus_mode="alignment")`:

1. **Graph construction**:
   - 15 read nodes.
   - 30 qualifying edges (each group of 5 reads forms an intra-cluster clique of $5 \times 4 / 2 = 10$ edges; 0 inter-cluster edges).
   - 3 connected components with cluster sizes `[5, 5, 5]`.
2. **Sequence-level consensus recovery**:
   - `alignment` mode reconstructs 3 candidate consensus reads.
   - All 3 consensus reads match the 3 ground-truth reference strands exactly (**3 / 3, 100.0%**).
3. **End-to-end archive verification**:
   - **Direct read recovery (without graph reconstruction)**: `False`. Unaligned indels and substitutions break direct frame decoding and Reed-Solomon error correction.
   - **Graph-assisted recovery**: `True`. Appending the reconstructed consensus reads allows `recover_bytes()` to decode the frames and verify the exact expected SHA-256 hash.
4. **Ablation contrast (medoid vs alignment)**:
   - When `consensus_mode="medoid"` is evaluated on the exact same graph, consensus recovery drops to **0 / 3**, and archive recovery fails (`False`). This confirms that alignment-based indel voting is necessary to repair trace length shifts.

## Reproduction commands

### Via pytest (recommended for CI and testing)
Run the focused deterministic baseline test with stdout reporting enabled (`-s`):

```bash
pytest tests/test_reconstruct.py -k test_graph_reconstruction_baseline_fixture -s
```

### Via standalone benchmark script
Run the standalone reproducible script directly with Python:

```bash
python benchmarks/run_graph_reconstruction_baseline.py
```

For JSON output:

```bash
python benchmarks/run_graph_reconstruction_baseline.py --json
```

## Expected output

When running either command, the expected human-readable output is:

```text
================== GRAPH RECONSTRUCTION BASELINE ==================
Benchmark:               graph-reconstruction-baseline
Reference:               OligoArk #43 (referencing #18)

--- Simulation Parameters ---
Random Seed:             42
Substitution Rate:       1.5%
Insertion Rate:          0.5%
Deletion Rate:           0.5%
Dropout Rate:            0.0%
Coverage (per strand):   5 traces
Ground Truth Strands:    3
Simulated Read Count:    15

--- Graph Construction & Clustering ---
Scorer:                  LevenshteinEdgeScorer (threshold=0.9)
Graph Nodes (Reads):     15
Candidate Pairs:         30
Graph Edges:             30
Connected Components:    3
Cluster Sizes:           [5, 5, 5]

--- Sequence-Level Consensus Recovery ---
Consensus Mode:          alignment
Candidate Consensuses:   3
Exact Ground-Truth Match:3 / 3 (100.0%)

--- End-to-End Archival Verification ---
Direct Read Recovery:    False (expected failure)
Graph-Assisted Recovery: True (SHA-256 verified)
Archive SHA-256 Match:   True

--- Ablation Contrast (Medoid vs Alignment) ---
Medoid Matches:          0 / 3
Medoid Archive Recovery: False (fails without alignment repair)

Reconstruction Runtime:  < 5.0 s (CPU-only)
===================================================================
```

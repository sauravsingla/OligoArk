# Benchmarking and reproducibility

Run:

```bash
python benchmarks/run_benchmark.py
```

The benchmark generates machine-readable CSV/JSON results, reproducibility metadata, and optional plots. The payload is generated deterministically from seed `2026`.

## Metrics

- `recovered`: exact payload equality after SHA-256-verified recovery.
- `encoded_nucleotides`: measured length of generated software DNA strings.
- `redundancy_ratio_vs_2bit_ideal`: generated nucleotide count divided by the theoretical nucleotide count of a raw two-bit mapping; this includes framing/ECC/parity overhead and is not biochemical efficiency.
- `strand_count`: number of encoded OligoArk strands.
- `read_count`: reads emitted by the software channel.
- `runtime_seconds`: wall-clock software runtime in the current environment.
- `graph_reconstruction_used`: whether successful recovery required the graph-consensus fallback.

## Interpretation

Do not compare runtime between machines without reporting environment metadata. Do not convert nucleotide counts into synthesis price unless a dated external pricing model is explicitly supplied. Do not describe simulated error rates as observed sequencing-platform error rates.

A failure is a useful result. In particular, the current baseline deliberately exposes insertion/deletion and some erasure regimes that need stronger alignment-aware or graph-based reconstruction.

# Matched DNA-storage codec baselines

OligoArk keeps storage-codec comparisons separate from external physical-read reconstruction.

## Baselines

### DNA Fountain clean-room reference

The DNA Fountain-style baseline follows the published Erlich & Zielinski design at a high
level: robust-soliton LT droplets, 32-bit seeds, XOR payloads, Reed-Solomon protection,
2-bit DNA mapping, and GC/homopolymer screening. It is independently implemented and is
not claimed bit-compatible with the historical GPL implementation.

### Goldman-style rotating ternary reference

The second baseline is an independently implemented Goldman-style rotating ternary codec.
Binary data are converted to trits and each trit selects one of the three DNA bases different
from the previous base, preventing homopolymers by construction. For matched dropout
experiments, a simple XOR erasure layer is added at the same nominal redundancy budget.
This is a reference implementation, not a reproduction of the original 2013 archive format.

## Fair comparison

The runner is:

```bash
python benchmarks/run_dna_fountain_baseline.py --profile full --trials 5
```

Profiles:

- `ci`: 1 KiB, clean + 5% dropout, three trials.
- `full`: 1 KiB, 64 KiB, and 1 MiB across clean, 1%/5% dropout, substitution,
  insertion/deletion, and mixed noise.
- `scale`: attempts 1 KiB, 64 KiB, 1 MiB, and 10 MiB. Each worker has an explicit
  timeout; a timeout or failure is retained in the raw results rather than removed.

Every method receives the same payload bytes, 152-nt strand ceiling, nominal redundancy
budget, channel rates, trial seeds, and SHA-256 exact-recovery definition. Results include:

- exact recovery rate and Wilson 95% confidence interval;
- logical bits/nt and encoded nucleotide count;
- measured redundancy and strand count;
- encode/decode runtime and throughput;
- peak RSS;
- timeouts and negative results.

The largest size at which all methods complete is therefore measured rather than assumed.

## External physical reads

Microsoft CNR, Grass et al., LCRC HFS, and DNAformer Pilot contain strands encoded by their
respective studies. Applying DNA Fountain or the rotating-ternary decoder to those strands
would not be a valid codec comparison.

Those datasets remain a separate **reference-strand reconstruction** benchmark against
pinned Bidirectional Beam Search (BBS). They are not evidence of end-to-end physical
OligoArk archive storage.

## Provenance

- DNA Fountain: Erlich & Zielinski, *Science* (2017), DOI `10.1126/science.aaj2038`.
- Rotating ternary concept: Goldman et al., *Nature* (2013), DOI
  `10.1038/nature11875`.

The repository records provenance but does not vendor historical third-party source code.

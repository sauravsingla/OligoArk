# DNA Fountain comparison and external-baseline policy

OligoArk includes an independent **DNA Fountain-style research baseline** for controlled
codec comparisons. It follows the published design at a high level: LT droplets, a robust
soliton degree distribution, a 32-bit seed, XOR payloads, Reed-Solomon protection, 2-bit DNA
mapping, and GC/homopolymer screening. Decoding uses ordinary LT peeling first and an exact
GF(2) elimination fallback when the ripple stalls but the accepted equation set has sufficient
rank. It is a clean-room implementation and is **not** claimed to be bit-compatible with the
historical TeamErlich implementation.

The comparison runner is:

```bash
python benchmarks/run_dna_fountain_baseline.py --profile full
```

The controlled comparison uses the same input payload, a 152-nt maximum strand target, the
same nominal redundancy budget, the same dropout rates and trial seeds, and the same
SHA-256 exact-payload recovery gate. It reports recovery rate with Wilson 95% intervals,
bits/nt, strand count, encoded nt, encode/decode runtime, and peak RSS. This is a software
codec comparison; it is not a reproduction of the 2017 wet-lab experiment.

## Why external physical benchmarks use BBS rather than forcing DNA Fountain

Published physical datasets such as Microsoft CNR, Grass et al., LCRC HFS, and DNAformer
Pilot contain strands encoded by their respective studies. Applying the DNA Fountain decoder
to strands that were not DNA-Fountain encoded would not be a fair baseline.

OligoArk therefore keeps the comparisons separated:

- **codec comparison:** OligoArk versus the clean-room DNA Fountain-style baseline on payloads
  encoded under controlled, matched software conditions;
- **physical-read reconstruction:** OligoArk confidence fusion versus pinned Bidirectional
  Beam Search (BBS) on CNR and the other supported published read datasets;
- **end-to-end archive recovery:** success only when OligoArk-generated archive bytes are
  reconstructed and the original payload SHA-256 matches.

See `external-cnr-benchmark.md`, `external-grass-benchmark.md`,
`external-lcrc-benchmark.md`, and `external-dnaformer-pilot-benchmark.md` for the
physical-read evidence and claim boundaries.

## Original DNA Fountain physical-data provenance

The 2017 DNA Fountain work is identified by DOI `10.1126/science.aaj2038`. The historical
implementation points to European Nucleotide Archive projects `PRJEB19305` and
`PRJEB19307`. OligoArk records this provenance but does not vendor third-party sequencing
data or GPL-licensed source code. A true external reproduction should pin downloaded
accessions/checksums, preserve the original codec parameters, and report it separately from
the clean-room baseline.

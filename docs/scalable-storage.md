# Scalable archival storage

OligoArk exposes two archive surfaces with deliberately different goals:

- `oligoark-archive-v1`: human-readable JSON research format;
- `oligoark-stream-v2`: compact length-prefixed binary container for bounded-memory
  file archival and scale experiments.

The streaming container reuses OligoArk frame semantics, CRC, optional Reed-Solomon
protection, XOR/fountain redundancy and the final SHA-256 integrity gate. Encoding and
recovery operate sequentially and write recovered chunks directly to logical file offsets
instead of holding the whole archive in RAM.

## Storage profiles

| Profile | Full strand | Inner RS | Purpose |
| --- | ---: | ---: | --- |
| `scale-1024` | 1024 nt max | 0 | systems scaling, throughput, memory, erasure recovery |
| `oligoark-152` | 152 nt | 8 symbols | realistic short-strand software experiments |
| `oligoark-200` | 200 nt | 8 symbols | realistic physical-design software experiments |
| `oligoark-248` | 248 nt | 8 symbols | RS-enabled noise/fault matrix |

The `scale-1024` profile intentionally relaxes biochemical constraints and inner RS so
large-file measurements reflect streaming-system behavior rather than pure-Python RS/mask
search cost. It must not be presented as a realistic physical oligo profile.

## Commands

Quick CI profile:

```bash
python benchmarks/run_storage_scale.py --profile ci
```

1 GiB scale acceptance:

```bash
python benchmarks/run_storage_scale.py --profile acceptance
```

RS-enabled 248-nt fault matrix:

```bash
python benchmarks/run_storage_scale.py --profile physical
```

Combined storage evidence:

```bash
python benchmarks/run_storage_scale.py --profile storage
```

Exhaustive research sweep:

```bash
python benchmarks/run_storage_scale.py --profile full
```

## 1 GiB acceptance gate

The acceptance profile verifies clean XOR recovery at:

```text
1 KiB -> 64 KiB -> 1 MiB -> 10 MiB -> 100 MiB -> 1 GiB
```

At 1 GiB it additionally requires exact SHA-256 recovery under:

- clean channel;
- 1% controlled data-strand dropout;
- 5% controlled data-strand dropout.

The milestone is achieved only if all three 1 GiB cases pass SHA-256 exactly. The benchmark
records logical bits/nt, total encoded nt, strand count, redundancy, binary archive overhead,
encode/decode throughput, wall time, peak RSS, interpreter baseline RSS, RSS growth,
constraint sampling and detailed recovery counters.

## Realistic-strand noise matrix

The `physical` profile uses `oligoark-248` with 8 Reed-Solomon symbols and runs a
64 KiB heterogeneous payload across:

- none, XOR, fountain and hybrid redundancy;
- clean;
- 1% and 5% strand dropout;
- substitutions;
- insertion/deletion errors;
- mixed noise.

Negative results are retained. RS corrects bounded symbol substitutions, but insertions and
deletions can still destroy frame alignment without multi-read reconstruction; the matrix
therefore measures the real limitation rather than converting every failure into an erasure.

## Heterogeneous payload

The deterministic fixture combines UTF-8 text, JSONL, CSV, Python source, image bytes,
random binary, compressed binary and a ZIP mixed-file archive. ZIP timestamps are fixed so
independent workers operate on identical bytes.

## Matched codec comparison

Use:

```bash
python benchmarks/run_dna_fountain_baseline.py --profile full --trials 5
```

The matched comparison uses the same payload bytes, 152-nt ceiling, nominal redundancy,
fault rates, trial seeds and SHA-256 success definition for:

- OligoArk fountain;
- clean-room DNA Fountain;
- Goldman-style rotating ternary + XOR erasure reference.

The runner preserves 95% Wilson intervals, density, redundancy, runtime, throughput,
peak RSS, failures and worker timeouts. The `scale` profile attempts through 10 MiB; the
largest size at which all methods complete is reported by the raw results rather than assumed.

## Claim boundary

These are software archive and software-channel results. Controlled strand dropout is not
measured synthesis loss, and RS-enabled software mutation tests are not sequencing evidence.
External CNR, Grass, LCRC and DNAformer experiments remain reference-strand reconstruction
benchmarks. A true physical OligoArk archive claim requires OligoArk-generated strands to be
synthesized, sequenced, reconstructed and verified against the original file SHA-256.

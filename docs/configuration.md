# Configuration

OligoArk keeps archive-format configuration separate from runtime/service settings, optimization preferences, and caller-supplied lifecycle assumptions.

## Archive configuration

`ArchiveConfig` values affect generated strands and are stored in archive metadata:

- `chunk_size`
- `rs_nsym`
- `parity_group_size`
- `adaptive_masks`
- `redundancy_scheme`: `none`, `xor`, `fountain`, or `hybrid`
- `fountain_redundancy`, `fountain_seed`, `fountain_max_degree`
- `min_gc_fraction`, `max_gc_fraction`, `max_homopolymer`
- `mask_search_limit`

GC and homopolymer values are hard software constraints. Adaptive masking deterministically searches reversible candidates; if none satisfy the configured bounds, encoding fails explicitly.

Older v0.1-v0.4 configuration mappings remain readable. Missing later fields receive the same documented defaults; the archive format remains `oligoark-archive-v1`.

## Runtime configuration

`RuntimeConfig` controls service behavior:

| Environment variable | Default | Meaning |
| --- | ---: | --- |
| `OLIGOARK_LOG_LEVEL` | `INFO` | Standard Python logging level |
| `OLIGOARK_RECONSTRUCTION_THRESHOLD` | `0.90` | Similarity threshold for graph reconstruction |
| `OLIGOARK_MAX_API_PAYLOAD_BYTES` | `10485760` | Reference API payload guard |

```bash
export OLIGOARK_LOG_LEVEL=DEBUG
export OLIGOARK_RECONSTRUCTION_THRESHOLD=0.92
export OLIGOARK_MAX_API_PAYLOAD_BYTES=5242880
uvicorn oligoark.api:app
```

The payload setting is a research-service guard, not a substitute for authentication, rate limiting, reverse-proxy limits, quotas or deployment hardening.

## Normalized economic assumptions

Storage-tier normalized economics are deliberately not built in as current prices. Supply indices from `0` (lower assumed cost) to `1` (higher assumed cost) for every tier:

```json
{
  "storage_cost_index": {
    "ssd": 0.8,
    "object_archive": 0.3,
    "tape": 0.2,
    "dna_future": 0.5
  },
  "retrieval_cost_index": {
    "ssd": 0.1,
    "object_archive": 0.4,
    "tape": 0.7,
    "dna_future": 0.9
  }
}
```

These values are user assumptions, not vendor prices asserted by OligoArk.

## Explicit lifecycle assumptions

For lifecycle estimates, provide all five fields for every tier:

```json
{
  "tiers": {
    "ssd": {
      "storage_cost_per_gb_year": 2.0,
      "retrieval_cost_per_gb": 0.0,
      "idle_energy_kwh_per_tb_year": 20.0,
      "retrieval_energy_kwh_per_gb": 0.001,
      "retrieval_latency_hours": 0.001
    },
    "object_archive": {
      "storage_cost_per_gb_year": 0.5,
      "retrieval_cost_per_gb": 0.05,
      "idle_energy_kwh_per_tb_year": 5.0,
      "retrieval_energy_kwh_per_gb": 0.002,
      "retrieval_latency_hours": 1.0
    },
    "tape": {
      "storage_cost_per_gb_year": 0.2,
      "retrieval_cost_per_gb": 0.02,
      "idle_energy_kwh_per_tb_year": 0.5,
      "retrieval_energy_kwh_per_gb": 0.005,
      "retrieval_latency_hours": 4.0
    },
    "dna_future": {
      "storage_cost_per_gb_year": 1.0,
      "retrieval_cost_per_gb": 1.0,
      "idle_energy_kwh_per_tb_year": 0.1,
      "retrieval_energy_kwh_per_gb": 0.1,
      "retrieval_latency_hours": 24.0
    }
  }
}
```

These are **illustrative schema values only**, not vendor quotes or OligoArk estimates. Replace every value with a sourced, dated assumption before interpreting lifecycle output.

Each tier response includes a `score_breakdowns` entry whose contributions reproduce that tier's final score. When lifecycle assumptions exist, the selected tier's lifecycle totals can also enter the codec optimizer. Without lifecycle inputs, physical lifecycle objective terms are omitted.

## Search-based optimization

`oligoark optimize-plan` performs real encode/simulate/recover/SHA-256-verify trials. The v0.5 search controls are:

- `--max-candidates` — evaluation budget;
- `--search-method balanced|full_grid`;
- `--search-seed` — deterministic balanced-sampling seed;
- `--seeds` — backward-compatible alias for calibration seeds;
- `--calibration-seeds` — explicit optimizer calibration/training realizations;
- `--evaluation-seeds` — optional disjoint held-out realizations;
- `--duplicate` — duplicate-read probability used by held-out evaluation;
- `--weights-json` — objective-weight override.

Example:

```bash
oligoark optimize-plan payload.bin \
  --retention-years 100 \
  --substitution 0.005 \
  --dropout 0.02 \
  --calibration-seeds 9201,9202,9203,9204 \
  --evaluation-seeds 2026,2027,2028,2029 \
  --search-method balanced \
  --search-seed 5050 \
  --max-candidates 24
```

Calibration and evaluation sets must be disjoint.

## Optimization weights

A JSON file can override any `OptimizationWeights` field:

```json
{
  "recovery": 0.40,
  "overhead": 0.15,
  "redundancy": 0.10,
  "runtime": 0.08,
  "retrieval": 0.07,
  "durability": 0.10,
  "lifecycle_storage_cost": 0.04,
  "lifecycle_retrieval_cost": 0.02,
  "lifecycle_energy": 0.02,
  "lifecycle_latency": 0.02
}
```

Weights are software research preferences, not physical measurements. Recovery/overhead/redundancy/runtime/retrieval/durability terms use normalized software quantities. Lifecycle terms are activated only when caller lifecycle assumptions supply the corresponding physical quantities.

## Reconstruction diagnostics

To compare direct recovery, graph+medoid consensus and graph+alignment consensus on the same reads:

```bash
oligoark diagnose-reconstruction \
  archive.oligoark.json reads.txt \
  --similarity-threshold 0.90
```

The response includes graph node/pair/edge/component counts, cluster sizes, consensus lengths, runtime, verification result, and whether alignment reconstruction rescued a direct failure.

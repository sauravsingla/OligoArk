# Configuration

OligoArk keeps archive-format configuration separate from runtime/service and lifecycle assumptions.

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

The GC and homopolymer settings are hard software constraints. With adaptive masking enabled, OligoArk deterministically searches reversible mask candidates. If no candidate satisfies the configured bounds, encoding fails explicitly rather than silently emitting a nonconforming strand.

Example:

```python
from oligoark import ArchiveConfig

config = ArchiveConfig(
    chunk_size=64,
    rs_nsym=16,
    redundancy_scheme="hybrid",
    fountain_redundancy=0.4,
    min_gc_fraction=0.40,
    max_gc_fraction=0.60,
    max_homopolymer=4,
    mask_search_limit=128,
)
```

## Runtime configuration

`RuntimeConfig` controls operational behavior and can be loaded from JSON or environment variables:

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

The payload setting is a research-service guard, not a substitute for reverse-proxy limits, authentication, rate limiting, or production hardening.

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

Use with `--economics-json` or the equivalent REST/Python API object. These are user assumptions, not prices asserted by OligoArk.

## Explicit lifecycle assumptions

For lifecycle estimates, provide all five values for every tier. Units are explicit and values must come from the researcher:

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

The numbers above are **illustrative schema values only**, not vendor quotes or OligoArk estimates. Replace every value with sourced, dated assumptions before interpreting lifecycle output.

CLI usage:

```bash
oligoark recommend --retention-years 100 --data-size-gb 100 --lifecycle-json lifecycle.json
```

## Search-based optimization

`oligoark optimize-plan` runs actual encode/simulate/recover/verify trials over a candidate grid. `--max-candidates` bounds search cost and `--seeds` controls deterministic trials. The Python API exposes `CodecSearchSpace` and `OptimizationWeights` for full control.

# Configuration

OligoArk keeps archive-format configuration separate from runtime/service configuration.

## Archive configuration

`ArchiveConfig` controls values that affect encoded strands:

- `chunk_size`
- `rs_nsym`
- `parity_group_size`
- `adaptive_masks`

These settings are stored inside archive metadata so decoding is self-describing.

## Runtime configuration

`RuntimeConfig` controls operational behavior and can be loaded from JSON or environment variables:

| Environment variable | Default | Meaning |
| --- | ---: | --- |
| `OLIGOARK_LOG_LEVEL` | `INFO` | Standard Python logging level |
| `OLIGOARK_RECONSTRUCTION_THRESHOLD` | `0.90` | Similarity threshold for graph reconstruction |
| `OLIGOARK_MAX_API_PAYLOAD_BYTES` | `10485760` | Reference API payload guard |

Example:

```bash
export OLIGOARK_LOG_LEVEL=DEBUG
export OLIGOARK_RECONSTRUCTION_THRESHOLD=0.92
export OLIGOARK_MAX_API_PAYLOAD_BYTES=5242880
uvicorn oligoark.api:app
```

The payload setting is a research-service guard, not a substitute for reverse-proxy limits, authentication, rate limiting, or other production hardening.

## Economic assumptions

Storage-tier economics are deliberately not built in as current prices. Supply normalized cost indices from `0` (lower assumed cost) to `1` (higher assumed cost) for every tier:

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

Use it from the CLI with `oligoark recommend --retention-years 100 --economics-json economics.json`, or pass the equivalent object to the REST `/recommend` endpoint. Values should be sourced and dated by the researcher for the scenario being studied.

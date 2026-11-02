# Engineering hardening and acceptance gates

This document records what was changed and what still needs independent validation.

## Input integrity

The codec, simulator, runtime config, storage tier inputs, and optimizer reject NaN and
infinite numeric values rather than silently accepting an invalid score or rate.
See `tests/test_numeric_guards.py` for the regression matrix. This does not turn
the reference API into a publicly exposed production service.

## Graph candidate enumeration

When the maximum possible edit budget proves that every qualifying pair must share
a q-gram, the similarity graph can iterate index-present candidate pairs in sorted
order instead of scanning every absent pair. For shorter reads or looser thresholds
it retains exhaustive pair enumeration. Exact weighted edges and candidate diagnostics
must match the prior implementation: `tests/test_graph_sparse_candidates.py`.

The inverted q-gram index still has worst-case quadratic size on high-overlap inputs.
This change is not evidence of universal linear scaling, a measured speedup, or a fix
for the 10 MiB decoder timeout. Follow up with independent profiling.

## Evidence gates before publication

1. Keep the entire cross-version CI matrix green, including strict mypy, Ruff,
   package builds, unit/property tests, security scans, and SHA-256 recovery.
2. Preserve realistic 152-nt, constraint-checked, fixed-budget comparisons with
   unchanged timeouts and independent seeds; include every failure.
3. Reproduce improvements on pinned external baselines, not only software
   approximations, and report runtime alongside exact recovery.
4. Require synthesis records, raw FASTQ, controls and exact source SHA-256 proof
   before describing any OligoArk archive as physically validated.
5. Archive benchmark data and release source immutably, with DOI and dependency lock,
   before making durable peer-review claims.

## Deployment checklist

Use authenticated access, API gateway body limits, CPU/memory quotas, request
deadlines, audit logs and rate limits if exposing the research API. Protect `main`
with required reviews/checks in GitHub settings; this cannot be achieved through
a code-only change. See `SECURITY.md` for existing threat boundaries.

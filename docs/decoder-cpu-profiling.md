# Optional decoder CPU profiling for unresolved indel/mixed timeouts

The matched codec runner now supports an **opt-in** per-trial `cProfile`
snapshot. It records decoder calls, cumulative CPU time and call counts,
including partial profiles when a decoder raises an exception or reaches its
existing deadline. No timeout, success criterion, trial seed or decoder
algorithm is changed.

Run a focused reproducible experiment on a capable runner:

```bash
python -m pip install -e ".[dev,bench]"
mkdir -p performance-results/decoder-cpu
OLIGOARK_CODEC_CPROFILE_DIR=performance-results/decoder-cpu \
  python benchmarks/run_matched_codec_scale.py \
    --profile research-10 --method oligoark-compact-v3-hybrid \
    --size 10485760 --conditions indel-low,mixed \
    --trials 3 --condition-timeout-seconds 300 \
    --output performance-results/decoder-research-10
python -m pstats performance-results/decoder-cpu/oligoark-compact-v3-hybrid-indel-low-trial-0.prof
```

**Interpretation:** Python deterministic profiler overhead can materially
increase wall time and cause additional deadline failures. Never compare
profiled wall times or recovery rates with unprofiled results as a speedup
claim. Compare hot-function cumulative CPU times and call counts first,
then validate any targeted change on identical *unprofiled* trial seeds
with the original deadline, SHA-256 exact-recovery criteria and complete
failure/timeout accounting. This is software profiling, not wet-lab evidence.

Profile files may be large; archive them as research artifacts only when
needed. Do not claim a decoder fix until measured end-to-end tests pass.

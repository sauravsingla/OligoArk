# 100 MiB codec worker progress and memory evidence (2026-11-01)

## Motivation and observed result

The 100 MiB `dna-fountain-cleanroom` job in [run 41230201569](https://github.com/sauravsingla/OligoArk/actions/runs/41230201569) hit its 2,400-second GNU timeout (exit 124). The retained artifact contained `interruption.json` and `resource-usage.txt` but no `results.json` or verified trial rows. The roughly 26 MiB maximum RSS reported by the outer `/usr/bin/time` invocation cannot be assumed to measure the Python worker and its descendants. The benchmark launches a separate worker with captured output, so a killed worker loses its unreturned JSON.

## Opt-in instrumentation

Set `OLIGOARK_CODEC_PROGRESS_FILE=/path/to/worker-progress.jsonl` before invoking `benchmarks/run_matched_codec_scale.py`. The worker appends JSON Lines records for startup, payload creation, encoding, condition and trial boundaries, completion/error, plus a heartbeat every 15 seconds. Every record contains a Unix timestamp, worker PID, and the worker process's high-water resident memory from `resource.getrusage(RUSAGE_SELF)`. The 100 MiB fountain workflow enables this automatically, placing the file alongside other always-uploaded diagnostics.

For a small reproducible smoke run:

```bash
mkdir -p matched-codec-progress-smoke
export OLIGOARK_CODEC_PROGRESS_FILE="$PWD/matched-codec-progress-smoke/worker-progress.jsonl"
python benchmarks/run_matched_codec_scale.py --profile ci \
  --method dna-fountain-cleanroom --trials 1 \
  --condition-timeout-seconds 75 \
  --output matched-codec-progress-smoke
```

This is diagnostic instrumentation, **not a codec speedup**. The heartbeat can sample memory during an encoding stall; its maximum observed value is a lower bound on actual process peak RSS between samples. It excludes descendant processes and operating-system page cache. An event indicating a trial started or exited is **not** a verified trial result: only retained `results.json`, complete trial accounting, and SHA-256 exact-recovery verdicts support such claims. When a worker is terminated, incomplete trials remain unverified and validation must fail. A missing/unwritable optional progress file emits a warning without changing codec outcomes.

The existing 10 MiB indel-low and mixed-noise timeout failures remain unresolved. The 100 MiB baseline has not demonstrated complete exact recovery, and no 1 GiB exact-recovery or physical DNA storage claim is established. Use phase and RSS evidence to prioritize further profiling before any algorithm or deadline changes.

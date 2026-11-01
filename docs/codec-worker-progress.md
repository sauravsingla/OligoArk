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

## Encoding-screening diagnostics

The clean-room DNA Fountain encoder can spend its full 2,400-second 100 MiB
watchdog in encoding before a single channel trial starts. When the worker
progress file is enabled, the baseline now records `fountain-encode-progress`
events at the start, every 65,536 screened droplet attempts, and on completion.
Each event includes accepted droplets, attempted droplets, target droplets,
acceptance rate, and sampled worker peak RSS. These are **incomplete progress
measurements**, not verified recoveries or evidence of a speedup. The callback
is opt-in, and the parity test compares the exact encoded archive with and
without instrumentation, followed by exact decoding. The existing watchdog,
negative-result accounting, screening constraints and SHA-256 gate remain
unchanged. No event is expected if the worker is killed before initialization.


## Measured 100 MiB screening failure and 10 MiB recovery limits (2026-11-01)

Evidence: [main workflow run 41270450630](https://github.com/sauravsingla/OligoArk/actions/runs/41270450630), commit
`c40ae00540d5711e7607af92c9e7a0e80a64f0e6`, artifact
`matched-codec-100mib-dna-fountain-cleanroom-41270450630`
([artifact 13069669297](https://github.com/sauravsingla/OligoArk/actions/runs/41270450630/artifacts/13069669297)).
The benchmark used a deterministic 100 MiB payload, the unmodified clean-room
Fountain screening policy, and a 2,400-second external watchdog. It exited
124 during **encoding** and never produced `results.json` or any verified
channel trial. The workflow correctly marks this job failed.

The retained `worker-progress.jsonl` contains 393 JSON Lines records:
230 `fountain-encode-progress` events, 159 worker heartbeats, and four
startup/payload/encoding phase events. At the last progress event, about
2,397.45 seconds after worker start, screening had attempted **15,007,744**
droplets and accepted **1,850,418** of **4,096,000** targeted droplets
(**12.3298%** acceptance; **45.18%** of the target). The sampled worker
high-water RSS was **1,048.605 MiB**. This is the Python worker's own
reported high-water RSS, not the outer `/usr/bin/time` wrapper's
26 MiB figure. The last sample precedes termination; actual peak RSS could
be higher. The encoder was incomplete, so encode throughput, decode
throughput, exact recovery, and a 100 MiB success rate are **not measured**.
Do not extrapolate a completed runtime or claim a speedup from these partial
progress records.

Reproduce the **artifact analysis** without re-running the expensive benchmark:

```python
import json
from pathlib import Path

path = Path("dna-fountain-cleanroom/worker-progress.jsonl")
events = [json.loads(line) for line in path.read_text().splitlines()]
progress = [row for row in events if row["event"] == "fountain-encode-progress"]
last = progress[-1]
assert last["accepted_droplets"] == 1_850_418
assert last["attempted_droplets"] == 15_007_744
assert last["target_droplets"] == 4_096_000
assert last["peak_worker_rss_mib"] == 1048.605
print(last["accepted_droplets"] / last["target_droplets"])
```

The independent [10 MiB codec validation run 41259949586](https://github.com/sauravsingla/OligoArk/actions/runs/41259949586)
(commit `55ee92499beade8e792421563e78a369bb324c72`, ten trials per
condition, 75-second per-trial deadline) retained complete negative and
positive trial outcomes:

| Codec | Clean | Dropout 1% | Dropout 5% | Substitution-low | Indel-low | Mixed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| OligoArk compact hybrid | 10/10 | 0/10 | 0/10 | 0/10 | 0/10 (10 timeouts) | 0/10 (10 timeouts) |
| OligoArk compact v3 hybrid | 10/10 | 10/10 | 10/10 | 9/10 | 0/10 (10 timeouts) | 9/10 |

These are exact-recovery successes from the retained `results.json` rows;
a timeout is a failed trial, not an omitted trial. The interrupted stacks
cluster across `dna.py`, `ecc.py`, and `framing.py`, but a single stack
per deadline is not enough to attribute CPU time to a specific function.
Before optimizing indel handling, collect repeated matched CPU/RSS profiles
on the same seeds and noise conditions, preserve the 75-second limit and
every failure, and validate exact SHA-256 recovery and randomized parity.
Do not introduce parallel reconstruction, caching, or relaxed CI timing
thresholds without controlled evidence of deterministic benefit.

This software benchmark is not a physical DNA-storage experiment, and the
1 GiB exact-recovery milestone remains unverified.

# Negative result: post-confidence-fusion CNR refinements

Date: 2026-10-06

This note preserves the negative research outcome from the next low-compute reconstruction round after the confidence-fusion CNR result. These experiments were **not** merged into the production reconstruction path because they did not improve the primary 10-read held-out metric beyond the existing 93/96 result.

## Protected development protocol

The existing 96-cluster held-out CNR benchmark (seed `20261005`) remained sealed during development. A new deterministic seed-`20261008` development split used 64 previously unused CNR clusters after excluding:

- the 96 historical held-out clusters;
- the earlier 48-cluster calibration split; and
- the prior 48-cluster seed-`20261007` development split.

The frozen confidence-fusion baseline scored **48/64 at 5 reads** and **61/64 at 10 reads** on this new development set.

## Reliability/stability fusion

Read-reliability weighting, reverse-consensus stability, alignment-confidence scoring and small local candidate search were evaluated first.

Result: the tested variants remained at **61/64 at 10 reads**, provided no 10-read rescues, and were materially slower than confidence fusion. They were rejected before any held-out evaluation.

## Homopolymer-balance repair

A bounded homopolymer gap-migration repair was then evaluated. The selected low-cost configuration improved the new development split to:

- **49/64 at 5 reads**, versus 48/64 baseline;
- **62/64 at 10 reads**, versus 61/64 baseline;
- **1 rescue and 0 regressions** at 10 reads.

The configuration was frozen before held-out validation.

## Held-out result

A first full validation job exceeded the 10-minute job budget while rerunning unnecessary historical methods and BBS repetitions. No parameter changes were made from that run. A lean technical retry used the same frozen configuration, evaluated only confidence fusion and homopolymer balance on the same 96 held-out clusters, and reused the already-authoritative pinned BBS result from workflow run `37456317557` because the dataset, held-out seed, read ordering and BBS commit were unchanged.

The lean held-out result was:

| Reads / strand | Confidence fusion | Homopolymer balance |
| ---: | ---: | ---: |
| 5 | 73/96 (76.0%) | **75/96 (78.1%)** |
| 10 | **93/96 (96.9%)** | **93/96 (96.9%)** |

At 5 reads, homopolymer balance had **2 rescues and 0 regressions** relative to confidence fusion (two-sided exact McNemar p = 0.5). Mean edit distance improved from **0.4896 to 0.4271**.

At 10 reads, the methods were identical on exact outcomes: **0 rescues, 0 regressions**, the same **93/96** exact recovery, the same **0.0417** mean edit distance, and the same remaining **2 one-edit + 1 two-edit** failures. Homopolymer balance was also slower: about **102.8 s** versus **92.7 s** for confidence fusion across the 96 10-read clusters.

## Decision

The primary acceptance criterion required held-out 10-read exact recovery to be **strictly greater than 93/96**. The candidate did not satisfy that requirement, so it was rejected and was not merged as an accuracy improvement. No retuning was performed against the held-out failures.

The authoritative benchmark therefore remains the existing confidence-fusion result: **73/96 at 5 reads and 93/96 at 10 reads**, matching the pinned BBS exact-recovery result at 10 reads on the fixed CNR subset.

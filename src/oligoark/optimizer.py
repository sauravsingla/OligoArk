"""Deterministic search-based adaptive codec optimisation for OligoArk."""

from __future__ import annotations

import itertools
import time
from dataclasses import asdict, dataclass

from .archive import ArchiveConfig, archive_bytes, recover_bytes, recover_from_reads
from .dna import SequenceConstraintError, SequenceConstraints
from .policy import ChannelProfile
from .simulator import SimulationConfig, simulate_channel
from .tiering import WorkloadProfile


@dataclass(frozen=True)
class OptimizationWeights:
    recovery: float = 0.45
    overhead: float = 0.20
    runtime: float = 0.10
    retrieval: float = 0.10
    durability: float = 0.15

    def validate(self) -> None:
        values = vars(self)
        if any(value < 0 for value in values.values()):
            raise ValueError("optimization weights must be non-negative")
        if sum(values.values()) <= 0:
            raise ValueError("at least one optimization weight must be positive")


@dataclass(frozen=True)
class CodecSearchSpace:
    chunk_sizes: tuple[int, ...] = (64, 96)
    rs_nsyms: tuple[int, ...] = (8, 16)
    redundancy_schemes: tuple[str, ...] = ("xor", "fountain", "hybrid")
    parity_group_sizes: tuple[int, ...] = (4,)
    fountain_redundancies: tuple[float, ...] = (0.25,)
    constraints: tuple[SequenceConstraints, ...] = (
        SequenceConstraints(0.35, 0.65, 4),
        SequenceConstraints(0.40, 0.60, 4),
    )
    reconstruction_modes: tuple[str, ...] = ("direct", "graph")
    mask_search_limit: int = 96
    max_candidates: int = 48

    def validate(self) -> None:
        if not self.chunk_sizes or not self.rs_nsyms or not self.redundancy_schemes:
            raise ValueError("search space dimensions must not be empty")
        if self.max_candidates < 1:
            raise ValueError("max_candidates must be positive")
        if any(mode not in {"direct", "graph"} for mode in self.reconstruction_modes):
            raise ValueError("reconstruction modes must be direct and/or graph")
        for constraints in self.constraints:
            constraints.validate()


@dataclass(frozen=True)
class CandidateSpec:
    config: ArchiveConfig
    reconstruction_mode: str

    def key(self) -> tuple[object, ...]:
        return (*asdict(self.config).values(), self.reconstruction_mode)


@dataclass(frozen=True)
class _RawEvaluation:
    spec: CandidateSpec
    successes: int
    rate: float
    encoded: int
    overhead: float
    runtime: float
    graph_count: int
    rejected: str | None


@dataclass(frozen=True)
class CandidateEvaluation:
    config: ArchiveConfig
    reconstruction_mode: str
    trials: int
    verified_successes: int
    recovery_rate: float
    encoded_nucleotides: int
    overhead_ratio: float
    mean_runtime_seconds: float
    graph_recovery_count: int
    score: float
    rejected_reason: str | None = None

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class OptimizationResult:
    best_config: ArchiveConfig
    reconstruction_mode: str
    best_score: float
    evaluations: tuple[CandidateEvaluation, ...]
    rationale: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _candidate_specs(space: CodecSearchSpace) -> list[CandidateSpec]:
    space.validate()
    specs: list[CandidateSpec] = []
    seen: set[tuple[object, ...]] = set()
    product = itertools.product(
        space.chunk_sizes,
        space.rs_nsyms,
        space.redundancy_schemes,
        space.parity_group_sizes,
        space.fountain_redundancies,
        space.constraints,
        space.reconstruction_modes,
    )
    for chunk_size, rs_nsym, scheme, parity_group, fountain_ratio, constraints, mode in product:
        ratio = fountain_ratio if scheme in {"fountain", "hybrid"} else 0.0
        config = ArchiveConfig(
            chunk_size=chunk_size,
            rs_nsym=rs_nsym,
            parity_group_size=parity_group,
            adaptive_masks=True,
            redundancy_scheme=scheme,
            fountain_redundancy=ratio,
            min_gc_fraction=constraints.min_gc_fraction,
            max_gc_fraction=constraints.max_gc_fraction,
            max_homopolymer=constraints.max_homopolymer,
            mask_search_limit=space.mask_search_limit,
        )
        try:
            config.validate()
        except ValueError:
            continue
        spec = CandidateSpec(config, mode)
        key = spec.key()
        if key not in seen:
            seen.add(key)
            specs.append(spec)
        if len(specs) >= space.max_candidates:
            break
    if not specs:
        raise ValueError("search space produced no valid candidate configurations")
    return specs


def _simulate_once(
    payload: bytes,
    spec: CandidateSpec,
    channel: ChannelProfile,
    *,
    seed: int,
) -> tuple[bool, int, float, bool]:
    started = time.perf_counter()
    archive = archive_bytes(payload, spec.config)
    reads = simulate_channel(
        archive.strands,
        SimulationConfig(
            substitution_rate=channel.substitution_rate,
            insertion_rate=channel.insertion_rate,
            deletion_rate=channel.deletion_rate,
            dropout_rate=channel.dropout_rate,
            duplicate_rate=min(0.5, max(0.0, channel.dropout_rate * 2.0)),
            seed=seed,
        ),
    )
    graph_used = False
    recovered = False
    try:
        if spec.reconstruction_mode == "graph":
            decoded, report = recover_from_reads(archive, reads)
            graph_used = report.graph_reconstruction_used
        else:
            decoded = recover_bytes(archive, reads)
        recovered = decoded == payload
    except ValueError:
        recovered = False
    measured = archive.metadata.get("measured")
    if not isinstance(measured, dict):
        raise ValueError("archive measured metadata must be an object")
    encoded_nucleotides = int(measured["encoded_nucleotides"])
    elapsed = time.perf_counter() - started
    return recovered, encoded_nucleotides, elapsed, graph_used


def optimize_codec(
    payload: bytes,
    channel: ChannelProfile,
    workload: WorkloadProfile,
    *,
    search_space: CodecSearchSpace | None = None,
    weights: OptimizationWeights | None = None,
    seeds: tuple[int, ...] = (2026, 2027),
) -> OptimizationResult:
    """Evaluate a deterministic candidate grid and return the highest-scoring policy.

    Every success is based on exact decoded bytes after OligoArk's SHA-256 verification.
    No physical DNA performance is inferred.
    """
    if not payload:
        raise ValueError("optimizer payload must not be empty")
    channel.validate()
    workload.validate()
    if not seeds:
        raise ValueError("at least one deterministic seed is required")
    space = search_space or CodecSearchSpace()
    resolved_weights = weights or OptimizationWeights()
    resolved_weights.validate()
    specs = _candidate_specs(space)

    raw: list[_RawEvaluation] = []
    for spec in specs:
        successes = 0
        graph_count = 0
        runtimes: list[float] = []
        encoded_values: list[int] = []
        rejected_reason: str | None = None
        for seed in seeds:
            try:
                success, encoded, elapsed, graph_used = _simulate_once(
                    payload,
                    spec,
                    channel,
                    seed=seed,
                )
            except (SequenceConstraintError, ValueError) as exc:
                rejected_reason = str(exc)
                successes = 0
                runtimes = []
                encoded_values = []
                break
            successes += int(success)
            graph_count += int(graph_used)
            runtimes.append(elapsed)
            encoded_values.append(encoded)
        if rejected_reason is not None:
            raw.append(
                _RawEvaluation(
                    spec=spec,
                    successes=0,
                    rate=0.0,
                    encoded=0,
                    overhead=0.0,
                    runtime=0.0,
                    graph_count=0,
                    rejected=rejected_reason,
                )
            )
            continue
        encoded = max(encoded_values)
        ideal = max(1, len(payload) * 4)
        raw.append(
            _RawEvaluation(
                spec=spec,
                successes=successes,
                rate=successes / len(seeds),
                encoded=encoded,
                overhead=encoded / ideal,
                runtime=sum(runtimes) / len(runtimes),
                graph_count=graph_count,
                rejected=None,
            )
        )

    finite = [item for item in raw if item.rejected is None]
    if not finite:
        raise ValueError("all optimization candidates were rejected by configured constraints")
    max_overhead = max(item.overhead for item in finite) or 1.0
    max_runtime = max(item.runtime for item in finite) or 1.0
    weight_total = sum(vars(resolved_weights).values())

    evaluations: list[CandidateEvaluation] = []
    for item in raw:
        spec = item.spec
        rejected = item.rejected
        if rejected is not None:
            evaluations.append(
                CandidateEvaluation(
                    config=spec.config,
                    reconstruction_mode=spec.reconstruction_mode,
                    trials=len(seeds),
                    verified_successes=0,
                    recovery_rate=0.0,
                    encoded_nucleotides=0,
                    overhead_ratio=0.0,
                    mean_runtime_seconds=0.0,
                    graph_recovery_count=0,
                    score=-1_000_000_000.0,
                    rejected_reason=str(rejected),
                )
            )
            continue

        recovery_rate = item.rate
        overhead_norm = item.overhead / max_overhead
        runtime_norm = item.runtime / max_runtime
        retrieval_penalty = runtime_norm * workload.retrieval_urgency
        durability_reward = recovery_rate * workload.durability_priority
        score = (
            resolved_weights.recovery * recovery_rate
            - resolved_weights.overhead * overhead_norm
            - resolved_weights.runtime * runtime_norm
            - resolved_weights.retrieval * retrieval_penalty
            + resolved_weights.durability * durability_reward
        ) / weight_total
        evaluations.append(
            CandidateEvaluation(
                config=spec.config,
                reconstruction_mode=spec.reconstruction_mode,
                trials=len(seeds),
                verified_successes=item.successes,
                recovery_rate=recovery_rate,
                encoded_nucleotides=item.encoded,
                overhead_ratio=round(item.overhead, 6),
                mean_runtime_seconds=round(item.runtime, 6),
                graph_recovery_count=item.graph_count,
                score=round(score, 8),
            )
        )

    valid_evaluations = [
        evaluation for evaluation in evaluations if evaluation.rejected_reason is None
    ]
    winner = max(
        valid_evaluations,
        key=lambda evaluation: (
            evaluation.score,
            evaluation.recovery_rate,
            -evaluation.overhead_ratio,
            -evaluation.mean_runtime_seconds,
        ),
    )
    rationale = (
        f"evaluated {len(evaluations)} candidate configurations across {len(seeds)} seed(s)",
        "success means exact SHA-256-verified byte recovery",
        (
            f"winner recovery_rate={winner.recovery_rate:.3f}, "
            f"overhead_ratio={winner.overhead_ratio:.3f}, "
            f"runtime={winner.mean_runtime_seconds:.6f}s"
        ),
        (
            f"selected redundancy={winner.config.redundancy_scheme}, "
            f"rs_nsym={winner.config.rs_nsym}, "
            f"reconstruction={winner.reconstruction_mode}"
        ),
        "objective weights are explicit software-research preferences, not physical DNA claims",
    )
    return OptimizationResult(
        best_config=winner.config,
        reconstruction_mode=winner.reconstruction_mode,
        best_score=winner.score,
        evaluations=tuple(evaluations),
        rationale=rationale,
    )

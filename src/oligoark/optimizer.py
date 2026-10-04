"""Deterministic, explainable search-based codec optimisation for OligoArk."""

from __future__ import annotations

import hashlib
import itertools
import time
from collections.abc import Mapping
from dataclasses import asdict, dataclass

from .archive import ArchiveConfig, archive_bytes, recover_bytes, recover_from_reads
from .dna import SequenceConstraintError, SequenceConstraints
from .policy import ChannelProfile
from .simulator import SimulationConfig, simulate_channel
from .tiering import WorkloadProfile


@dataclass(frozen=True)
class OptimizationWeights:
    """Explicit objective weights. Physical terms are used only when supplied."""

    recovery: float = 0.40
    overhead: float = 0.15
    redundancy: float = 0.10
    runtime: float = 0.08
    retrieval: float = 0.07
    durability: float = 0.10
    lifecycle_storage_cost: float = 0.04
    lifecycle_retrieval_cost: float = 0.02
    lifecycle_energy: float = 0.02
    lifecycle_latency: float = 0.02

    def validate(self) -> None:
        if any(value < 0 for value in vars(self).values()):
            raise ValueError("optimization weights must be non-negative")
        if sum(vars(self).values()) <= 0:
            raise ValueError("at least one optimization weight must be positive")

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> OptimizationWeights:
        defaults = cls()
        allowed = set(vars(defaults))
        unknown = sorted(set(values) - allowed)
        if unknown:
            raise ValueError(f"Unknown optimization weight(s): {unknown}")

        def number(name: str) -> float:
            value = values.get(name, getattr(defaults, name))
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"Optimization weight {name} must be numeric")
            return float(value)

        result = cls(
            recovery=number("recovery"),
            overhead=number("overhead"),
            redundancy=number("redundancy"),
            runtime=number("runtime"),
            retrieval=number("retrieval"),
            durability=number("durability"),
            lifecycle_storage_cost=number("lifecycle_storage_cost"),
            lifecycle_retrieval_cost=number("lifecycle_retrieval_cost"),
            lifecycle_energy=number("lifecycle_energy"),
            lifecycle_latency=number("lifecycle_latency"),
        )
        result.validate()
        return result


@dataclass(frozen=True)
class LifecycleObjectiveInputs:
    """Optional caller-derived lifecycle values for the selected storage tier.

    Values are physical/user-supplied quantities. Candidate overhead scales storage,
    retrieval and energy quantities; retrieval latency is treated as an operation-level
    quantity and therefore does not scale with nucleotide overhead.
    """

    storage_cost: float
    retrieval_cost: float
    energy_kwh: float
    retrieval_latency_hours: float

    def validate(self) -> None:
        if any(value < 0 for value in vars(self).values()):
            raise ValueError("lifecycle objective inputs must be non-negative")


@dataclass(frozen=True)
class CodecSearchSpace:
    chunk_sizes: tuple[int, ...] = (48, 64, 96)
    rs_nsyms: tuple[int, ...] = (8, 16, 24)
    redundancy_schemes: tuple[str, ...] = ("xor", "fountain", "hybrid")
    parity_group_sizes: tuple[int, ...] = (3, 5)
    fountain_redundancies: tuple[float, ...] = (0.25, 0.40)
    constraints: tuple[SequenceConstraints, ...] = (
        SequenceConstraints(0.35, 0.65, 4),
        SequenceConstraints(0.40, 0.60, 4),
    )
    reconstruction_modes: tuple[str, ...] = ("direct", "graph")
    mask_search_limit: int = 96
    max_candidates: int = 48
    search_method: str = "balanced"
    search_seed: int = 5050

    def validate(self) -> None:
        if not self.chunk_sizes or not self.rs_nsyms or not self.redundancy_schemes:
            raise ValueError("search space dimensions must not be empty")
        if not self.parity_group_sizes or not self.fountain_redundancies:
            raise ValueError("search space redundancy dimensions must not be empty")
        if not self.constraints or not self.reconstruction_modes:
            raise ValueError("search space constraint/reconstruction dimensions must not be empty")
        if self.max_candidates < 1:
            raise ValueError("max_candidates must be positive")
        if self.search_method not in {"balanced", "full_grid"}:
            raise ValueError("search_method must be 'balanced' or 'full_grid'")
        if any(mode not in {"direct", "graph"} for mode in self.reconstruction_modes):
            raise ValueError("reconstruction modes must be direct and/or graph")
        for constraints in self.constraints:
            constraints.validate()


@dataclass(frozen=True)
class CandidateSpec:
    config: ArchiveConfig
    reconstruction_mode: str

    def key(self) -> tuple[object, ...]:
        return (
            self.config.chunk_size,
            self.config.rs_nsym,
            self.config.redundancy_scheme,
            self.config.parity_group_size,
            self.config.fountain_redundancy,
            self.config.min_gc_fraction,
            self.config.max_gc_fraction,
            self.config.max_homopolymer,
            self.reconstruction_mode,
        )


@dataclass(frozen=True)
class ObjectiveBreakdown:
    recovery_reward: float
    overhead_penalty: float
    redundancy_penalty: float
    runtime_penalty: float
    retrieval_penalty: float
    durability_reward: float
    lifecycle_storage_cost_penalty: float = 0.0
    lifecycle_retrieval_cost_penalty: float = 0.0
    lifecycle_energy_penalty: float = 0.0
    lifecycle_latency_penalty: float = 0.0

    @property
    def total(self) -> float:
        return round(
            self.recovery_reward
            - self.overhead_penalty
            - self.redundancy_penalty
            - self.runtime_penalty
            - self.retrieval_penalty
            + self.durability_reward
            - self.lifecycle_storage_cost_penalty
            - self.lifecycle_retrieval_cost_penalty
            - self.lifecycle_energy_penalty
            - self.lifecycle_latency_penalty,
            10,
        )


@dataclass(frozen=True)
class CandidateEvaluation:
    config: ArchiveConfig
    reconstruction_mode: str
    trials: int
    verified_successes: int
    recovery_rate: float
    encoded_nucleotides: int
    overhead_ratio: float
    redundancy_ratio: float
    mean_runtime_seconds: float
    graph_recovery_count: int
    score: float
    objective: ObjectiveBreakdown
    rejected_reason: str | None = None

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class OptimizationResult:
    best_config: ArchiveConfig
    reconstruction_mode: str
    best_score: float
    evaluations: tuple[CandidateEvaluation, ...]
    total_possible_candidates: int
    evaluated_candidates: int
    rejected_candidates: int
    search_method: str
    search_seed: int
    calibration_seeds: tuple[int, ...]
    rationale: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class _RawEvaluation:
    spec: CandidateSpec
    successes: int
    encoded: int
    overhead: float
    redundancy: float
    runtime: float
    graph_count: int
    rejected: str | None


def _canonical_values(space: CodecSearchSpace) -> tuple[
    tuple[int, ...],
    tuple[int, ...],
    tuple[str, ...],
    tuple[int, ...],
    tuple[float, ...],
    tuple[SequenceConstraints, ...],
    tuple[str, ...],
]:
    constraints = tuple(
        sorted(
            set(space.constraints),
            key=lambda c: (c.min_gc_fraction, c.max_gc_fraction, c.max_homopolymer),
        )
    )
    return (
        tuple(sorted(set(space.chunk_sizes))),
        tuple(sorted(set(space.rs_nsyms))),
        tuple(sorted(set(space.redundancy_schemes))),
        tuple(sorted(set(space.parity_group_sizes))),
        tuple(sorted(set(space.fountain_redundancies))),
        constraints,
        tuple(sorted(set(space.reconstruction_modes))),
    )


def enumerate_candidate_specs(space: CodecSearchSpace) -> tuple[CandidateSpec, ...]:
    """Return canonical valid candidates before search-budget sampling."""
    space.validate()
    (
        chunks,
        rs_values,
        schemes,
        groups,
        fountain_values,
        constraints_values,
        modes,
    ) = _canonical_values(space)
    specs: dict[tuple[object, ...], CandidateSpec] = {}
    for chunk, rs, scheme, group, fountain, constraints, mode in itertools.product(
        chunks,
        rs_values,
        schemes,
        groups,
        fountain_values,
        constraints_values,
        modes,
    ):
        ratio = fountain if scheme in {"fountain", "hybrid"} else 0.0
        config = ArchiveConfig(
            chunk_size=chunk,
            rs_nsym=rs,
            parity_group_size=group,
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
        specs[spec.key()] = spec
    return tuple(specs[key] for key in sorted(specs, key=lambda value: repr(value)))


def select_candidate_specs(space: CodecSearchSpace) -> tuple[CandidateSpec, ...]:
    """Select a deterministic order-independent subset from the canonical full grid."""
    all_specs = list(enumerate_candidate_specs(space))
    if not all_specs:
        raise ValueError("search space produced no valid candidate configurations")
    if space.search_method == "full_grid" or len(all_specs) <= space.max_candidates:
        return tuple(all_specs)

    by_scheme_mode: dict[tuple[str, str], list[CandidateSpec]] = {}
    for spec in all_specs:
        by_scheme_mode.setdefault(
            (spec.config.redundancy_scheme, spec.reconstruction_mode), []
        ).append(spec)

    selected: list[CandidateSpec] = []
    groups = sorted(by_scheme_mode)
    round_index = 0
    while len(selected) < space.max_candidates:
        progressed = False
        for group in groups:
            candidates = by_scheme_mode[group]
            if round_index >= len(candidates):
                continue
            seed_material = f"{space.search_seed}|{group}|{round_index}".encode()
            digest = int(hashlib.sha256(seed_material).hexdigest()[:16], 16)
            index = digest % len(candidates)
            candidate = candidates[index]
            if candidate in selected:
                remaining = [item for item in candidates if item not in selected]
                if not remaining:
                    continue
                candidate = remaining[digest % len(remaining)]
            selected.append(candidate)
            progressed = True
            if len(selected) >= space.max_candidates:
                break
        if not progressed:
            break
        round_index += 1
    return tuple(sorted(selected, key=lambda spec: repr(spec.key())))


def _simulate_once(
    payload: bytes,
    spec: CandidateSpec,
    channel: ChannelProfile,
    *,
    seed: int,
    duplicate_rate: float | None,
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
            duplicate_rate=(
                duplicate_rate
                if duplicate_rate is not None
                else min(0.65, max(0.0, channel.dropout_rate * 2.0))
            ),
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
    elapsed = time.perf_counter() - started
    return recovered, int(measured["encoded_nucleotides"]), elapsed, graph_used


def _redundancy_ratio(config: ArchiveConfig) -> float:
    """Return an explicit software redundancy proxy used by the objective.

    This combines frame-level Reed-Solomon parity with archive-level XOR/fountain
    redundancy. It is an optimization feature, not a physical synthesis overhead metric.
    """
    protected_payload_bytes = max(1, config.chunk_size + 18)
    rs_component = config.rs_nsym / protected_payload_bytes
    xor_component = 0.0
    if config.redundancy_scheme in {"xor", "hybrid"}:
        xor_component = 1.0 / config.parity_group_size
    fountain_component = (
        config.fountain_redundancy
        if config.redundancy_scheme in {"fountain", "hybrid"}
        else 0.0
    )
    return rs_component + xor_component + fountain_component


def optimize_codec(
    payload: bytes,
    channel: ChannelProfile,
    workload: WorkloadProfile,
    *,
    search_space: CodecSearchSpace | None = None,
    weights: OptimizationWeights | None = None,
    seeds: tuple[int, ...] = (2026, 2027),
    lifecycle: LifecycleObjectiveInputs | None = None,
    duplicate_rate: float | None = None,
) -> OptimizationResult:
    """Search candidates on calibration seeds only and return an explainable winner."""
    if not payload:
        raise ValueError("optimizer payload must not be empty")
    channel.validate()
    workload.validate()
    if not seeds:
        raise ValueError("at least one calibration seed is required")
    if len(set(seeds)) != len(seeds):
        raise ValueError("calibration seeds must be unique")
    if duplicate_rate is not None and not 0 <= duplicate_rate <= 1:
        raise ValueError("duplicate_rate must be between 0 and 1")
    space = search_space or CodecSearchSpace()
    resolved_weights = weights or OptimizationWeights()
    resolved_weights.validate()
    if lifecycle is not None:
        lifecycle.validate()

    all_specs = enumerate_candidate_specs(space)
    specs = select_candidate_specs(space)
    raw: list[_RawEvaluation] = []
    for spec in specs:
        successes = 0
        graph_count = 0
        runtimes: list[float] = []
        encoded_values: list[int] = []
        rejected: str | None = None
        for seed in seeds:
            try:
                success, encoded, elapsed, graph_used = _simulate_once(
                    payload,
                    spec,
                    channel,
                    seed=seed,
                    duplicate_rate=duplicate_rate,
                )
            except (SequenceConstraintError, ValueError) as exc:
                rejected = str(exc)
                break
            successes += int(success)
            graph_count += int(graph_used)
            runtimes.append(elapsed)
            encoded_values.append(encoded)
        if rejected is not None:
            raw.append(_RawEvaluation(spec, 0, 0, 0.0, 0.0, 0.0, 0, rejected))
            continue
        encoded = max(encoded_values)
        ideal = max(1, len(payload) * 4)
        raw.append(
            _RawEvaluation(
                spec,
                successes,
                encoded,
                encoded / ideal,
                _redundancy_ratio(spec.config),
                sum(runtimes) / len(runtimes),
                graph_count,
                None,
            )
        )

    valid = [item for item in raw if item.rejected is None]
    if not valid:
        raise ValueError("all optimization candidates were rejected")

    max_overhead = max(item.overhead for item in valid) or 1.0
    min_overhead = min(item.overhead for item in valid) or 1.0
    max_redundancy = max(item.redundancy for item in valid) or 1.0
    max_runtime = max(item.runtime for item in valid) or 1.0
    max_footprint_scale = max_overhead / min_overhead
    lifecycle_storage_max = (
        max(1e-12, lifecycle.storage_cost * max_footprint_scale) if lifecycle else 1.0
    )
    lifecycle_retrieval_max = (
        max(1e-12, lifecycle.retrieval_cost * max_footprint_scale) if lifecycle else 1.0
    )
    lifecycle_energy_max = (
        max(1e-12, lifecycle.energy_kwh * max_footprint_scale) if lifecycle else 1.0
    )
    lifecycle_latency_max = (
        max(1e-12, lifecycle.retrieval_latency_hours) if lifecycle else 1.0
    )

    active_weight_names = [
        "recovery",
        "overhead",
        "redundancy",
        "runtime",
        "retrieval",
        "durability",
    ]
    if lifecycle is not None:
        active_weight_names.extend(
            [
                "lifecycle_storage_cost",
                "lifecycle_retrieval_cost",
                "lifecycle_energy",
                "lifecycle_latency",
            ]
        )
    weight_total = sum(getattr(resolved_weights, name) for name in active_weight_names) or 1.0

    evaluations: list[CandidateEvaluation] = []
    for item in raw:
        if item.rejected is not None:
            empty = ObjectiveBreakdown(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            evaluations.append(
                CandidateEvaluation(
                    item.spec.config,
                    item.spec.reconstruction_mode,
                    len(seeds),
                    0,
                    0.0,
                    0,
                    0.0,
                    0.0,
                    0.0,
                    0,
                    -1_000_000_000.0,
                    empty,
                    item.rejected,
                )
            )
            continue

        recovery_rate = item.successes / len(seeds)
        overhead_norm = item.overhead / max_overhead
        redundancy_norm = item.redundancy / max_redundancy
        runtime_norm = item.runtime / max_runtime
        footprint_scale = item.overhead / min_overhead
        lifecycle_storage = lifecycle.storage_cost * footprint_scale if lifecycle else 0.0
        lifecycle_retrieval = lifecycle.retrieval_cost * footprint_scale if lifecycle else 0.0
        lifecycle_energy = lifecycle.energy_kwh * footprint_scale if lifecycle else 0.0
        lifecycle_latency = lifecycle.retrieval_latency_hours if lifecycle else 0.0

        breakdown = ObjectiveBreakdown(
            recovery_reward=resolved_weights.recovery * recovery_rate / weight_total,
            overhead_penalty=resolved_weights.overhead * overhead_norm / weight_total,
            redundancy_penalty=resolved_weights.redundancy * redundancy_norm / weight_total,
            runtime_penalty=resolved_weights.runtime * runtime_norm / weight_total,
            retrieval_penalty=(
                resolved_weights.retrieval
                * runtime_norm
                * workload.retrieval_urgency
                / weight_total
            ),
            durability_reward=(
                resolved_weights.durability
                * recovery_rate
                * workload.durability_priority
                / weight_total
            ),
            lifecycle_storage_cost_penalty=(
                resolved_weights.lifecycle_storage_cost
                * lifecycle_storage
                / lifecycle_storage_max
                / weight_total
                if lifecycle
                else 0.0
            ),
            lifecycle_retrieval_cost_penalty=(
                resolved_weights.lifecycle_retrieval_cost
                * lifecycle_retrieval
                / lifecycle_retrieval_max
                / weight_total
                if lifecycle
                else 0.0
            ),
            lifecycle_energy_penalty=(
                resolved_weights.lifecycle_energy
                * lifecycle_energy
                / lifecycle_energy_max
                / weight_total
                if lifecycle
                else 0.0
            ),
            lifecycle_latency_penalty=(
                resolved_weights.lifecycle_latency
                * lifecycle_latency
                / lifecycle_latency_max
                / weight_total
                if lifecycle
                else 0.0
            ),
        )
        evaluations.append(
            CandidateEvaluation(
                config=item.spec.config,
                reconstruction_mode=item.spec.reconstruction_mode,
                trials=len(seeds),
                verified_successes=item.successes,
                recovery_rate=round(recovery_rate, 6),
                encoded_nucleotides=item.encoded,
                overhead_ratio=round(item.overhead, 6),
                redundancy_ratio=round(item.redundancy, 6),
                mean_runtime_seconds=round(item.runtime, 6),
                graph_recovery_count=item.graph_count,
                score=breakdown.total,
                objective=breakdown,
            )
        )

    usable = [item for item in evaluations if item.rejected_reason is None]
    winner = max(
        usable,
        key=lambda item: (
            item.score,
            item.recovery_rate,
            -item.overhead_ratio,
            -item.redundancy_ratio,
            -item.mean_runtime_seconds,
        ),
    )
    rationale = (
        (
            f"evaluated {len(evaluations)} of {len(all_specs)} candidates "
            f"using search_method={space.search_method}"
        ),
        f"search_seed={space.search_seed}",
        f"calibration_seeds={list(seeds)}",
        "success requires normal recovery and original SHA-256 verification",
        (
            f"winner recovery={winner.recovery_rate:.3f}, overhead={winner.overhead_ratio:.3f}, "
            f"redundancy={winner.redundancy_ratio:.3f}, score={winner.score:.6f}"
        ),
        "physical lifecycle terms are omitted unless caller-supplied inputs are provided",
    )
    return OptimizationResult(
        best_config=winner.config,
        reconstruction_mode=winner.reconstruction_mode,
        best_score=winner.score,
        evaluations=tuple(evaluations),
        total_possible_candidates=len(all_specs),
        evaluated_candidates=len(evaluations),
        rejected_candidates=sum(item.rejected_reason is not None for item in evaluations),
        search_method=space.search_method,
        search_seed=space.search_seed,
        calibration_seeds=tuple(seeds),
        rationale=rationale,
    )

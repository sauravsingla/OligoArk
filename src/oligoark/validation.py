"""Controlled reconstruction validation and ablation helpers."""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass

from .archive import DNAArchive, recover_bytes
from .reconstruct import (
    GraphConsensusReconstructor,
    ReconstructionResult,
    TraceConsensusReconstructor,
)


@dataclass(frozen=True)
class RescueModeResult:
    mode: str
    recovered: bool
    verified_sha256: bool
    node_count: int
    candidate_pairs: int
    edge_count: int
    component_count: int
    cluster_sizes: tuple[int, ...]
    consensus_lengths: tuple[int, ...]
    consensus_count: int
    reconstruction_runtime_seconds: float
    total_runtime_seconds: float
    failure: str | None = None

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class GraphRescueComparison:
    direct_recovered: bool
    direct_failure: str | None
    medoid: RescueModeResult
    alignment: RescueModeResult
    trace: RescueModeResult
    alignment_rescued_direct_failure: bool
    trace_rescued_direct_failure: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _attempt_reconstruction(
    archive: DNAArchive,
    reads: list[str],
    *,
    mode: str,
    threshold: float,
) -> RescueModeResult:
    started = time.perf_counter()
    reconstructor = GraphConsensusReconstructor(
        threshold=threshold,
        consensus_mode=mode,
    )
    reconstruction: ReconstructionResult = reconstructor.reconstruct(reads)
    recovered = False
    failure: str | None = None
    try:
        recover_bytes(archive, reads + reconstruction.consensus_reads)
        recovered = True
    except ValueError as exc:
        failure = str(exc)
    return RescueModeResult(
        mode=mode,
        recovered=recovered,
        verified_sha256=recovered,
        node_count=reconstruction.node_count,
        candidate_pairs=reconstruction.candidate_pairs,
        edge_count=reconstruction.edge_count,
        component_count=reconstruction.component_count,
        cluster_sizes=tuple(reconstruction.cluster_sizes),
        consensus_lengths=reconstruction.consensus_lengths,
        consensus_count=len(reconstruction.consensus_reads),
        reconstruction_runtime_seconds=reconstruction.runtime_seconds,
        total_runtime_seconds=round(time.perf_counter() - started, 8),
        failure=failure,
    )


def _attempt_trace_reconstruction(
    archive: DNAArchive,
    reads: list[str],
    *,
    threshold: float,
) -> RescueModeResult:
    started = time.perf_counter()
    thresholds = tuple(
        value
        for value in (
            min(0.99, threshold + 0.04),
            threshold,
            max(0.0, threshold - 0.04),
        )
    )
    reconstructor = TraceConsensusReconstructor(thresholds=thresholds)
    reconstruction = reconstructor.reconstruct(reads)
    recovered = False
    failure: str | None = None
    try:
        recover_bytes(archive, reads + reconstruction.consensus_reads)
        recovered = True
    except ValueError as exc:
        failure = str(exc)
    return RescueModeResult(
        mode="trace",
        recovered=recovered,
        verified_sha256=recovered,
        node_count=reconstruction.node_count,
        candidate_pairs=reconstruction.candidate_pairs,
        edge_count=reconstruction.edge_count,
        component_count=reconstruction.component_count,
        cluster_sizes=tuple(reconstruction.cluster_sizes),
        consensus_lengths=reconstruction.consensus_lengths,
        consensus_count=len(reconstruction.consensus_reads),
        reconstruction_runtime_seconds=reconstruction.runtime_seconds,
        total_runtime_seconds=round(time.perf_counter() - started, 8),
        failure=failure,
    )


def compare_reconstruction_modes(
    archive: DNAArchive,
    reads: list[str],
    *,
    threshold: float = 0.90,
) -> GraphRescueComparison:
    """Compare direct, graph-medoid and graph-alignment recovery on the same reads."""
    archive.validate()
    clean_reads = [read.strip().upper() for read in reads if read.strip()]
    if not clean_reads:
        raise ValueError("reads must not be empty")

    direct_recovered = False
    direct_failure: str | None = None
    try:
        recover_bytes(archive, clean_reads)
        direct_recovered = True
    except ValueError as exc:
        direct_failure = str(exc)

    medoid = _attempt_reconstruction(
        archive,
        clean_reads,
        mode="medoid",
        threshold=threshold,
    )
    alignment = _attempt_reconstruction(
        archive,
        clean_reads,
        mode="alignment",
        threshold=threshold,
    )
    trace = _attempt_trace_reconstruction(
        archive,
        clean_reads,
        threshold=threshold,
    )
    return GraphRescueComparison(
        direct_recovered=direct_recovered,
        direct_failure=direct_failure,
        medoid=medoid,
        alignment=alignment,
        trace=trace,
        alignment_rescued_direct_failure=(
            not direct_recovered and alignment.recovered and alignment.verified_sha256
        ),
        trace_rescued_direct_failure=(
            not direct_recovered and trace.recovered and trace.verified_sha256
        ),
    )

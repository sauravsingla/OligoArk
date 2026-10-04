"""FastAPI application for OligoArk."""

from __future__ import annotations

import base64
import json

try:
    from fastapi import FastAPI, HTTPException
    from pydantic import BaseModel, Field
except ImportError as exc:  # pragma: no cover
    raise RuntimeError("Install the API extra with: pip install 'oligoark[api]'") from exc

from . import __version__
from .archive import ArchiveConfig, DNAArchive, archive_bytes, recover_bytes, recover_from_reads
from .config import RuntimeConfig
from .intelligence import (
    evaluate_optimized_archive_plan,
    optimize_archive_plan,
    plan_archive,
)
from .logging_utils import configure_logging
from .optimizer import CodecSearchSpace, OptimizationWeights
from .policy import ChannelProfile, PolicyObjective, recommend_codec_policy
from .simulator import SimulationConfig, simulate_channel
from .tiering import (
    EconomicAssumptions,
    LifecycleAssumptions,
    WorkloadProfile,
    recommend_storage_tier,
)
from .validation import compare_reconstruction_modes

RUNTIME_CONFIG = RuntimeConfig.from_environment()
configure_logging(RUNTIME_CONFIG.log_level)
app = FastAPI(title="OligoArk", version=__version__)


class EncodeRequest(BaseModel):
    data_b64: str
    chunk_size: int = Field(default=96, ge=8)
    rs_nsym: int = Field(default=8, ge=0, le=64)
    parity_group_size: int = Field(default=8, ge=2)
    redundancy_scheme: str = "xor"
    fountain_redundancy: float = Field(default=0.25, ge=0, le=5)
    min_gc_fraction: float = Field(default=0.35, ge=0, le=1)
    max_gc_fraction: float = Field(default=0.65, ge=0, le=1)
    max_homopolymer: int = Field(default=4, ge=1)
    mask_search_limit: int = Field(default=64, ge=1, le=256)


class RecoverRequest(BaseModel):
    archive: dict[str, object]


class RecoverReadsRequest(BaseModel):
    archive: dict[str, object]
    reads: list[str]
    similarity_threshold: float = Field(
        default=RUNTIME_CONFIG.reconstruction_threshold,
        ge=0,
        le=1,
    )


class ReconstructionDiagnosticsRequest(RecoverReadsRequest):
    pass


class SimulateRequest(BaseModel):
    archive: dict[str, object]
    substitution_rate: float = Field(default=0.0, ge=0, le=1)
    insertion_rate: float = Field(default=0.0, ge=0, le=1)
    deletion_rate: float = Field(default=0.0, ge=0, le=1)
    dropout_rate: float = Field(default=0.0, ge=0, le=1)
    duplicate_rate: float = Field(default=0.0, ge=0, le=1)
    seed: int = 7


class EconomicRequest(BaseModel):
    storage_cost_index: dict[str, float]
    retrieval_cost_index: dict[str, float]


class LifecycleTierRequest(BaseModel):
    storage_cost_per_gb_year: float = Field(ge=0)
    retrieval_cost_per_gb: float = Field(ge=0)
    idle_energy_kwh_per_tb_year: float = Field(ge=0)
    retrieval_energy_kwh_per_gb: float = Field(ge=0)
    retrieval_latency_hours: float = Field(ge=0)


class LifecycleRequest(BaseModel):
    tiers: dict[str, LifecycleTierRequest]


class TierRequest(BaseModel):
    retention_years: float = Field(gt=0)
    accesses_per_year: float = Field(default=0.0, ge=0)
    mutability: float = Field(default=0.0, ge=0, le=1)
    retrieval_urgency: float = Field(default=0.0, ge=0, le=1)
    durability_priority: float = Field(default=1.0, ge=0, le=1)
    energy_priority: float = Field(default=0.5, ge=0, le=1)
    redundancy_priority: float = Field(default=0.5, ge=0, le=1)
    cost_priority: float = Field(default=0.5, ge=0, le=1)
    data_size_gb: float = Field(default=1.0, gt=0)
    expected_access_probability: float | None = Field(default=None, ge=0, le=1)
    economics: EconomicRequest | None = None
    lifecycle: LifecycleRequest | None = None


class PolicyRequest(BaseModel):
    substitution_rate: float = Field(default=0.0, ge=0, le=1)
    insertion_rate: float = Field(default=0.0, ge=0, le=1)
    deletion_rate: float = Field(default=0.0, ge=0, le=1)
    dropout_rate: float = Field(default=0.0, ge=0, le=1)
    durability_priority: float = Field(default=0.7, ge=0, le=1)
    storage_overhead_priority: float = Field(default=0.2, ge=0, le=1)
    retrieval_speed_priority: float = Field(default=0.1, ge=0, le=1)


class PlanRequest(TierRequest):
    substitution_rate: float = Field(default=0.0, ge=0, le=1)
    insertion_rate: float = Field(default=0.0, ge=0, le=1)
    deletion_rate: float = Field(default=0.0, ge=0, le=1)
    dropout_rate: float = Field(default=0.0, ge=0, le=1)


class OptimizationWeightsRequest(BaseModel):
    recovery: float = Field(default=0.40, ge=0)
    overhead: float = Field(default=0.15, ge=0)
    redundancy: float = Field(default=0.10, ge=0)
    runtime: float = Field(default=0.08, ge=0)
    retrieval: float = Field(default=0.07, ge=0)
    durability: float = Field(default=0.10, ge=0)
    lifecycle_storage_cost: float = Field(default=0.04, ge=0)
    lifecycle_retrieval_cost: float = Field(default=0.02, ge=0)
    lifecycle_energy: float = Field(default=0.02, ge=0)
    lifecycle_latency: float = Field(default=0.02, ge=0)


class OptimizePlanRequest(PlanRequest):
    data_b64: str
    seeds: list[int] | None = None
    calibration_seeds: list[int] = Field(default_factory=lambda: [2026, 2027])
    evaluation_seeds: list[int] | None = None
    duplicate_rate: float = Field(default=0.0, ge=0, le=1)
    max_candidates: int = Field(default=24, ge=1, le=256)
    search_method: str = "balanced"
    search_seed: int = 5050
    weights: OptimizationWeightsRequest | None = None


def _archive_from_mapping(value: dict[str, object]) -> DNAArchive:
    encoded = json.dumps(value).encode("utf-8")
    if len(encoded) > RUNTIME_CONFIG.max_api_payload_bytes:
        raise ValueError("Archive request exceeds configured API payload limit")
    return DNAArchive.from_json(encoded.decode("utf-8"))


def _decode_payload(data_b64: str) -> bytes:
    raw = base64.b64decode(data_b64, validate=True)
    if len(raw) > RUNTIME_CONFIG.max_api_payload_bytes:
        raise ValueError("Decoded input exceeds configured API payload limit")
    return raw


def _economics(value: EconomicRequest | None) -> EconomicAssumptions | None:
    if value is None:
        return None
    return EconomicAssumptions.from_mapping(value.model_dump())


def _lifecycle(value: LifecycleRequest | None) -> LifecycleAssumptions | None:
    if value is None:
        return None
    raw = {tier: assumptions.model_dump() for tier, assumptions in value.tiers.items()}
    return LifecycleAssumptions.from_mapping(raw)


def _workload(req: TierRequest) -> WorkloadProfile:
    return WorkloadProfile(
        retention_years=req.retention_years,
        accesses_per_year=req.accesses_per_year,
        mutability=req.mutability,
        retrieval_urgency=req.retrieval_urgency,
        durability_priority=req.durability_priority,
        energy_priority=req.energy_priority,
        redundancy_priority=req.redundancy_priority,
        cost_priority=req.cost_priority,
        data_size_gb=req.data_size_gb,
        expected_access_probability=req.expected_access_probability,
    )


def _channel(req: PlanRequest) -> ChannelProfile:
    return ChannelProfile(
        req.substitution_rate,
        req.insertion_rate,
        req.deletion_rate,
        req.dropout_rate,
    )


def _weights(req: OptimizationWeightsRequest | None) -> OptimizationWeights | None:
    if req is None:
        return None
    return OptimizationWeights.from_mapping(req.model_dump())


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "oligoark", "version": __version__}


@app.post("/encode")
def encode(req: EncodeRequest) -> dict[str, object]:
    try:
        raw = _decode_payload(req.data_b64)
        archive = archive_bytes(
            raw,
            ArchiveConfig(
                chunk_size=req.chunk_size,
                rs_nsym=req.rs_nsym,
                parity_group_size=req.parity_group_size,
                adaptive_masks=True,
                redundancy_scheme=req.redundancy_scheme,
                fountain_redundancy=req.fountain_redundancy,
                min_gc_fraction=req.min_gc_fraction,
                max_gc_fraction=req.max_gc_fraction,
                max_homopolymer=req.max_homopolymer,
                mask_search_limit=req.mask_search_limit,
            ),
        )
        return {"metadata": archive.metadata, "strands": archive.strands}
    except (ValueError, TypeError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/recover")
def recover(req: RecoverRequest) -> dict[str, str]:
    try:
        raw = recover_bytes(_archive_from_mapping(req.archive))
        return {"data_b64": base64.b64encode(raw).decode("ascii")}
    except (ValueError, TypeError, KeyError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/recover-reads")
def recover_reads(req: RecoverReadsRequest) -> dict[str, object]:
    try:
        archive = _archive_from_mapping(req.archive)
        raw, report = recover_from_reads(
            archive,
            req.reads,
            similarity_threshold=req.similarity_threshold,
        )
        return {
            "data_b64": base64.b64encode(raw).decode("ascii"),
            "report": report.to_dict(),
        }
    except (ValueError, TypeError, KeyError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/reconstruction-diagnostics")
def reconstruction_diagnostics(req: ReconstructionDiagnosticsRequest) -> dict[str, object]:
    try:
        return compare_reconstruction_modes(
            _archive_from_mapping(req.archive),
            req.reads,
            threshold=req.similarity_threshold,
        ).to_dict()
    except (ValueError, TypeError, KeyError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/simulate")
def simulate(req: SimulateRequest) -> dict[str, object]:
    try:
        archive = _archive_from_mapping(req.archive)
        config = SimulationConfig(
            req.substitution_rate,
            req.insertion_rate,
            req.deletion_rate,
            req.dropout_rate,
            req.duplicate_rate,
            req.seed,
        )
        return {"reads": simulate_channel(archive.strands, config)}
    except (ValueError, TypeError, KeyError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/recommend")
def recommend(req: TierRequest) -> dict[str, object]:
    return recommend_storage_tier(
        _workload(req),
        _economics(req.economics),
        _lifecycle(req.lifecycle),
    ).to_dict()


@app.post("/policy")
def policy(req: PolicyRequest) -> dict[str, object]:
    return recommend_codec_policy(
        ChannelProfile(
            req.substitution_rate,
            req.insertion_rate,
            req.deletion_rate,
            req.dropout_rate,
        ),
        PolicyObjective(
            req.durability_priority,
            req.storage_overhead_priority,
            req.retrieval_speed_priority,
        ),
    ).to_dict()


@app.post("/plan")
def plan(req: PlanRequest) -> dict[str, object]:
    return plan_archive(
        _workload(req),
        _channel(req),
        economics=_economics(req.economics),
        lifecycle=_lifecycle(req.lifecycle),
    ).to_dict()


@app.post("/optimize-plan")
def optimize_plan(req: OptimizePlanRequest) -> dict[str, object]:
    try:
        payload = _decode_payload(req.data_b64)
        calibration = tuple(req.seeds or req.calibration_seeds)
        search = CodecSearchSpace(
            max_candidates=req.max_candidates,
            search_method=req.search_method,
            search_seed=req.search_seed,
        )
        economics = _economics(req.economics)
        lifecycle = _lifecycle(req.lifecycle)
        weights = _weights(req.weights)
        if req.evaluation_seeds:
            return evaluate_optimized_archive_plan(
                payload,
                _workload(req),
                _channel(req),
                calibration_seeds=calibration,
                evaluation_seeds=tuple(req.evaluation_seeds),
                duplicate_rate=req.duplicate_rate,
                economics=economics,
                lifecycle=lifecycle,
                search_space=search,
                weights=weights,
            ).to_dict()
        return optimize_archive_plan(
            payload,
            _workload(req),
            _channel(req),
            seeds=calibration,
            economics=economics,
            lifecycle=lifecycle,
            search_space=search,
            weights=weights,
            duplicate_rate=req.duplicate_rate,
        ).to_dict()
    except (ValueError, TypeError, KeyError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

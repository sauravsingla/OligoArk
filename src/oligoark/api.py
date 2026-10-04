"""FastAPI application for OligoArk."""

from __future__ import annotations

import base64
import json

try:
    from fastapi import FastAPI, HTTPException
    from pydantic import BaseModel, Field
except ImportError as exc:  # pragma: no cover
    raise RuntimeError("Install the API extra with: pip install 'oligoark[api]'") from exc

from .archive import ArchiveConfig, DNAArchive, archive_bytes, recover_bytes, recover_from_reads
from .config import RuntimeConfig
from .logging_utils import configure_logging
from .policy import ChannelProfile, PolicyObjective, recommend_codec_policy
from .simulator import SimulationConfig, simulate_channel
from .tiering import EconomicAssumptions, WorkloadProfile, recommend_storage_tier

RUNTIME_CONFIG = RuntimeConfig.from_environment()
configure_logging(RUNTIME_CONFIG.log_level)
app = FastAPI(title="OligoArk", version="0.2.0")


class EncodeRequest(BaseModel):
    data_b64: str
    chunk_size: int = Field(default=96, ge=8)
    rs_nsym: int = Field(default=8, ge=0, le=64)
    parity_group_size: int = Field(default=8, ge=2)


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


class TierRequest(BaseModel):
    retention_years: float = Field(gt=0)
    accesses_per_year: float = Field(ge=0)
    mutability: float = Field(ge=0, le=1)
    retrieval_urgency: float = Field(ge=0, le=1)
    durability_priority: float = Field(ge=0, le=1)
    energy_priority: float = Field(ge=0, le=1)
    redundancy_priority: float = Field(default=0.5, ge=0, le=1)
    cost_priority: float = Field(default=0.5, ge=0, le=1)
    economics: EconomicRequest | None = None


class PolicyRequest(BaseModel):
    substitution_rate: float = Field(default=0.0, ge=0, le=1)
    insertion_rate: float = Field(default=0.0, ge=0, le=1)
    deletion_rate: float = Field(default=0.0, ge=0, le=1)
    dropout_rate: float = Field(default=0.0, ge=0, le=1)
    durability_priority: float = Field(default=0.7, ge=0, le=1)
    storage_overhead_priority: float = Field(default=0.2, ge=0, le=1)
    retrieval_speed_priority: float = Field(default=0.1, ge=0, le=1)


def _archive_from_mapping(value: dict[str, object]) -> DNAArchive:
    encoded = json.dumps(value).encode("utf-8")
    if len(encoded) > RUNTIME_CONFIG.max_api_payload_bytes:
        raise ValueError("Archive request exceeds configured API payload limit")
    return DNAArchive.from_json(encoded.decode("utf-8"))


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "oligoark"}


@app.post("/encode")
def encode(req: EncodeRequest) -> dict[str, object]:
    try:
        raw = base64.b64decode(req.data_b64, validate=True)
        if len(raw) > RUNTIME_CONFIG.max_api_payload_bytes:
            raise ValueError("Decoded input exceeds configured API payload limit")
        archive = archive_bytes(
            raw,
            ArchiveConfig(req.chunk_size, req.rs_nsym, req.parity_group_size, True),
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
    economics = None
    if req.economics is not None:
        economics = EconomicAssumptions.from_mapping(req.economics.model_dump())
    profile = WorkloadProfile(
        retention_years=req.retention_years,
        accesses_per_year=req.accesses_per_year,
        mutability=req.mutability,
        retrieval_urgency=req.retrieval_urgency,
        durability_priority=req.durability_priority,
        energy_priority=req.energy_priority,
        redundancy_priority=req.redundancy_priority,
        cost_priority=req.cost_priority,
    )
    return recommend_storage_tier(profile, economics).to_dict()


@app.post("/policy")
def policy(req: PolicyRequest) -> dict[str, object]:
    channel = ChannelProfile(
        req.substitution_rate,
        req.insertion_rate,
        req.deletion_rate,
        req.dropout_rate,
    )
    objective = PolicyObjective(
        req.durability_priority,
        req.storage_overhead_priority,
        req.retrieval_speed_priority,
    )
    return recommend_codec_policy(channel, objective).to_dict()

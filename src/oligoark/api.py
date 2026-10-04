"""FastAPI application for OligoArk."""

from __future__ import annotations
import base64
import json

try:
    from fastapi import FastAPI, HTTPException
    from pydantic import BaseModel, Field
except ImportError as exc:
    raise RuntimeError("Install the API extra with: pip install 'oligoark[api]'") from exc

from .archive import ArchiveConfig, DNAArchive, archive_bytes, recover_bytes
from .tiering import WorkloadProfile, recommend_storage_tier

app = FastAPI(title="OligoArk", version="0.1.0")


class EncodeRequest(BaseModel):
    data_b64: str
    chunk_size: int = Field(default=96, ge=8)
    rs_nsym: int = Field(default=8, ge=0, le=64)
    parity_group_size: int = Field(default=8, ge=2)


class RecoverRequest(BaseModel):
    archive: dict[str, object]


class TierRequest(BaseModel):
    retention_years: float = Field(gt=0)
    accesses_per_year: float = Field(ge=0)
    mutability: float = Field(ge=0, le=1)
    retrieval_urgency: float = Field(ge=0, le=1)
    durability_priority: float = Field(ge=0, le=1)
    energy_priority: float = Field(ge=0, le=1)


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "oligoark"}


@app.post("/encode")
def encode(req: EncodeRequest) -> dict[str, object]:
    try:
        raw = base64.b64decode(req.data_b64, validate=True)
        archive = archive_bytes(raw,
            ArchiveConfig(req.chunk_size, req.rs_nsym, req.parity_group_size, True))
        return {"metadata": archive.metadata, "strands": archive.strands}
    except (ValueError, TypeError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/recover")
def recover(req: RecoverRequest) -> dict[str, str]:
    try:
        archive = DNAArchive.from_json(json.dumps(req.archive))
        raw = recover_bytes(archive)
        return {"data_b64": base64.b64encode(raw).decode("ascii")}
    except (ValueError, TypeError, KeyError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/recommend")
def recommend(req: TierRequest) -> dict[str, object]:
    return recommend_storage_tier(WorkloadProfile(**req.model_dump())).to_dict()

import base64

from fastapi.testclient import TestClient

from oligoark.api import app

client = TestClient(app)


def test_api_encode_recover_and_health() -> None:
    health = client.get("/health").json()
    assert health["status"] == "ok"
    assert health["version"] == "0.6.0"

    payload = b"api-roundtrip"
    encoded = client.post(
        "/encode",
        json={
            "data_b64": base64.b64encode(payload).decode("ascii"),
            "redundancy_scheme": "hybrid",
            "fountain_redundancy": 0.5,
        },
    )
    assert encoded.status_code == 200
    archive = encoded.json()
    recovered = client.post("/recover", json={"archive": archive})
    assert recovered.status_code == 200
    assert base64.b64decode(recovered.json()["data_b64"]) == payload


def test_api_simulate_policy_tiering_and_plan() -> None:
    payload = base64.b64encode(b"simulation").decode("ascii")
    archive = client.post("/encode", json={"data_b64": payload}).json()
    simulated = client.post("/simulate", json={"archive": archive, "seed": 11})
    assert simulated.status_code == 200
    assert simulated.json()["reads"]

    policy = client.post("/policy", json={"substitution_rate": 0.01})
    assert policy.status_code == 200
    assert policy.json()["rs_nsym"] >= 8

    tier_payload = {
        "retention_years": 100,
        "accesses_per_year": 0.1,
        "mutability": 0,
        "retrieval_urgency": 0.1,
        "durability_priority": 1,
        "energy_priority": 0.8,
        "data_size_gb": 0.1,
        "economics": {
            "storage_cost_index": {
                "ssd": 0.8,
                "object_archive": 0.2,
                "tape": 0.3,
                "dna_future": 0.5,
            },
            "retrieval_cost_index": {
                "ssd": 0.1,
                "object_archive": 0.4,
                "tape": 0.6,
                "dna_future": 0.8,
            },
        },
    }
    tier = client.post("/recommend", json=tier_payload)
    assert tier.status_code == 200
    assert tier.json()["recommended_tier"] in tier.json()["scores"]

    plan = client.post(
        "/plan",
        json={
            **tier_payload,
            "substitution_rate": 0.01,
            "dropout_rate": 0.05,
        },
    )
    assert plan.status_code == 200
    result = plan.json()
    assert result["policy_source"] == "deterministic-heuristic"
    assert result["tier"]["recommended_tier"] in result["tier"]["scores"]
    assert result["codec_policy"]["rs_nsym"] >= 16


def test_api_measured_optimizer_endpoint() -> None:
    payload = base64.b64encode(b"optimizer api payload").decode("ascii")
    response = client.post(
        "/optimize-plan",
        json={
            "data_b64": payload,
            "retention_years": 100,
            "substitution_rate": 0.002,
            "max_candidates": 4,
            "seeds": [2026],
        },
    )
    assert response.status_code == 200
    result = response.json()
    assert result["optimization"]["evaluations"]
    assert result["selected_redundancy_scheme"] in {"xor", "fountain", "hybrid"}


def test_api_held_out_optimizer_and_reconstruction_diagnostics() -> None:
    payload = b"api held-out validation"
    response = client.post(
        "/optimize-plan",
        json={
            "data_b64": base64.b64encode(payload).decode("ascii"),
            "retention_years": 50,
            "substitution_rate": 0.001,
            "max_candidates": 4,
            "calibration_seeds": [9001],
            "evaluation_seeds": [2026],
            "search_method": "balanced",
            "search_seed": 77,
            "weights": {"recovery": 1.0, "overhead": 0.1},
        },
    )
    assert response.status_code == 200
    held_out = response.json()
    assert held_out["calibration_seeds"] == [9001]
    assert held_out["evaluation_seeds"] == [2026]
    assert len(held_out["trials"]) == 1
    assert held_out["plan"]["optimization"]["search_method"] == "balanced"

    encoded = client.post(
        "/encode",
        json={
            "data_b64": base64.b64encode(b"diagnostic rescue").decode("ascii"),
            "rs_nsym": 0,
            "redundancy_scheme": "none",
            "mask_search_limit": 128,
        },
    ).json()
    strand = encoded["strands"][0]
    reads = []
    for position in (64, 76, 88, 100, 112):
        if position >= len(strand):
            continue
        current = strand[position]
        replacement = next(base for base in "ACGT" if base != current)
        reads.append(strand[:position] + replacement + strand[position + 1 :])
    diagnostic = client.post(
        "/reconstruction-diagnostics",
        json={
            "archive": encoded,
            "reads": reads,
            "similarity_threshold": 0.96,
        },
    )
    assert diagnostic.status_code == 200
    result = diagnostic.json()
    assert result["direct_recovered"] is False
    assert result["alignment"]["verified_sha256"] is True
    assert result["alignment"]["node_count"] == len(reads)

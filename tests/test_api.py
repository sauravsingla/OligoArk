import base64

from fastapi.testclient import TestClient

from oligoark.api import app

client = TestClient(app)


def test_api_encode_recover_and_health() -> None:
    assert client.get("/health").json()["status"] == "ok"
    payload = b"api-roundtrip"
    encoded = client.post(
        "/encode",
        json={"data_b64": base64.b64encode(payload).decode("ascii")},
    )
    assert encoded.status_code == 200
    archive = encoded.json()
    recovered = client.post("/recover", json={"archive": archive})
    assert recovered.status_code == 200
    assert base64.b64decode(recovered.json()["data_b64"]) == payload


def test_api_simulate_policy_and_tiering() -> None:
    payload = base64.b64encode(b"simulation").decode("ascii")
    archive = client.post("/encode", json={"data_b64": payload}).json()
    simulated = client.post("/simulate", json={"archive": archive, "seed": 11})
    assert simulated.status_code == 200
    assert simulated.json()["reads"]

    policy = client.post("/policy", json={"substitution_rate": 0.01})
    assert policy.status_code == 200
    assert policy.json()["rs_nsym"] >= 8

    tier = client.post(
        "/recommend",
        json={
            "retention_years": 100,
            "accesses_per_year": 0.1,
            "mutability": 0,
            "retrieval_urgency": 0.1,
            "durability_priority": 1,
            "energy_priority": 0.8,
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
        },
    )
    assert tier.status_code == 200
    assert tier.json()["recommended_tier"] in tier.json()["scores"]

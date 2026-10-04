from pathlib import Path

import pytest

from oligoark.config import RuntimeConfig
from oligoark.learning import EmpiricalPolicyModel, PolicyObservation
from oligoark.policy import ChannelProfile, CodecPolicy


def test_runtime_config_roundtrip(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "config.json"
    path.write_text(
        '{"log_level":"DEBUG","reconstruction_threshold":0.85,"max_api_payload_bytes":4096}',
        encoding="utf-8",
    )
    config = RuntimeConfig.from_json_file(path)
    assert config.log_level == "DEBUG"
    assert config.reconstruction_threshold == pytest.approx(0.85)

    monkeypatch.setenv("OLIGOARK_LOG_LEVEL", "ERROR")
    monkeypatch.setenv("OLIGOARK_RECONSTRUCTION_THRESHOLD", "0.8")
    env_config = RuntimeConfig.from_environment()
    assert env_config.log_level == "ERROR"
    assert env_config.reconstruction_threshold == pytest.approx(0.8)


def test_runtime_config_rejects_unknown_fields() -> None:
    with pytest.raises(ValueError, match="Unknown runtime configuration"):
        RuntimeConfig.from_mapping({"mystery": 1})


def test_empirical_policy_model_learns_nearby_success() -> None:
    conservative = CodecPolicy(48, 24, 3, True, ("conservative",))
    lean = CodecPolicy(96, 8, 8, True, ("lean",))
    observations = [
        PolicyObservation(ChannelProfile(0.03, 0, 0, 0.10), conservative, True, 2000, 0.2),
        PolicyObservation(ChannelProfile(0.03, 0, 0, 0.10), lean, False, 1200, 0.1),
        PolicyObservation(ChannelProfile(0.0, 0, 0, 0.0), lean, True, 1200, 0.1),
    ]
    recommendation = EmpiricalPolicyModel().fit(observations).recommend(
        ChannelProfile(0.028, 0, 0, 0.09)
    )
    assert recommendation.policy == conservative
    assert 0 <= recommendation.confidence <= 1

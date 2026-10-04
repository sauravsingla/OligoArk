import pytest

from oligoark.policy import ChannelProfile, recommend_codec_policy
from oligoark.tiering import WorkloadProfile, recommend_storage_tier


def test_adaptive_policy_strengthens_on_noisy_channel() -> None:
    clean = recommend_codec_policy(ChannelProfile())
    noisy = recommend_codec_policy(ChannelProfile(substitution_rate=0.03, dropout_rate=0.12))
    assert noisy.rs_nsym > clean.rs_nsym
    assert noisy.chunk_size < clean.chunk_size
    assert noisy.parity_group_size < clean.parity_group_size


def test_tiering_is_explainable_and_deterministic() -> None:
    p = WorkloadProfile(200, 0.01, 0.0, 0.0, 1.0, 1.0)
    first = recommend_storage_tier(p)
    assert first == recommend_storage_tier(p)
    assert first.recommended_tier in first.scores
    assert first.rationale


def test_invalid_profile_rejected() -> None:
    with pytest.raises(ValueError):
        recommend_storage_tier(WorkloadProfile(0, 0, 0, 0, 1, 1))

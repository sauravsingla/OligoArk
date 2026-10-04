import pytest

from oligoark.policy import ChannelProfile, PolicyObjective, recommend_codec_policy
from oligoark.tiering import EconomicAssumptions, WorkloadProfile, recommend_storage_tier


def test_adaptive_policy_strengthens_on_noisy_channel() -> None:
    clean = recommend_codec_policy(ChannelProfile())
    noisy = recommend_codec_policy(ChannelProfile(substitution_rate=0.03, dropout_rate=0.12))
    assert noisy.rs_nsym > clean.rs_nsym
    assert noisy.chunk_size < clean.chunk_size
    assert noisy.parity_group_size < clean.parity_group_size


def test_durability_objective_strengthens_policy() -> None:
    channel = ChannelProfile(substitution_rate=0.002)
    baseline = recommend_codec_policy(channel)
    durable = recommend_codec_policy(channel, PolicyObjective(durability_priority=0.95))
    assert durable.rs_nsym >= baseline.rs_nsym
    assert durable.parity_group_size <= baseline.parity_group_size


def test_tiering_is_explainable_and_deterministic() -> None:
    profile = WorkloadProfile(200, 0.01, 0.0, 0.0, 1.0, 1.0)
    first = recommend_storage_tier(profile)
    assert first == recommend_storage_tier(profile)
    assert first.recommended_tier in first.scores
    assert first.rationale


def test_economic_assumptions_can_change_ranking() -> None:
    profile = WorkloadProfile(20, 1, 0.2, 0.2, 0.5, 0.5, cost_priority=1.0)
    neutral = recommend_storage_tier(profile)
    assumptions = EconomicAssumptions(
        storage_cost_index={
            "ssd": 1.0,
            "object_archive": 0.0,
            "tape": 1.0,
            "dna_future": 1.0,
        },
        retrieval_cost_index={
            "ssd": 1.0,
            "object_archive": 0.0,
            "tape": 1.0,
            "dna_future": 1.0,
        },
    )
    custom = recommend_storage_tier(profile, assumptions)
    assert custom.scores["object_archive"] > neutral.scores["object_archive"]
    reconstructed = EconomicAssumptions.from_mapping(
        {
            "storage_cost_index": assumptions.storage_cost_index,
            "retrieval_cost_index": assumptions.retrieval_cost_index,
        }
    )
    assert reconstructed == assumptions


def test_invalid_profile_rejected() -> None:
    with pytest.raises(ValueError):
        recommend_storage_tier(WorkloadProfile(0, 0, 0, 0, 1, 1))


def test_invalid_channel_rejected() -> None:
    with pytest.raises(ValueError):
        recommend_codec_policy(ChannelProfile(substitution_rate=1.1))

"""Non-finite numeric inputs must fail before codec, optimizer or policy work."""

import math

import pytest

from oligoark.archive import ArchiveConfig
from oligoark.config import RuntimeConfig
from oligoark.dna import SequenceConstraints
from oligoark.optimizer import CodecSearchSpace, LifecycleObjectiveInputs, OptimizationWeights
from oligoark.policy import ChannelProfile, PolicyObjective
from oligoark.simulator import SimulationConfig
from oligoark.streaming import StreamingFaultProfile
from oligoark.tiering import (
    EconomicAssumptions,
    TierLifecycleAssumption,
    WorkloadProfile,
)


@pytest.mark.parametrize("invalid", [math.nan, math.inf, -math.inf])
def test_nonfinite_archive_and_encoding_inputs_are_rejected(invalid: float) -> None:
    with pytest.raises(ValueError):
        SequenceConstraints(min_gc_fraction=invalid).validate()
    with pytest.raises(ValueError):
        SequenceConstraints(max_gc_fraction=invalid).validate()
    with pytest.raises(ValueError):
        ArchiveConfig(fountain_redundancy=invalid).validate()
    with pytest.raises(ValueError):
        ArchiveConfig.from_mapping({"fountain_redundancy": invalid})
    with pytest.raises(ValueError):
        CodecSearchSpace(fountain_redundancies=(invalid,)).validate()


@pytest.mark.parametrize("invalid", [math.nan, math.inf, -math.inf])
def test_nonfinite_runtime_optimizer_and_noise_inputs_are_rejected(invalid: float) -> None:
    with pytest.raises(ValueError):
        RuntimeConfig(reconstruction_threshold=invalid).validate()
    with pytest.raises(ValueError):
        OptimizationWeights(recovery=invalid).validate()
    with pytest.raises(ValueError):
        LifecycleObjectiveInputs(invalid, 1.0, 1.0, 1.0).validate()
    with pytest.raises(ValueError):
        ChannelProfile(substitution_rate=invalid).validate()
    with pytest.raises(ValueError):
        PolicyObjective(durability_priority=invalid).validate()
    with pytest.raises(ValueError):
        SimulationConfig(substitution_rate=invalid).validate()
    with pytest.raises(ValueError):
        StreamingFaultProfile(substitution_rate=invalid).validate()


@pytest.mark.parametrize("invalid", [math.nan, math.inf, -math.inf])
def test_nonfinite_lifecycle_and_workload_inputs_are_rejected(invalid: float) -> None:
    valid_workload = dict(
        retention_years=5.0,
        accesses_per_year=0.5,
        mutability=0.1,
        retrieval_urgency=0.2,
        durability_priority=0.8,
        energy_priority=0.5,
    )
    for field in (
        "retention_years",
        "accesses_per_year",
        "mutability",
        "data_size_gb",
        "expected_access_probability",
    ):
        with pytest.raises(ValueError):
            WorkloadProfile(**{**valid_workload, field: invalid}).validate()

    economics = EconomicAssumptions()
    economics.storage_cost_index["ssd"] = invalid
    with pytest.raises(ValueError):
        economics.validate()
    with pytest.raises(ValueError):
        TierLifecycleAssumption(invalid, 1.0, 1.0, 1.0, 1.0).validate()

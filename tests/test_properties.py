import pytest

from oligoark.archive import ArchiveConfig, archive_bytes, recover_bytes
from oligoark.dna import bytes_to_dna, dna_to_bytes
from oligoark.ecc import rs_decode, rs_encode

hypothesis = pytest.importorskip("hypothesis")
st = pytest.importorskip("hypothesis.strategies")
given = hypothesis.given
settings = hypothesis.settings


@given(st.binary(max_size=512))
@settings(max_examples=60, deadline=None)
def test_binary_dna_property_roundtrip(payload: bytes) -> None:
    assert dna_to_bytes(bytes_to_dna(payload)) == payload


@given(st.binary(min_size=1, max_size=128))
@settings(max_examples=40, deadline=None)
def test_archive_property_roundtrip(payload: bytes) -> None:
    config = ArchiveConfig(chunk_size=48, rs_nsym=12, parity_group_size=4)
    assert recover_bytes(archive_bytes(payload, config)) == payload


@given(st.binary(min_size=1, max_size=100))
@settings(max_examples=40, deadline=None)
def test_rs_property_roundtrip(payload: bytes) -> None:
    assert rs_decode(rs_encode(payload, 12), 12) == payload

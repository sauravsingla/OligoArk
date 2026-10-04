import pytest

from oligoark.archive import ArchiveConfig, archive_bytes, recover_bytes
from oligoark.dna import bytes_to_dna, dna_to_bytes
from oligoark.ecc import rs_decode, rs_encode
from oligoark.optimizer import CodecSearchSpace, select_candidate_specs
from oligoark.reconstruct import build_similarity_graph

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


@given(
    st.permutations((32, 48, 64)),
    st.permutations((8, 12, 16)),
    st.permutations(("xor", "fountain", "hybrid")),
)
@settings(max_examples=20, deadline=None)
def test_balanced_search_subset_is_property_order_invariant(
    chunks: tuple[int, ...],
    rs_values: tuple[int, ...],
    schemes: tuple[str, ...],
) -> None:
    baseline = CodecSearchSpace(
        chunk_sizes=(32, 48, 64),
        rs_nsyms=(8, 12, 16),
        redundancy_schemes=("xor", "fountain", "hybrid"),
        parity_group_sizes=(3,),
        fountain_redundancies=(0.25,),
        reconstruction_modes=("direct", "graph"),
        max_candidates=8,
        search_method="balanced",
        search_seed=123,
    )
    permuted = CodecSearchSpace(
        chunk_sizes=chunks,
        rs_nsyms=rs_values,
        redundancy_schemes=schemes,
        parity_group_sizes=(3,),
        fountain_redundancies=(0.25,),
        reconstruction_modes=("graph", "direct"),
        constraints=tuple(reversed(baseline.constraints)),
        max_candidates=8,
        search_method="balanced",
        search_seed=123,
    )
    assert [item.key() for item in select_candidate_specs(baseline)] == [
        item.key() for item in select_candidate_specs(permuted)
    ]


@given(st.text(alphabet="ACGT", min_size=24, max_size=80))
@settings(max_examples=30, deadline=None)
def test_qgram_prefilter_preserves_exact_qualifying_edges(sequence: str) -> None:
    position = len(sequence) // 2
    replacement = {"A": "C", "C": "G", "G": "T", "T": "A"}[sequence[position]]
    mutated = sequence[:position] + replacement + sequence[position + 1 :]
    reads = [sequence, mutated, sequence[::-1]]
    indexed = build_similarity_graph(
        reads,
        threshold=0.86,
        qgram_width=5,
        use_qgram_prefilter=True,
    )
    exhaustive = build_similarity_graph(
        reads,
        threshold=0.86,
        qgram_width=5,
        use_qgram_prefilter=False,
    )
    assert {(edge.left, edge.right) for edge in indexed.edges} == {
        (edge.left, edge.right) for edge in exhaustive.edges
    }

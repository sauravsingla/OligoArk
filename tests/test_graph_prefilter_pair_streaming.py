"""Regression for graph q-gram prefilter after streaming pair enumeration."""

from oligoark.reconstruct import build_similarity_graph


def test_qgram_pair_streaming_preserves_qualified_graph_edges():
    reads = [
        "AAAAAAACCC",
        "AAAAAAACCT",
        "GGGGGGGTTA",
        "GGGGGGGTTC",
        "ACGTACGTAC",
    ]
    fast = build_similarity_graph(
        reads, threshold=0.85, qgram_width=5, use_qgram_prefilter=True
    )
    exhaustive = build_similarity_graph(
        reads, threshold=0.85, qgram_width=5, use_qgram_prefilter=False
    )
    assert fast.nodes == exhaustive.nodes
    assert fast.edges == exhaustive.edges
    assert fast.connected_components() == exhaustive.connected_components()
    assert fast.candidate_pairs <= exhaustive.candidate_pairs
    assert exhaustive.candidate_pairs == 10

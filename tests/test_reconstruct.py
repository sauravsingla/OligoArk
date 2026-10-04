from oligoark.reconstruct import (
    GraphConsensusReconstructor,
    ReconstructionResult,
    edit_distance,
    graph_cluster_consensus,
    normalized_similarity,
)


class ConstantEdgeScorer:
    def score(self, left: str, right: str) -> float:
        return 1.0


def test_edit_distance() -> None:
    assert edit_distance("ACGT", "ACCT") == 1
    assert normalized_similarity("ACGT", "ACCT") == 0.75


def test_graph_cluster_consensus_groups_similar_reads() -> None:
    reads = ["ACGTACGT", "ACGTACGA", "TTTTGGGG", "TTTTGGGA"]
    result = graph_cluster_consensus(reads, threshold=0.80)
    assert sorted(result.cluster_sizes) == [2, 2]
    assert len(result.consensus_reads) == 2


def test_custom_edge_scorer_extension_point() -> None:
    reads = ["AAAA", "TTTT"]
    reconstructor = GraphConsensusReconstructor(
        threshold=0.9,
        scorer=ConstantEdgeScorer(),
        qgram_width=2,
        use_qgram_prefilter=False,
    )
    result = reconstructor.reconstruct(reads)
    assert result.cluster_sizes == [2]
    assert isinstance(result, ReconstructionResult)


def test_invalid_graph_configuration_rejected() -> None:
    reconstructor = GraphConsensusReconstructor(threshold=1.1)
    try:
        reconstructor.reconstruct(["ACGT"])
    except ValueError as exc:
        assert "threshold" in str(exc)
    else:
        raise AssertionError("invalid threshold should fail")

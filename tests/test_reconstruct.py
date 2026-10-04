from oligoark.reconstruct import edit_distance, graph_cluster_consensus, normalized_similarity


def test_edit_distance() -> None:
    assert edit_distance("ACGT", "ACCT") == 1
    assert normalized_similarity("ACGT", "ACCT") == 0.75


def test_graph_cluster_consensus_groups_similar_reads() -> None:
    reads = ["ACGTACGT", "ACGTACGA", "TTTTGGGG", "TTTTGGGA"]
    result = graph_cluster_consensus(reads, threshold=0.80)
    assert sorted(result.cluster_sizes) == [2, 2]
    assert len(result.consensus_reads) == 2

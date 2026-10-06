from oligoark.reconstruct import (
    GraphConsensusReconstructor,
    ReconstructionResult,
    alignment_consensus,
    build_similarity_graph,
    confidence_fusion_trace_consensus,
    edit_distance,
    global_align,
    graph_cluster_consensus,
    homopolymer_balance_trace_consensus,
    multistart_trace_consensus,
    normalized_similarity,
    stability_fusion_trace_consensus,
    targeted_trace_consensus,
)


class ConstantEdgeScorer:
    def score(self, left: str, right: str) -> float:
        return 1.0


def test_edit_distance() -> None:
    assert edit_distance("ACGT", "ACCT") == 1
    assert normalized_similarity("ACGT", "ACCT") == 0.75


def test_explicit_graph_has_edges_and_components() -> None:
    reads = ["ACGTACGT", "ACGTACGA", "TTTTGGGG", "TTTTGGGA"]
    graph = build_similarity_graph(reads, threshold=0.80)
    assert len(graph.nodes) == 4
    assert len(graph.edges) == 2
    assert sorted(len(component) for component in graph.connected_components()) == [2, 2]


def test_graph_cluster_consensus_groups_similar_reads() -> None:
    reads = ["ACGTACGT", "ACGTACGA", "TTTTGGGG", "TTTTGGGA"]
    result = graph_cluster_consensus(reads, threshold=0.80)
    assert sorted(result.cluster_sizes) == [2, 2]
    assert result.edge_count == 2
    assert result.component_count == 2


def test_alignment_consensus_repairs_simple_indels() -> None:
    original = "ACGTACGT"
    reads = [
        original,
        "ACGTTACGT",
        "ACGACGT",
        original,
        original,
    ]
    assert alignment_consensus(reads) == original
    aligned_reference, aligned_query = global_align(original, "ACGTTACGT")
    assert len(aligned_reference) == len(aligned_query)
    assert "-" in aligned_reference


def test_custom_edge_scorer_extension_point() -> None:
    reads = ["AAAA", "TTTT"]
    reconstructor = GraphConsensusReconstructor(
        threshold=0.9,
        scorer=ConstantEdgeScorer(),
        qgram_width=2,
        use_qgram_prefilter=False,
        consensus_mode="medoid",
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


def test_qgram_index_reduces_unrelated_graph_candidates() -> None:
    reads = [
        "A" * 80 + "C" * 20,
        "A" * 79 + "G" + "C" * 20,
        "CGTAC" * 20,
        "TGCAT" * 20,
    ]
    indexed = build_similarity_graph(
        reads,
        threshold=0.90,
        qgram_width=5,
        use_qgram_prefilter=True,
    )
    exhaustive = build_similarity_graph(
        reads,
        threshold=0.90,
        qgram_width=5,
        use_qgram_prefilter=False,
    )
    assert indexed.candidate_pairs < exhaustive.candidate_pairs
    assert {(edge.left, edge.right) for edge in indexed.edges} == {
        (edge.left, edge.right) for edge in exhaustive.edges
    }

def test_multistart_trace_consensus_repairs_anchor_bias() -> None:
    original = (
        "TTAGTTGTGCCGCAGCGAAGTAGTGCTTGAAATATGCGACCCCTAAGTAGGAGCGTATGCGCCC"
        "AGTAACCAATGCCTGTTGAGATGCCAGACGCGTAACCAAAACATAG"
    )
    reads = [
        "TTAGTTGTGCCGTAGTGAAGTAGTGCTTAGAAATAATGCGACCCCTAAGTAGGAGCGTATGCGCCCAGTAACCAATGCTTGTTGGGATGCCAGCCGCGTAACCAAAACATAG",
        "TTAGTTGTGCCGCAGCTAGAGTAGTGCTTGAAATTGCGACCCCTAAGTGGGAGCGTATGCGCCCAGAACAATGCCTGTTGAGATGCCAAGACGCGTAACCAAAACATAG",
        "TTAGTTGTGCCGCAGCGAAGTCAGTGCTTGAAGATGCACCCCTAAGTAGGAGCTATGCGCCCAGTAACCATTGCCTGTTGAGATGCAGACGCGTAGCCAAAACATAG",
        "TTAGTTGTGCCGCAGCAAGTAGTGCTTGAAATATGCGAGCCCCTAAGGGAAGCGTATGCGCCCAGTAACCAATGCCTGTTGAGATGCCATACGCGTAACCAAAACATAG",
        "TTAGTTGTGCCGGCAGCGAAGTAGTGCTTGAAATATAGCGACCCCCTAAGTAGGAGCTGTATGCGCCCAGTAACCAAATGCCTGTTGATATCCCAGACGCGTAACCAATAACATAG",
    ]
    result = multistart_trace_consensus(
        reads,
        target_length=len(original),
        anchors=2,
        rounds=1,
        bidirectional=True,
        length_penalty=1.0,
    )
    assert result == original


def test_multistart_trace_consensus_validates_configuration() -> None:
    reads = ["ACGT", "ACGA"]
    for kwargs in (
        {"anchors": 0},
        {"rounds": 0},
        {"target_length": 0},
        {"length_penalty": -0.1},
    ):
        try:
            multistart_trace_consensus(reads, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"configuration should fail: {kwargs}")

def test_targeted_trace_consensus_restores_known_length_homopolymer_deletion() -> None:
    original = "AAAACCCC"
    reads = ["AAACCCC", "AAACCCC", original]
    assert targeted_trace_consensus(
        reads,
        target_length=len(original),
        anchors=2,
        rounds=1,
        bidirectional=True,
        homopolymer_weight=1.0,
    ) == original


def test_targeted_trace_consensus_high_coverage_homopolymer_substitution() -> None:
    baseline = "ACAGGGTA"
    original = "ACGGGGTA"
    reads = [baseline] * 6 + [original] * 4
    assert targeted_trace_consensus(
        reads,
        target_length=len(original),
        anchors=3,
        rounds=1,
        bidirectional=True,
        homopolymer_weight=0.75,
        substitution_homopolymer_weight=2.0,
        substitution_homopolymer_min_reads=10,
    ) == original


def test_targeted_trace_consensus_keeps_low_coverage_prior_disabled() -> None:
    baseline = "ACAGGGTA"
    original = "ACGGGGTA"
    reads = [baseline] * 3 + [original] * 2
    assert targeted_trace_consensus(
        reads,
        target_length=len(original),
        anchors=3,
        rounds=1,
        bidirectional=True,
        homopolymer_weight=0.75,
        substitution_homopolymer_weight=2.0,
        substitution_homopolymer_min_reads=10,
    ) == baseline


def test_targeted_trace_consensus_validates_configuration() -> None:
    reads = ["ACGT", "ACGA"]
    for kwargs in (
        {"target_length": 0},
        {"target_length": 4, "anchors": 0},
        {"target_length": 4, "rounds": 0},
        {"target_length": 4, "length_penalty": -0.1},
        {"target_length": 4, "homopolymer_weight": -0.1},
        {"target_length": 4, "min_homopolymer_run": 0},
        {"target_length": 4, "substitution_min_gain": -0.1},
        {"target_length": 4, "substitution_homopolymer_weight": -0.1},
        {"target_length": 4, "substitution_homopolymer_min_reads": 0},
    ):
        try:
            targeted_trace_consensus(reads, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"configuration should fail: {kwargs}")

def test_confidence_fusion_trace_consensus_preserves_supported_consensus() -> None:
    original = "AAAACCCC"
    reads = ["AAACCCC", "AAACCCC", original, original, original]
    assert confidence_fusion_trace_consensus(
        reads,
        target_length=len(original),
        anchors=2,
        rounds=1,
        top_positions=2,
        max_candidates=12,
        trim_farthest=0,
        qgram_width=3,
        qgram_weight=0.25,
        minimum_score_gain=0.01,
    ) == original


def test_confidence_fusion_trace_consensus_validates_configuration() -> None:
    reads = ["ACGT", "ACGA"]
    invalid = (
        {"target_length": 0},
        {"target_length": 4, "top_positions": 0},
        {"target_length": 4, "max_candidates": 0},
        {"target_length": 4, "trim_farthest": 2},
        {"target_length": 4, "qgram_width": 0},
        {"target_length": 4, "qgram_weight": -0.1},
        {"target_length": 4, "minimum_score_gain": -0.1},
    )
    for kwargs in invalid:
        try:
            confidence_fusion_trace_consensus(reads, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"configuration should fail: {kwargs}")

def test_stability_fusion_is_deterministic_under_read_order() -> None:
    reads = [
        "AAAACCCC",
        "AAACCCC",
        "AAAACCCC",
        "AAAACCCC",
        "AAAACCCA",
    ]
    forward = stability_fusion_trace_consensus(
        reads,
        target_length=8,
        top_positions=2,
        max_candidates=10,
        reliability_power=1.0,
        support_weight=0.25,
        confidence_weight=0.5,
        qgram_width=3,
        qgram_weight=0.25,
        minimum_score_gain=0.05,
        skip_margin=3,
    )
    reverse_order = stability_fusion_trace_consensus(
        list(reversed(reads)),
        target_length=8,
        top_positions=2,
        max_candidates=10,
        reliability_power=1.0,
        support_weight=0.25,
        confidence_weight=0.5,
        qgram_width=3,
        qgram_weight=0.25,
        minimum_score_gain=0.05,
        skip_margin=3,
    )
    assert forward == reverse_order


def test_stability_fusion_preserves_strong_supported_consensus() -> None:
    original = "ACGTACGT"
    reads = [original] * 8 + ["ACGTTCGT", "ACGTACG"]
    assert stability_fusion_trace_consensus(
        reads,
        target_length=len(original),
        top_positions=2,
        max_candidates=10,
        reliability_power=1.0,
        support_weight=0.25,
        confidence_weight=0.5,
        qgram_width=4,
        qgram_weight=0.25,
        minimum_score_gain=0.05,
        skip_margin=3,
    ) == original


def test_stability_fusion_validates_configuration() -> None:
    reads = ["ACGT", "ACGA"]
    invalid = (
        {"target_length": 0},
        {"target_length": 4, "top_positions": 0},
        {"target_length": 4, "max_candidates": 1},
        {"target_length": 4, "reliability_power": -0.1},
        {"target_length": 4, "support_weight": -0.1},
        {"target_length": 4, "confidence_weight": -0.1},
        {"target_length": 4, "qgram_width": 0},
        {"target_length": 4, "qgram_weight": -0.1},
        {"target_length": 4, "minimum_score_gain": -0.1},
        {"target_length": 4, "skip_margin": -1},
    )
    for kwargs in invalid:
        try:
            stability_fusion_trace_consensus(reads, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"configuration should fail: {kwargs}")

def test_homopolymer_balance_preserves_strong_supported_consensus() -> None:
    original = "AAACCCCC"
    reads = [original] * 8 + ["AAAACCCC", "AAACCCCC"]
    assert homopolymer_balance_trace_consensus(
        reads,
        target_length=len(original),
        top_runs=2,
        beam_width=6,
        depth=1,
        min_gap_support=2,
        homopolymer_weight=0.15,
        qgram_width=4,
        qgram_weight=0.25,
        minimum_score_gain=0.05,
    ) == original


def test_homopolymer_balance_is_deterministic_under_read_order() -> None:
    reads = [
        "AAACCCCC",
        "AAAACCCC",
        "AAACCCCC",
        "AAACCCCC",
        "AAACCCCCC",
    ]
    kwargs = {
        "target_length": 8,
        "top_runs": 2,
        "beam_width": 6,
        "depth": 1,
        "min_gap_support": 1,
        "homopolymer_weight": 0.15,
        "qgram_width": 4,
        "qgram_weight": 0.25,
        "minimum_score_gain": 0.05,
    }
    assert homopolymer_balance_trace_consensus(
        reads,
        **kwargs,
    ) == homopolymer_balance_trace_consensus(
        list(reversed(reads)),
        **kwargs,
    )


def test_homopolymer_balance_validates_configuration() -> None:
    reads = ["ACGT", "ACGA"]
    invalid = (
        {"target_length": 0},
        {"target_length": 4, "top_runs": 0},
        {"target_length": 4, "beam_width": 0},
        {"target_length": 4, "depth": 3},
        {"target_length": 4, "min_gap_support": -1},
        {"target_length": 4, "homopolymer_weight": -0.1},
        {"target_length": 4, "qgram_width": 0},
        {"target_length": 4, "qgram_weight": -0.1},
        {"target_length": 4, "minimum_score_gain": -0.1},
    )
    for kwargs in invalid:
        try:
            homopolymer_balance_trace_consensus(reads, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"configuration should fail: {kwargs}")


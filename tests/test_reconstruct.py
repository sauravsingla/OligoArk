from dataclasses import dataclass

import pytest

from oligoark.archive import ArchiveConfig, DNAArchive, archive_bytes, recover_bytes
from oligoark.reconstruct import (
    GraphConsensusReconstructor,
    LevenshteinEdgeScorer,
    ReconstructionResult,
    alignment_consensus,
    build_similarity_graph,
    confidence_fusion_trace_consensus,
    edit_distance,
    global_align,
    graph_cluster_consensus,
    multistart_trace_consensus,
    normalized_similarity,
    targeted_trace_consensus,
)
from oligoark.simulator import SimulationConfig, simulate_channel


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


@dataclass(frozen=True)
class SyntheticGraphBaselineFixture:
    """Tiny reproducible synthetic read fixture for graph reconstruction baseline (#43, #18)."""

    payload: bytes
    archive: DNAArchive
    ground_truth_strands: tuple[str, ...]
    config: SimulationConfig
    reads: tuple[str, ...]


def generate_synthetic_graph_baseline_fixture() -> SyntheticGraphBaselineFixture:
    """Generate a tiny deterministic synthetic-read fixture.

    Parameters:
    - Payload: 69 bytes synthetic ASCII text, creating 3 distinct archival strands
      with standard Reed-Solomon framing and zero inter-strand cross-clustering.
    - SimulationConfig:
      - substitution_rate = 0.015 (1.5% substitution noise)
      - insertion_rate = 0.005 (0.5% insertion noise)
      - deletion_rate = 0.005 (0.5% deletion noise)
      - dropout_rate = 0.0 (0.0% strand loss; all strands retained to provide trace evidence)
      - duplicate_rate = 0.0 (0.0% duplicates)
      - copies_per_strand = 5 (5 simulated reads per strand, total 15 reads)
      - seed = 42 (fixed random seed for full determinism)
    """
    payload = b"OligoArk reproducible graph-reconstruction baseline fixture (#43, #18)"
    archive = archive_bytes(
        payload,
        ArchiveConfig(chunk_size=32, rs_nsym=8, redundancy_scheme="none"),
    )
    ground_truth = tuple(archive.strands)
    config = SimulationConfig(
        substitution_rate=0.015,
        insertion_rate=0.005,
        deletion_rate=0.005,
        dropout_rate=0.0,
        duplicate_rate=0.0,
        copies_per_strand=5,
        seed=42,
    )
    reads = tuple(simulate_channel(list(ground_truth), config))
    return SyntheticGraphBaselineFixture(
        payload=payload,
        archive=archive,
        ground_truth_strands=ground_truth,
        config=config,
        reads=reads,
    )


@pytest.fixture
def synthetic_graph_baseline() -> SyntheticGraphBaselineFixture:
    """Pytest fixture providing the reproducible graph reconstruction baseline dataset."""
    return generate_synthetic_graph_baseline_fixture()


def test_graph_reconstruction_baseline_fixture(
    synthetic_graph_baseline: SyntheticGraphBaselineFixture,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Deterministic graph reconstruction baseline verification (#43, #18).

    Verifies:
    1. Structural graph invariants: nodes (15), edges (30), and connected components (3).
    2. Component clustering: 3 clusters of size 5, corresponding to the 3 ground-truth strands.
    3. Consensus recovery: 100% (3/3) exact match to synthetic ground-truth strands
       when using alignment consensus. Ground truth is used strictly for post-hoc scoring.
    4. Recovery distinction:
       - Direct frame decoding without graph reconstruction fails (False).
       - Alignment-consensus-assisted recovery achieves SHA-256 verified recovery (True).
       - Medoid non-alignment consensus fails to verify (False), proving the value
         of alignment indel repair.
    """
    fixture = synthetic_graph_baseline

    # Reconstruction strictly takes reads without access to ground truth or archive
    reconstructor = GraphConsensusReconstructor(
        threshold=0.90,
        scorer=LevenshteinEdgeScorer(),
        consensus_mode="alignment",
    )
    result = reconstructor.reconstruct(list(fixture.reads))

    # 1. Structural graph invariants
    assert result.node_count == len(fixture.reads) == 15
    assert result.edge_count == 30
    assert result.component_count == len(fixture.ground_truth_strands) == 3
    assert sorted(result.cluster_sizes) == [5, 5, 5]
    assert sum(result.cluster_sizes) == result.node_count

    # 2. Sequence-level consensus recovery (post-hoc evaluation only)
    ground_truth_set = set(fixture.ground_truth_strands)
    consensus_set = set(result.consensus_reads)
    exact_matches = len(consensus_set & ground_truth_set)
    assert exact_matches == len(ground_truth_set) == 3
    assert consensus_set == ground_truth_set

    # 3. Direct archive recovery fails without reconstruction due to noise
    direct_recovered = False
    try:
        direct_recovered = (
            recover_bytes(fixture.archive, list(fixture.reads)) == fixture.payload
        )
    except (ValueError, Exception):
        direct_recovered = False
    assert direct_recovered is False

    # 4. Graph-assisted recovery succeeds with verified SHA-256
    augmented_reads = list(fixture.reads) + result.consensus_reads
    archive_recovered = (
        recover_bytes(fixture.archive, augmented_reads) == fixture.payload
    )
    assert archive_recovered is True

    # 5. Medoid consensus ablation contrast (non-alignment fails to repair indels)
    medoid_reconstructor = GraphConsensusReconstructor(
        threshold=0.90,
        scorer=LevenshteinEdgeScorer(),
        consensus_mode="medoid",
    )
    medoid_result = medoid_reconstructor.reconstruct(list(fixture.reads))
    medoid_matches = len(set(medoid_result.consensus_reads) & ground_truth_set)
    assert medoid_matches < len(ground_truth_set)
    medoid_recovered = False
    try:
        medoid_recovered = (
            recover_bytes(
                fixture.archive,
                list(fixture.reads) + medoid_result.consensus_reads,
            )
            == fixture.payload
        )
    except (ValueError, Exception):
        medoid_recovered = False
    assert medoid_recovered is False

    # Report baseline summary for plain-language transparency when running with -s
    report = (
        "\n"
        "================ GRAPH RECONSTRUCTION BASELINE REPORT ================\n"
        f"Input nodes (reads):        {result.node_count}\n"
        f"Candidate pairs evaluated:  {result.candidate_pairs}\n"
        f"Filtered graph edges:       {result.edge_count}\n"
        f"Connected components:       {result.component_count}\n"
        f"Component cluster sizes:    {result.cluster_sizes}\n"
        f"Consensus reads generated:  {len(result.consensus_reads)}\n"
        f"Ground-truth exact match:   {exact_matches} / {len(ground_truth_set)} (100.0%)\n"
        f"Direct recovery (no graph): {direct_recovered}\n"
        f"Graph-assisted recovery:    {archive_recovered} (SHA-256 verified)\n"
        f"Medoid ablation recovery:   {medoid_recovered} (fails due to unaligned indels)\n"
        "======================================================================\n"
    )
    print(report)


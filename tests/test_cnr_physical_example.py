import importlib.util
import json
from pathlib import Path
from types import ModuleType

import pytest

from oligoark.physical import (
    PhysicalDatasetManifest,
    evaluate_physical_reconstruction,
    evaluate_supplied_clusters,
    read_sequences,
)
from oligoark.reconstruct import normalized_similarity

ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "benchmarks"
FIXTURE = Path(__file__).resolve().parent / "fixtures" / "cnr_tiny"
MANIFEST = ROOT / "datasets" / "cnr.json"


def _load_converter(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    monkeypatch.syspath_prepend(str(BENCHMARKS))
    spec = importlib.util.spec_from_file_location(
        "convert_cnr_to_physical", BENCHMARKS / "convert_cnr_to_physical.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fixture(converter: ModuleType) -> tuple[list[str], list[list[str]]]:
    centers = [
        line.strip()
        for line in (FIXTURE / "Centers.txt").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    clusters = converter.load_clusters(FIXTURE / "Clusters.txt", expected_count=len(centers))
    return centers, clusters


def test_cnr_manifest_records_provenance_checksums_and_limitation() -> None:
    manifest = PhysicalDatasetManifest.load(MANIFEST)
    assert manifest.doi == "10.48550/arXiv.2107.06440"
    assert manifest.run_accessions == (
        "microsoft/clustered-nanopore-reads-dataset@6938f44796185902a08381943c2895782886c5c3",
    )
    assert "MIT" in manifest.data_restrictions
    assert "nearest reference" in manifest.reference_mapping_status
    assert "not the authoritative" in manifest.reference_mapping_status

    raw = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert "long-range dependencies" in raw["known_limitations"]
    assert "malformed" in raw["known_limitations"]
    files = raw["input_files"]
    assert files["Centers.txt"]["git_blob"] == "2d72f78cd86f2172578aaae297113a4baa14975c"
    assert files["Clusters.txt"]["git_blob"] == "95a71c77cb82f6d501d6aa99c7415591565c5ae3"
    for entry in files.values():
        assert len(entry["sha256"]) == 64


def test_fixture_selection_keeps_cnr_association(monkeypatch: pytest.MonkeyPatch) -> None:
    converter = _load_converter(monkeypatch)
    centers, clusters = _fixture(converter)
    assert [len(reads) for reads in clusters] == [3, 3, 0, 3]

    selected = converter.select_clusters(centers, clusters, limit=10, max_reads_per_cluster=3)

    assert [index for index, _, _ in selected] == [0, 1, 3]
    for index, center, reads in selected:
        assert center == centers[index]
        assert reads == clusters[index]


def test_fixture_nearest_reference_agrees_with_cnr_clusters(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    converter = _load_converter(monkeypatch)
    centers, clusters = _fixture(converter)
    selected = converter.select_clusters(centers, clusters, limit=10, max_reads_per_cluster=3)
    references = [center for _, center, _ in selected]

    for expected, (_, _, reads) in enumerate(selected):
        for read in reads:
            scores = [normalized_similarity(reference, read) for reference in references]
            assert scores.index(max(scores)) == expected
            assert max(scores) >= 0.70


def test_fixture_both_evaluation_paths(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    converter = _load_converter(monkeypatch)
    centers, clusters = _fixture(converter)
    selected = converter.select_clusters(centers, clusters, limit=10, max_reads_per_cluster=2)
    manifest = PhysicalDatasetManifest.load(MANIFEST)

    supplied = converter.evaluate_selection(manifest, selected)
    assert supplied["assignment"] == "supplied_clusters"
    assert supplied["cluster_indices"] == [0, 1, 3]
    assert supplied["reference_count"] == 3
    assert supplied["total_reads"] == supplied["assigned_reads"] == 6
    assert supplied["unassigned_reads"] == 0

    reads_out = tmp_path / "reads.fasta"
    references_out = tmp_path / "references.fasta"
    assert converter.write_physical_inputs(
        selected, reads_out=reads_out, references_out=references_out
    ) == (3, 6)
    assert ">cnr_00002\n" not in references_out.read_text(encoding="utf-8")
    assert read_sequences(references_out) == [centers[0], centers[1], centers[3]]

    nearest = evaluate_physical_reconstruction(
        manifest, read_sequences(reads_out), read_sequences(references_out)
    ).to_dict()
    assert nearest["assignment"] == "nearest_reference"
    assert nearest["assigned_reads"] == 6
    assert nearest["unassigned_reads"] == 0
    # Two raw nanopore reads per cluster are too noisy for any baseline to recover a
    # 110-base center exactly, so both paths report zero exact reconstructions here.
    for result in (supplied, nearest):
        assert result["exact_single_read_reference_matches"] == 0
        assert result["medoid_exact_reference_matches"] == 0
        assert result["alignment_exact_reference_matches"] == 0
        assert result["trace_exact_reference_matches"] == 0


def test_fixture_empty_cluster_is_scored_as_a_miss(monkeypatch: pytest.MonkeyPatch) -> None:
    converter = _load_converter(monkeypatch)
    centers, clusters = _fixture(converter)
    manifest = PhysicalDatasetManifest.load(MANIFEST)

    result = evaluate_supplied_clusters(manifest, centers, clusters)

    assert result.reference_count == 4
    assert result.assigned_reads == result.total_reads == 9
    assert result.clusters_with_multiple_reads == 3
    assert result.trace_exact_reference_matches == 0
    assert result.trace_exact_rate == 0.0


def test_supplied_clusters_recover_a_known_reference() -> None:
    manifest = PhysicalDatasetManifest.load(MANIFEST)
    reference = "ACGTTGCAAGCTTCGAGGATCCATGCAGTCGACTTAAGCTAGCATCG"
    # Each read carries one different error, so no read equals the reference but the
    # consensus does.
    reads = [
        reference[:5] + "A" + reference[6:],
        reference[:20] + "A" + reference[21:],
        reference[:35] + "C" + reference[36:],
        reference[:12] + reference[13:],
    ]
    assert reference not in reads

    result = evaluate_supplied_clusters(manifest, [reference], [reads])

    assert result.exact_single_read_reference_matches == 0
    assert result.alignment_exact_reference_matches == 1
    assert result.trace_exact_reference_matches == 1
    assert result.trace_exact_rate == 1.0


def test_selection_stops_at_limit_and_validates(monkeypatch: pytest.MonkeyPatch) -> None:
    converter = _load_converter(monkeypatch)
    selected = converter.select_clusters(
        ["AAAA", "CCCC", "GGGG"],
        [["AAAA"], [], ["GGGG", "GGGA"]],
        limit=1,
        max_reads_per_cluster=5,
    )
    assert selected == [(0, "AAAA", ["AAAA"])]
    with pytest.raises(ValueError):
        converter.select_clusters(["AAAA"], [], limit=1, max_reads_per_cluster=1)
    with pytest.raises(ValueError):
        converter.select_clusters(["AAAA"], [["AAAA"]], limit=0, max_reads_per_cluster=1)


def test_verify_inputs_rejects_files_that_differ_from_manifest(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    converter = _load_converter(monkeypatch)
    with pytest.raises(ValueError, match="Centers.txt SHA-256"):
        converter.verify_inputs(MANIFEST, FIXTURE / "Centers.txt", FIXTURE / "Clusters.txt")


def test_verify_inputs_accepts_matching_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    converter = _load_converter(monkeypatch)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "input_files": {
                    "Centers.txt": {"sha256": converter.sha256_file(FIXTURE / "Centers.txt")},
                    "Clusters.txt": {"sha256": converter.sha256_file(FIXTURE / "Clusters.txt")},
                }
            }
        ),
        encoding="utf-8",
    )
    hashes = converter.verify_inputs(manifest, FIXTURE / "Centers.txt", FIXTURE / "Clusters.txt")
    assert set(hashes) == {"Centers.txt", "Clusters.txt"}


def test_supplied_clusters_are_not_reassigned() -> None:
    manifest = PhysicalDatasetManifest.load(MANIFEST)
    first = "ACGT" * 10
    second = "TTGCA" * 8
    # The second cluster's read is closer to the first reference, but the dataset
    # says it belongs to the second one, so it must stay there.
    result = evaluate_supplied_clusters(manifest, [first, second], [[first, first], [first]])
    assert result.assignment == "supplied_clusters"
    assert result.assignment_threshold is None
    assert result.assigned_reads == 3
    assert result.exact_single_read_reference_matches == 2
    assert result.trace_exact_reference_matches == 1

    nearest = evaluate_physical_reconstruction(manifest, [first, first, first], [first, second])
    # Nearest-reference assignment moves the third read to the first reference, where
    # it counts as an exact single-read match; the supplied path keeps it in place.
    assert nearest.exact_single_read_reference_matches == 3
    assert nearest.trace_exact_reference_matches == 1
    assert nearest.clusters_with_multiple_reads == 1


def test_supplied_clusters_validate_input() -> None:
    manifest = PhysicalDatasetManifest.load(MANIFEST)
    with pytest.raises(ValueError, match="same length"):
        evaluate_supplied_clusters(manifest, ["ACGT"], [])
    with pytest.raises(ValueError, match="must not be empty"):
        evaluate_supplied_clusters(manifest, [], [])
    with pytest.raises(ValueError, match="empty sequence"):
        evaluate_supplied_clusters(manifest, [" "], [["ACGT"]])

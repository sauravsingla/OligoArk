import json
from pathlib import Path

from oligoark.physical import (
    PhysicalDatasetManifest,
    evaluate_physical_reconstruction,
    read_sequences,
)


def test_dna_aeon_manifest_records_public_provenance() -> None:
    manifest = PhysicalDatasetManifest.load("datasets/dna_aeon.json")
    assert manifest.doi == "10.1038/s41467-023-36297-3"
    assert manifest.bioproject == "PRJNA855029"
    assert "SRR19954693" in manifest.run_accessions
    assert "reference" in manifest.reference_mapping_status.lower()


def test_fasta_fastq_parsing_and_reference_reconstruction(tmp_path: Path) -> None:
    reference = "ACGT" * 30
    fasta = tmp_path / "refs.fasta"
    fasta.write_text(f">ref\n{reference}\n", encoding="utf-8")
    fastq = tmp_path / "reads.fastq"
    reads = [
        reference[:20] + "A" + reference[20:],
        reference[:44] + reference[45:],
        reference[:70] + "C" + reference[70:],
        reference,
    ]
    fastq.write_text(
        "".join(
            f"@r{index}\n{read}\n+\n{'I' * len(read)}\n"
            for index, read in enumerate(reads)
        ),
        encoding="utf-8",
    )
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "dataset": "synthetic physical-adapter fixture",
                "doi": "test",
                "bioproject": "test",
                "run_accessions": ["TEST1"],
                "data_restrictions": "fixture",
                "reference_mapping_status": "explicit fixture mapping",
                "notes": "test only",
            }
        ),
        encoding="utf-8",
    )
    manifest = PhysicalDatasetManifest.load(manifest_path)
    parsed_reads = read_sequences(fastq)
    parsed_refs = read_sequences(fasta)
    result = evaluate_physical_reconstruction(
        manifest,
        parsed_reads,
        parsed_refs,
        assignment_threshold=0.8,
    )
    assert result.assigned_reads == len(reads)
    assert result.reference_count == 1
    assert result.trace_exact_reference_matches == 1
    assert result.trace_exact_rate == 1.0

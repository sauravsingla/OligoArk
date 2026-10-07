from __future__ import annotations

import json

from oligoark.wetlab import prepare_wetlab_bundle, recover_wetlab_reads


def test_wetlab_bundle_is_self_describing_and_clean_roundtrip(tmp_path) -> None:
    payload = (
        b"OligoArk wet-lab handoff smoke test\n"
        + bytes(range(128))
        + b"exact SHA-256 recovery"
    )
    source = tmp_path / "pilot.bin"
    bundle = tmp_path / "bundle"
    recovered = tmp_path / "recovered.bin"
    source.write_bytes(payload)

    summary = prepare_wetlab_bundle(
        source,
        bundle,
        profile_name="oligoark-200-compact",
        redundancy_scheme="hybrid",
        fountain_redundancy=0.125,
    )

    assert summary.claim_status == "prepared-not-executed"
    assert summary.physical_execution_status == "prepared, not physically executed"
    assert summary.source_bytes == len(payload)
    assert summary.strand_count > summary.data_strands
    assert (bundle / "archive.json").exists()
    assert (bundle / "oligos.fasta").exists()
    assert (bundle / "oligos.csv").exists()
    assert (bundle / "provider-metadata.template.json").exists()
    assert (bundle / "preprocessing-record.template.json").exists()
    assert (bundle / "read-depth-plan.json").exists()
    assert (bundle / "bundle-checksums.sha256").exists()

    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "prepared-not-executed"
    assert manifest["physical_execution_status"] == "prepared, not physically executed"
    assert manifest["success_criterion"].startswith("Recovered payload SHA-256")
    assert manifest["planned_read_depths_per_strand"] == [1, 5, 10, 20]

    result = recover_wetlab_reads(
        bundle,
        bundle / "oligos.fasta",
        recovered,
    )
    assert result.verified_sha256 is True
    assert len(result.input_reads_sha256) == 64
    assert result.claim_status == "archive-recovery-sha256-verified"
    assert recovered.read_bytes() == payload


def test_wetlab_manifest_does_not_claim_physical_execution(tmp_path) -> None:
    source = tmp_path / "pilot.bin"
    source.write_bytes(b"claim boundary")
    bundle = tmp_path / "bundle"
    prepare_wetlab_bundle(source, bundle)

    manifest_text = (bundle / "manifest.json").read_text(encoding="utf-8")
    assert "prepared-not-executed" in manifest_text
    assert "prepared, not physically executed" in manifest_text
    assert "End-to-end physical-storage evidence requires" in manifest_text


def test_efficient_152_profile_can_prepare_synthesis_bundle(tmp_path) -> None:
    source = tmp_path / "pilot.bin"
    source.write_bytes(b"efficient wet-lab pilot" * 32)
    bundle = tmp_path / "bundle"

    summary = prepare_wetlab_bundle(
        source,
        bundle,
        profile_name="oligoark-152-efficient-v1",
        redundancy_scheme="hybrid",
        fountain_redundancy=0.125,
    )

    assert summary.target_strand_nt == 152
    assert summary.profile == "oligoark-152-efficient-v1"
    assert summary.physical_execution_status == "prepared, not physically executed"
    assert max(
        len(line.strip())
        for line in (bundle / "oligos.fasta").read_text(encoding="utf-8").splitlines()
        if line and not line.startswith(">")
    ) <= 152

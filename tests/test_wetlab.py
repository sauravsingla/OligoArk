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
        profile_name="oligoark-200",
        redundancy_scheme="hybrid",
        fountain_redundancy=0.125,
    )

    assert summary.claim_status == "prepared-not-executed"
    assert summary.source_bytes == len(payload)
    assert summary.strand_count > summary.data_strands
    assert (bundle / "archive.json").exists()
    assert (bundle / "oligos.fasta").exists()
    assert (bundle / "oligos.csv").exists()

    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "prepared-not-executed"
    assert manifest["success_criterion"].startswith("Recovered payload SHA-256")
    assert manifest["planned_read_depths_per_strand"] == [1, 5, 10, 20]

    result = recover_wetlab_reads(
        bundle,
        bundle / "oligos.fasta",
        recovered,
    )
    assert result.verified_sha256 is True
    assert result.claim_status == "physical-read-archive-recovery-verified"
    assert recovered.read_bytes() == payload


def test_wetlab_manifest_does_not_claim_physical_execution(tmp_path) -> None:
    source = tmp_path / "pilot.bin"
    source.write_bytes(b"claim boundary")
    bundle = tmp_path / "bundle"
    prepare_wetlab_bundle(source, bundle)

    manifest_text = (bundle / "manifest.json").read_text(encoding="utf-8")
    assert "prepared-not-executed" in manifest_text
    assert "End-to-end physical-storage evidence requires" in manifest_text

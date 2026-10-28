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

    assert summary.claim_status == "prepared, not physically executed"
    assert summary.source_bytes == len(payload)
    assert summary.strand_count > summary.data_strands
    assert (bundle / "archive.json").exists()
    assert (bundle / "oligos.fasta").exists()
    assert (bundle / "oligos.csv").exists()

    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "prepared, not physically executed"
    assert manifest["success_criterion"].startswith("Recovered payload SHA-256")
    assert manifest["planned_read_depths_per_strand"] == [1, 5, 10, 20]

    result = recover_wetlab_reads(
        bundle,
        bundle / "oligos.fasta",
        recovered,
    )
    assert result.verified_sha256 is True
    assert result.claim_status == "archive-recovery-sha256-verified"
    assert recovered.read_bytes() == payload


def test_wetlab_manifest_does_not_claim_physical_execution(tmp_path) -> None:
    source = tmp_path / "pilot.bin"
    source.write_bytes(b"claim boundary")
    bundle = tmp_path / "bundle"
    prepare_wetlab_bundle(source, bundle)

    manifest_text = (bundle / "manifest.json").read_text(encoding="utf-8")
    assert "prepared, not physically executed" in manifest_text
    assert "End-to-end physical-storage evidence requires" in manifest_text


def test_recovery_rejects_inconsistent_preorder_source_hash(tmp_path) -> None:
    source = tmp_path / "pilot.bin"
    source.write_bytes(b"frozen pilot identity must match decoder archive")
    bundle = tmp_path / "bundle"
    prepare_wetlab_bundle(source, bundle)

    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["source"]["sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    from pytest import raises

    with raises(ValueError, match="manifest and archive SHA-256 disagree"):
        recover_wetlab_reads(bundle, bundle / "oligos.fasta", tmp_path / "invalid.bin")
    assert not (tmp_path / "invalid.bin").exists()


def test_failed_checksum_preserves_existing_recovery_output(tmp_path, monkeypatch) -> None:
    import pytest

    import oligoark.wetlab as wetlab_module

    source = tmp_path / "source.bin"
    source.write_bytes(b"physical pilot must protect prior recovery")
    bundle = tmp_path / "bundle"
    prepare_wetlab_bundle(source, bundle)

    output = tmp_path / "verified.bin"
    output.write_bytes(b"previously verified recovery")

    def corrupt_recovery(*_args, **_kwargs):
        return b"corrupt reconstruction", object()

    monkeypatch.setattr(wetlab_module, "recover_from_reads", corrupt_recovery)
    with pytest.raises(ValueError, match="failed SHA-256 verification"):
        recover_wetlab_reads(bundle, bundle / "oligos.fasta", output)
    assert output.read_bytes() == b"previously verified recovery"
    assert not list(tmp_path.glob(".verified.bin.*.part"))
    assert not (bundle / "recovery-report.json").exists()


def test_wetlab_output_must_not_destroy_sequence_input(tmp_path) -> None:
    import pytest

    source = tmp_path / "source.bin"
    source.write_bytes(b"preserve sequencing evidence")
    bundle = tmp_path / "bundle"
    prepare_wetlab_bundle(source, bundle)
    reads_path = bundle / "oligos.fasta"
    original_reads = reads_path.read_bytes()

    with pytest.raises(ValueError, match="distinct from evidence inputs"):
        recover_wetlab_reads(bundle, reads_path, reads_path)
    assert reads_path.read_bytes() == original_reads
    assert not list(bundle.glob(".oligos.fasta.*.part"))


def test_successful_wetlab_recovery_replaces_prior_output(tmp_path) -> None:
    source = tmp_path / "source.bin"
    source.write_bytes(b"successful archive recovery" * 2)
    bundle = tmp_path / "bundle"
    prepare_wetlab_bundle(source, bundle)
    output = tmp_path / "recovered.bin"
    output.write_bytes(b"outdated content")

    result = recover_wetlab_reads(bundle, bundle / "oligos.fasta", output)
    assert result.verified_sha256
    assert output.read_bytes() == source.read_bytes()
    assert not list(tmp_path.glob(".recovered.bin.*.part"))

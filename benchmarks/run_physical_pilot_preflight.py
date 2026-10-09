"""Deterministic, non-biological, software-only preflight for a DNA synthesis pilot.

Creates a reproducible 4 KiB payload, prepares OligoArk core oligos, validates
the synthesized-core design constraints, and performs a pristine software round
trip. It does NOT demonstrate or assert physical DNA synthesis or sequencing.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from oligoark.wetlab import prepare_wetlab_bundle, recover_wetlab_reads


def pilot_payload(length: int = 4096) -> bytes:
    """Fixed diverse non-biological bytes; no secrets, human data or sample material."""
    if length <= 0:
        raise ValueError("payload length must be positive")
    prefix = b"OligoArk physical pilot v1 - synthetic nonbiological data\n"
    chunks = [prefix, bytes(range(256))]
    i = 0
    while sum(map(len, chunks)) < length:
        chunks.append(hashlib.sha256(b"oligoark-pilot-v1:" + i.to_bytes(8, "big")).digest())
        i += 1
    return b"".join(chunks)[:length]


def preflight(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    source = output_dir / "pilot.bin"
    data = pilot_payload()
    source.write_bytes(data)
    frozen_sha = hashlib.sha256(data).hexdigest()

    bundle = output_dir / "bundle"
    summary = prepare_wetlab_bundle(
        source, bundle, profile_name="oligoark-152-compact-v3",
        redundancy_scheme="hybrid", fountain_redundancy=0.125,
    )
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["source"]["sha256"] == frozen_sha, "source hash mismatch"
    assert summary.source_sha256 == frozen_sha, "summary hash mismatch"
    assert manifest["status"] == "prepared, not physically executed"

    with (bundle / "oligos.csv").open(newline="", encoding="utf-8") as handle:
        oligos = list(csv.DictReader(handle))
    assert len(oligos) == summary.strand_count > 0, "strand count mismatch"
    ids = [row["oligo_id"] for row in oligos]
    sequences = [row["sequence"] for row in oligos]
    assert len(set(ids)) == len(ids), "duplicate oligo IDs"
    assert len(set(sequences)) == len(sequences), "duplicate oligo sequences"
    assert all(set(seq) <= set("ACGT") for seq in sequences), "invalid bases"
    assert all(0 < len(seq) <= 152 for seq in sequences), "core length exceeded"

    cfg = manifest["codec"]["config"]
    for row in oligos:
        seq = row["sequence"]
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        assert cfg["min_gc_fraction"] <= gc <= cfg["max_gc_fraction"], "GC constraint failure"
        assert int(row["length_nt"]) == len(seq), "CSV length mismatch"
        homopolymer = 1
        current = 1
        for prev, base in zip(seq, seq[1:], strict=False):
            current = current + 1 if prev == base else 1
            homopolymer = max(homopolymer, current)
        assert homopolymer <= cfg["max_homopolymer"], "homopolymer constraint failure"
    with (bundle / "oligos.fasta").open(encoding="utf-8") as handle:
        fasta = [line.strip() for line in handle if line.strip() and not line.startswith(">")]
    assert fasta == sequences, "FASTA/CSV disagreement"

    # Designed core strands are used as synthetic perfect reads: NOT lab-generated.
    recovered = output_dir / "recovered-software.bin"
    outcome = recover_wetlab_reads(bundle, bundle / "oligos.fasta", recovered)
    assert outcome.verified_sha256 and recovered.read_bytes() == data
    assert outcome.recovered_sha256 == frozen_sha

    report: dict[str, object] = {
        "schema_version": 1,
        "experiment": "oligoark-physical-pilot-v1",
        "status": "software-preflight-passed; physical work not performed",
        "source_bytes": len(data),
        "original_source_sha256": frozen_sha,
        "profile": "oligoark-152-compact-v3",
        "scheme": "hybrid",
        "total_core_oligos": len(oligos),
        "max_core_nt": max(map(len, sequences)),
        "max_planned_total_nt_with_2x20nt_flanks": max(map(len, sequences)) + 40,
        "software_pristine_roundtrip_sha256_passed": True,
        "planned_physical_replicates": 2,
        "planned_read_depths_per_strand": [1, 5, 10, 20],
        "not_yet_done": [
            "supplier sequence screening and confirmed primer/adapter design",
            "synthesis purchase and documented lot QC",
            "independent aliquot preparations, positive/negative controls",
            "physical sequencing with authentic FASTQ/FASTA and provider run records",
            "independent recovery from genuine sequencing reads and SHA-256 gate",
        ],
        "claim_scope": (
            "No wet-lab success claim: reads here were generated from "
            "designed core strands."
        ),
    }
    (output_dir / "preflight-report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("physical-pilot-results"))
    args = parser.parse_args()
    print(json.dumps(preflight(args.output), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()

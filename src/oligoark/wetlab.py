"""Wet-lab preparation and recovery helpers.

This module prepares an OligoArk archive for a laboratory handoff and verifies sequencing
outputs against the original SHA-256 digest. Preparing a bundle is not a wet-lab result:
physical synthesis, storage, sequencing, and final SHA-256 verification must all be completed
before an end-to-end physical-storage claim is made.
"""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import cast

from .archive import DNAArchive, archive_bytes, archive_statistics, recover_from_reads
from .dna import sequence_metrics
from .framing import decode_frame
from .physical import read_sequences
from .profiles import physical_strand_profile
from .reconstruct import TraceConsensusReconstructor


@dataclass(frozen=True)
class WetLabBundleSummary:
    source_bytes: int
    source_sha256: str
    profile: str
    strand_count: int
    data_strands: int
    parity_strands: int
    fountain_strands: int
    encoded_nucleotides: int
    logical_bits_per_nucleotide: float
    target_strand_nt: int
    claim_status: str = "prepared-not-executed"

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class WetLabRecoverySummary:
    input_reads: int
    source_bytes: int
    expected_sha256: str
    recovered_sha256: str
    verified_sha256: bool
    reconstruction_strategy: str | None
    consensus_reads: int
    claim_status: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_fasta(strands: list[str], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for ordinal, sequence in enumerate(strands):
            handle.write(f">OA-{ordinal:08d}\n{sequence}\n")


def _write_oligo_csv(archive: DNAArchive, path: Path) -> None:
    config_obj = archive.metadata["config"]
    if not isinstance(config_obj, dict):
        raise ValueError("archive config metadata must be an object")
    total = int(cast(int, archive.metadata["data_strands"]))
    config = cast(dict[str, object], config_obj)

    rows: list[dict[str, object]] = []
    for ordinal, sequence in enumerate(archive.strands):
        frame = decode_frame(
            sequence,
            rs_nsym=int(cast(int, config["rs_nsym"])),
            mask_search_limit=int(cast(int, config["mask_search_limit"])),
            compact_framing=bool(config.get("compact_framing", False)),
            compact_index_bytes=int(cast(int, config.get("compact_index_bytes", 3))),
            expected_total_data=total if bool(config.get("compact_framing", False)) else None,
        )
        metrics = sequence_metrics(sequence)
        rows.append(
            {
                "oligo_id": f"OA-{ordinal:08d}",
                "ordinal": ordinal,
                "frame_kind": frame.kind,
                "frame_index": frame.index,
                "length_nt": len(sequence),
                "gc_fraction": round(metrics.gc_fraction, 6),
                "max_homopolymer": metrics.max_homopolymer,
                "sequence": sequence,
            }
        )

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "oligo_id",
                "ordinal",
                "frame_kind",
                "frame_index",
                "length_nt",
                "gc_fraction",
                "max_homopolymer",
                "sequence",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def prepare_wetlab_bundle(
    source_path: str | Path,
    output_dir: str | Path,
    *,
    profile_name: str = "oligoark-200",
    redundancy_scheme: str = "hybrid",
    fountain_redundancy: float = 0.125,
) -> WetLabBundleSummary:
    """Create synthesis-ready core oligos plus all metadata required for recovery."""
    source = Path(source_path)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    payload = source.read_bytes()

    profile = physical_strand_profile(profile_name).with_scheme(redundancy_scheme)
    config = profile.to_archive_config(fountain_redundancy=fountain_redundancy)
    archive = archive_bytes(payload, config)
    stats = archive_statistics(archive)

    archive.save(output / "archive.json")
    _write_fasta(archive.strands, output / "oligos.fasta")
    _write_oligo_csv(archive, output / "oligos.csv")

    summary = WetLabBundleSummary(
        source_bytes=len(payload),
        source_sha256=_sha256(payload),
        profile=profile_name,
        strand_count=stats.strand_count,
        data_strands=stats.data_strands,
        parity_strands=stats.parity_strands,
        fountain_strands=stats.fountain_strands,
        encoded_nucleotides=stats.encoded_nucleotides,
        logical_bits_per_nucleotide=stats.logical_bits_per_nucleotide,
        target_strand_nt=profile.target_nucleotides,
    )

    manifest: dict[str, object] = {
        "schema_version": 1,
        "experiment_id": "replace-with-lab-experiment-id",
        "status": summary.claim_status,
        "claim_scope": (
            "OligoArk oligo design and recovery package only. End-to-end physical-storage "
            "evidence requires documented synthesis, storage, sequencing, preprocessing, "
            "and final SHA-256 verification of these OligoArk-generated strands."
        ),
        "source": {
            "file_name": source.name,
            "bytes": len(payload),
            "sha256": summary.source_sha256,
        },
        "codec": {
            "profile": profile_name,
            "target_strand_nt": profile.target_nucleotides,
            "config": asdict(config),
        },
        "outputs": {
            "archive": "archive.json",
            "synthesis_fasta": "oligos.fasta",
            "synthesis_csv": "oligos.csv",
        },
        "sequencing_input_contract": {
            "accepted_formats": ["FASTA", "FASTQ", "FASTA.gz", "FASTQ.gz"],
            "orientation": (
                "reads should be supplied in OligoArk strand orientation; document any "
                "adapter trimming, reverse-complement handling, and quality filtering"
            ),
            "required_record_content": "DNA sequence using A/C/G/T after preprocessing",
        },
        "planned_read_depths_per_strand": [1, 5, 10, 20],
        "controls": [
            "clean software round-trip before ordering",
            "synthesis vendor QC where available",
            "sequencing negative control with no OligoArk library",
            "known-sequence positive control",
            "at least one independently prepared technical replicate when feasible",
        ],
        "success_criterion": (
            "Recovered payload SHA-256 must exactly equal the source SHA-256. "
            "Partial strand reconstruction is not end-to-end archive success."
        ),
        "provider_note": (
            "Oligos in this bundle are codec core sequences. Any provider-specific primers, "
            "adapters, barcodes, or flanking sequences must be documented separately and "
            "removed during preprocessing before OligoArk recovery."
        ),
        "summary": summary.to_dict(),
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return summary


def recover_wetlab_reads(
    bundle_dir: str | Path,
    reads_path: str | Path,
    output_path: str | Path,
    *,
    similarity_threshold: float = 0.86,
    trace_rounds: int = 4,
) -> WetLabRecoverySummary:
    """Recover supplied sequencing reads and enforce the archive SHA-256 acceptance gate."""
    bundle = Path(bundle_dir)
    archive = DNAArchive.load(bundle / "archive.json")
    reads = read_sequences(reads_path)
    if not reads:
        raise ValueError("sequencing input contains no reads")

    reconstructor = TraceConsensusReconstructor(
        thresholds=(0.94, 0.90, similarity_threshold),
        rounds=trace_rounds,
    )
    recovered, report = recover_from_reads(
        archive,
        reads,
        similarity_threshold=similarity_threshold,
        reconstructor=reconstructor,
    )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(recovered)

    expected = str(archive.metadata["sha256"])
    actual = _sha256(recovered)
    verified = actual == expected
    if not verified:
        raise ValueError("wet-lab recovery output failed SHA-256 verification")

    summary = WetLabRecoverySummary(
        input_reads=len(reads),
        source_bytes=len(recovered),
        expected_sha256=expected,
        recovered_sha256=actual,
        verified_sha256=True,
        reconstruction_strategy=report.reconstruction_strategy,
        consensus_reads=report.consensus_reads,
        claim_status="physical-read-archive-recovery-verified",
    )
    (bundle / "recovery-report.json").write_text(
        json.dumps(summary.to_dict(), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return summary

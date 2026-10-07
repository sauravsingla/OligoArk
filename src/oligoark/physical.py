"""Adapters for reproducible evaluation on externally supplied physical DNA reads."""

from __future__ import annotations

import gzip
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TextIO, cast

from .reconstruct import (
    alignment_consensus,
    iterative_trace_consensus,
    medoid_consensus,
    normalized_similarity,
)


@dataclass(frozen=True)
class PhysicalDatasetManifest:
    dataset: str
    doi: str
    bioproject: str
    run_accessions: tuple[str, ...]
    data_restrictions: str
    reference_mapping_status: str
    notes: str

    @classmethod
    def load(cls, path: str | Path) -> PhysicalDatasetManifest:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(value, dict):
            raise ValueError("physical dataset manifest must be a JSON object")
        runs = value.get("run_accessions")
        if not isinstance(runs, list) or not all(isinstance(item, str) for item in runs):
            raise ValueError("run_accessions must be a list of strings")
        fields = {
            "dataset",
            "doi",
            "bioproject",
            "data_restrictions",
            "reference_mapping_status",
            "notes",
        }
        if any(not isinstance(value.get(name), str) for name in fields):
            raise ValueError("physical dataset manifest text fields must be strings")
        return cls(
            dataset=cast(str, value["dataset"]),
            doi=cast(str, value["doi"]),
            bioproject=cast(str, value["bioproject"]),
            run_accessions=tuple(cast(list[str], runs)),
            data_restrictions=cast(str, value["data_restrictions"]),
            reference_mapping_status=cast(str, value["reference_mapping_status"]),
            notes=cast(str, value["notes"]),
        )

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class PhysicalReconstructionSummary:
    dataset: str
    total_reads: int
    reference_count: int
    assigned_reads: int
    unassigned_reads: int
    clusters_with_multiple_reads: int
    exact_single_read_reference_matches: int
    medoid_exact_reference_matches: int
    alignment_exact_reference_matches: int
    trace_exact_reference_matches: int
    assignment_threshold: float | None
    assignment: str = "nearest_reference"

    @property
    def trace_exact_rate(self) -> float:
        return round(
            self.trace_exact_reference_matches / max(1, self.reference_count),
            6,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            **asdict(self),
            "trace_exact_rate": self.trace_exact_rate,
            "claim_scope": (
                "reference-reconstruction evaluation on supplied physical reads; "
                "not OligoArk archive decoding unless the supplied references are "
                "OligoArk-generated strands"
            ),
        }


def _open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_sequences(path: str | Path, *, max_sequences: int | None = None) -> list[str]:
    """Read FASTA or FASTQ sequences, optionally from gzip-compressed input."""
    if max_sequences is not None and max_sequences < 1:
        raise ValueError("max_sequences must be positive")
    resolved = Path(path)
    if not resolved.exists():
        raise ValueError(f"sequence file does not exist: {resolved}")

    sequences: list[str] = []
    with _open_text(resolved) as handle:
        first = handle.readline()
        if not first:
            return []
        handle.seek(0)
        if first.startswith(">"):
            current: list[str] = []
            for line in handle:
                line = line.strip()
                if line.startswith(">"):
                    if current:
                        sequences.append("".join(current).upper())
                        current = []
                        if max_sequences is not None and len(sequences) >= max_sequences:
                            break
                    continue
                current.append(line)
            if (
                current
                and (max_sequences is None or len(sequences) < max_sequences)
            ):
                sequences.append("".join(current).upper())
        elif first.startswith("@"):
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().strip().upper()
                plus = handle.readline()
                quality = handle.readline()
                if not sequence or not plus.startswith("+") or not quality:
                    raise ValueError("malformed FASTQ record")
                sequences.append(sequence)
                if max_sequences is not None and len(sequences) >= max_sequences:
                    break
        else:
            raise ValueError("sequence input must be FASTA or FASTQ")

    if any(not sequence for sequence in sequences):
        raise ValueError("sequence input contains an empty sequence")
    return sequences


def evaluate_physical_reconstruction(
    manifest: PhysicalDatasetManifest,
    reads: list[str],
    references: list[str],
    *,
    assignment_threshold: float = 0.70,
) -> PhysicalReconstructionSummary:
    """Assign reads to nearest supplied reference and compare consensus baselines.

    This is intentionally a reference-reconstruction benchmark. It does not infer that
    external DNA-storage references use OligoArk's frame/archive format.
    """
    if not 0 <= assignment_threshold <= 1:
        raise ValueError("assignment_threshold must be between 0 and 1")
    if not reads or not references:
        raise ValueError("reads and references must not be empty")

    normalized_reads = [read.strip().upper() for read in reads if read.strip()]
    normalized_references = [
        reference.strip().upper()
        for reference in references
        if reference.strip()
    ]
    clusters: list[list[str]] = [[] for _ in normalized_references]
    unassigned = 0
    exact_reads = 0

    for read in normalized_reads:
        scored = [
            (normalized_similarity(reference, read), index)
            for index, reference in enumerate(normalized_references)
        ]
        score, index = max(scored, key=lambda item: (item[0], -item[1]))
        if score < assignment_threshold:
            unassigned += 1
            continue
        clusters[index].append(read)
        if read == normalized_references[index]:
            exact_reads += 1

    return _score_clusters(
        manifest,
        normalized_references,
        clusters,
        total_reads=len(normalized_reads),
        unassigned_reads=unassigned,
        exact_single_reads=exact_reads,
        assignment_threshold=assignment_threshold,
        assignment="nearest_reference",
    )


def evaluate_supplied_clusters(
    manifest: PhysicalDatasetManifest,
    references: list[str],
    clusters: list[list[str]],
) -> PhysicalReconstructionSummary:
    """Compare consensus baselines on reads already grouped by reference.

    ``clusters[i]`` holds the reads the dataset itself associates with
    ``references[i]`` (for example CNR's cluster order), so no read is
    reassigned by similarity.
    """
    if len(references) != len(clusters):
        raise ValueError("references and clusters must have the same length")
    if not references:
        raise ValueError("references must not be empty")
    normalized_references = [reference.strip().upper() for reference in references]
    if any(not reference for reference in normalized_references):
        raise ValueError("references must not contain an empty sequence")
    normalized_clusters = [
        [read.strip().upper() for read in cluster if read.strip()]
        for cluster in clusters
    ]
    exact_reads = sum(
        read == reference
        for reference, cluster in zip(normalized_references, normalized_clusters, strict=True)
        for read in cluster
    )
    return _score_clusters(
        manifest,
        normalized_references,
        normalized_clusters,
        total_reads=sum(len(cluster) for cluster in normalized_clusters),
        unassigned_reads=0,
        exact_single_reads=exact_reads,
        assignment_threshold=None,
        assignment="supplied_clusters",
    )


def _score_clusters(
    manifest: PhysicalDatasetManifest,
    references: list[str],
    clusters: list[list[str]],
    *,
    total_reads: int,
    unassigned_reads: int,
    exact_single_reads: int,
    assignment_threshold: float | None,
    assignment: str,
) -> PhysicalReconstructionSummary:
    medoid_matches = 0
    alignment_matches = 0
    trace_matches = 0
    multi = 0
    for reference, cluster in zip(references, clusters, strict=True):
        if not cluster:
            continue
        if len(cluster) > 1:
            multi += 1
        medoid_matches += int(medoid_consensus(cluster) == reference)
        alignment_matches += int(alignment_consensus(cluster) == reference)
        trace_matches += int(iterative_trace_consensus(cluster, rounds=4) == reference)

    return PhysicalReconstructionSummary(
        dataset=manifest.dataset,
        total_reads=total_reads,
        reference_count=len(references),
        assigned_reads=sum(len(cluster) for cluster in clusters),
        unassigned_reads=unassigned_reads,
        clusters_with_multiple_reads=multi,
        exact_single_read_reference_matches=exact_single_reads,
        medoid_exact_reference_matches=medoid_matches,
        alignment_exact_reference_matches=alignment_matches,
        trace_exact_reference_matches=trace_matches,
        assignment_threshold=assignment_threshold,
        assignment=assignment,
    )

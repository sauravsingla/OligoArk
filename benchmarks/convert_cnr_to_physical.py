"""Prepare and evaluate a CNR subset for OligoArk's physical-read evaluation.

The selected clusters are written as reference and read FASTA files for
run_physical_dataset.py (the generic nearest-reference adapter). With
``--supplied-output`` the same subset is also scored using CNR's own
cluster-to-center association, which is the scientifically faithful reading
of this dataset.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from run_external_cnr_benchmark import load_centers, load_clusters

from oligoark.physical import PhysicalDatasetManifest, evaluate_supplied_clusters


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_inputs(manifest_path: Path, centers: Path, clusters: Path) -> dict[str, str]:
    """Check both inputs against the SHA-256 values recorded in the manifest."""
    recorded = json.loads(manifest_path.read_text(encoding="utf-8"))["input_files"]
    actual = {"Centers.txt": sha256_file(centers), "Clusters.txt": sha256_file(clusters)}
    for name, value in actual.items():
        expected = recorded[name]["sha256"]
        if value != expected:
            raise ValueError(f"{name} SHA-256 is {value}, manifest records {expected}")
    return actual


def select_clusters(
    centers: list[str],
    clusters: list[list[str]],
    *,
    limit: int,
    max_reads_per_cluster: int,
) -> list[tuple[int, str, list[str]]]:
    """Return (cluster index, center, reads) for the first ``limit`` non-empty clusters."""
    if len(centers) != len(clusters):
        raise ValueError("centers and clusters must have the same length")
    if limit < 1 or max_reads_per_cluster < 1:
        raise ValueError("limit and max_reads_per_cluster must be positive")
    selected: list[tuple[int, str, list[str]]] = []
    for index, (center, reads) in enumerate(zip(centers, clusters, strict=True)):
        if not reads:
            continue
        selected.append((index, center, reads[:max_reads_per_cluster]))
        if len(selected) == limit:
            break
    return selected


def write_physical_inputs(
    selected: list[tuple[int, str, list[str]]],
    *,
    reads_out: Path,
    references_out: Path,
) -> tuple[int, int]:
    """Write FASTA files for run_physical_dataset.py and return (references, reads)."""
    reference_lines: list[str] = []
    read_lines: list[str] = []
    for index, center, reads in selected:
        reference_lines.append(f">cnr_{index:05d}\n{center}\n")
        for read_index, read in enumerate(reads):
            read_lines.append(f">cnr_{index:05d}_r{read_index}\n{read}\n")
    references_out.write_text("".join(reference_lines), encoding="utf-8")
    reads_out.write_text("".join(read_lines), encoding="utf-8")
    return len(reference_lines), len(read_lines)


def evaluate_selection(
    manifest: PhysicalDatasetManifest,
    selected: list[tuple[int, str, list[str]]],
) -> dict[str, Any]:
    """Score the selection with CNR's own cluster association."""
    summary = evaluate_supplied_clusters(
        manifest,
        [center for _, center, _ in selected],
        [reads for _, _, reads in selected],
    )
    return {**summary.to_dict(), "cluster_indices": [index for index, _, _ in selected]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("datasets/cnr.json"))
    parser.add_argument("--centers", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=20)
    parser.add_argument("--max-reads-per-cluster", type=int, default=10)
    parser.add_argument("--reads-out", type=Path, default=Path("cnr-reads.fasta"))
    parser.add_argument("--references-out", type=Path, default=Path("cnr-references.fasta"))
    parser.add_argument("--supplied-output", type=Path)
    args = parser.parse_args()

    hashes = verify_inputs(args.manifest, args.centers, args.clusters)
    centers = load_centers(args.centers)
    clusters = load_clusters(args.clusters, expected_count=len(centers))
    selected = select_clusters(
        centers,
        clusters,
        limit=args.limit,
        max_reads_per_cluster=args.max_reads_per_cluster,
    )
    references, reads = write_physical_inputs(
        selected,
        reads_out=args.reads_out,
        references_out=args.references_out,
    )
    print(f"verified inputs; wrote {references} references and {reads} reads")

    if args.supplied_output is not None:
        result = evaluate_selection(PhysicalDatasetManifest.load(args.manifest), selected)
        result["input_sha256"] = hashes
        args.supplied_output.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

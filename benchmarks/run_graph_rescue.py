"""Reproducible graph-rescue evidence for substitution and indel cases."""

from __future__ import annotations

import json
from pathlib import Path

from oligoark.archive import ArchiveConfig, archive_bytes
from oligoark.validation import compare_reconstruction_modes


def _substitute(sequence: str, position: int) -> str:
    current = sequence[position]
    replacement = next(base for base in "ACGT" if base != current)
    return sequence[:position] + replacement + sequence[position + 1 :]


def _insert(sequence: str, position: int, base: str) -> str:
    return sequence[:position] + base + sequence[position:]


def _delete(sequence: str, position: int) -> str:
    return sequence[:position] + sequence[position + 1 :]


def _archive(payload: bytes):
    return archive_bytes(
        payload,
        ArchiveConfig(
            chunk_size=96,
            rs_nsym=0,
            parity_group_size=8,
            redundancy_scheme="none",
            adaptive_masks=True,
            mask_search_limit=128,
        ),
    )


def substitution_case() -> dict[str, object]:
    archive = _archive(b"controlled substitution graph rescue evidence")
    original = archive.strands[0]
    positions = (72, 84, 96, 108, 120, 132, 144)
    reads = [_substitute(original, position) for position in positions]
    comparison = compare_reconstruction_modes(archive, reads, threshold=0.96)
    return {
        "case": "substitution-duplicates",
        "read_count": len(reads),
        **comparison.to_dict(),
    }


def indel_case() -> dict[str, object]:
    archive = _archive(b"controlled insertion deletion graph rescue evidence")
    original = archive.strands[0]
    insertions = (
        _insert(original, 86, "A"),
        _insert(original, 104, "C"),
        _insert(original, 122, "G"),
        _insert(original, 140, "T"),
    )
    deletions = (
        _delete(original, 94),
        _delete(original, 116),
        _delete(original, 134),
    )
    reads = list(insertions + deletions)
    comparison = compare_reconstruction_modes(archive, reads, threshold=0.94)
    return {
        "case": "insertion-deletion",
        "read_count": len(reads),
        **comparison.to_dict(),
    }


def main() -> None:
    results = [substitution_case(), indel_case()]
    output = Path("graph-rescue-results.json")
    output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()

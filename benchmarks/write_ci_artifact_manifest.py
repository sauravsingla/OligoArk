"""Create a SHA-256 inventory for software benchmark CI evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
from pathlib import Path


def manifest(paths: list[Path], root: Path) -> dict[str, object]:
    root = root.resolve()
    files: list[Path] = []
    for path in paths:
        resolved = root / path
        if resolved.is_dir():
            files.extend(p for p in resolved.rglob("*") if p.is_file())
        elif resolved.is_file():
            files.append(resolved)
        else:
            raise ValueError(f"Missing evidence input: {path}")

    entries = []
    for path in sorted(set(files)):
        if path.is_symlink() or not path.resolve().is_relative_to(root):
            raise ValueError(f"Evidence must be a regular workspace file: {path}")
        checksum = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                checksum.update(block)
        entries.append(
            {"path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size,
             "sha256": checksum.hexdigest()}
        )
    if not entries:
        raise ValueError("No evidence files found")
    return {
        "schema": "oligoark-ci-evidence-v1",
        "commit": os.getenv("GITHUB_SHA"),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "scope": "software benchmark only; not physical DNA validation",
        "files": entries,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("paths", nargs="+", type=Path)
    args = parser.parse_args()
    target = Path.cwd() / args.output
    if target.exists() or target.is_symlink():
        raise ValueError("Evidence manifest output already exists")
    if not target.resolve().is_relative_to(Path.cwd().resolve()):
        raise ValueError("Evidence manifest output is outside workspace")
    target.write_text(
        json.dumps(manifest(args.paths, Path.cwd()), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()

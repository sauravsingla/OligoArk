# Contributing

Thank you for improving OligoArk.

1. Create a focused branch and include tests for behavior changes.
2. Run `pytest`, `ruff check .`, `mypy src/oligoark`, and `python -m build`.
3. Keep deterministic baselines available even when adding ML methods.
4. Label claims as one of: **measured software result**, **simulation result**, **external published result**, or **hypothesis**.
5. Never present simulated synthesis/sequencing results as wet-lab evidence.
6. Cite external algorithms/papers and avoid copying code from repositories with incompatible licenses.

For research contributions, include a reproducible command, seed, dataset/source description, and machine-readable output.

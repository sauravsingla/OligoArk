# Contributing

Thank you for improving OligoArk.

1. Create a focused branch and include tests for behavior changes.
2. Run `pytest`, `ruff check .`, `mypy src/oligoark`, and `python -m build`.
3. Keep deterministic baselines available even when adding ML methods.
4. Label claims as one of: **measured software result**, **simulation result**, **external published result**, or **hypothesis**.
5. Never present simulated synthesis/sequencing results as wet-lab evidence.
6. Cite external algorithms/papers and avoid copying code from repositories with incompatible licenses.

For research contributions, include a reproducible command, explicit calibration/training and evaluation/test seed sets, dataset/source description, search method/seed, machine-readable raw output, and the exact commit. Calibration/training data must be disjoint from final evaluation/test data unless the contribution explicitly studies resubstitution bias.


For physical-data contributions, include the dataset DOI/accession, license or data-availability statement, exact reference-oligo source, preprocessing commands, read-assignment rule, and a clear distinction between reference reconstruction and end-to-end OligoArk archive recovery. Do not commit third-party sequencing data unless redistribution rights are explicit.

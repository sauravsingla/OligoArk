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

## Join a research challenge

OligoArk welcomes reproducible contributions from software engineers, students, and researchers:

- **[DNA insertion/deletion recovery challenge](https://github.com/sauravsingla/OligoArk/issues/16)** — improve exact recovery under indel noise while keeping timeouts bounded.
- **[Graph AI benchmark challenge](https://github.com/sauravsingla/OligoArk/issues/18)** — compare optional edge scorers against deterministic baselines on held-out data.
- **[Good first issue: deterministic graph baseline fixture](https://github.com/sauravsingla/OligoArk/issues/43)** — a smaller starting task; no GPU or external data required.

### Your first contribution

1. Run `python -m pip install -e ".[dev]"` then `python examples/quick_start.py` from a checkout.
2. Comment on the issue with your proposed approach, especially before starting substantial research work.
3. Create a focused branch and a small PR that links the issue.
4. Include exact reproduction commands, fixed seeds, expected vs actual results, and tests. Report **all** failures, timeouts, and regressions.
5. Use held-out test data for performance claims; never weaken integrity verification or describe computational simulations as physical experiments.

You can also start by reproducing an existing failure, improving documentation, or reviewing benchmark methodology. Useful negative results are welcome.

### Discuss ideas

Use [GitHub Discussions](https://github.com/sauravsingla/OligoArk/discussions) to share hypotheses, compare methods, and ask questions before opening a large PR. For a concrete bug or implementation task, use the linked issues instead.

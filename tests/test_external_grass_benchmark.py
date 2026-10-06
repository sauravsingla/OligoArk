from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "benchmarks"
    / "run_external_grass_benchmark.py"
)
SPEC = spec_from_file_location("run_external_grass_benchmark", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

load_grass_binned = MODULE.load_grass_binned
select_split = MODULE.select_split


def test_load_grass_binned_parses_reference_and_reads(tmp_path: Path) -> None:
    dataset = tmp_path / "Grass.txt"
    dataset.write_text(
        "ACGT\n******************\nACGT\nACGA\n\n\n"
        "TGCA\n******************\nTGCA\nTGTA\n\n",
        encoding="utf-8",
    )
    records = load_grass_binned(dataset)
    assert [record["reference"] for record in records] == ["ACGT", "TGCA"]
    assert records[0]["reads"] == ["ACGT", "ACGA"]
    assert records[1]["reads"] == ["TGCA", "TGTA"]


def test_select_split_is_deterministic_and_disjoint() -> None:
    records = [
        {
            "cluster_index": index,
            "reference": "ACGT",
            "reads": ["ACGT"] * 12,
        }
        for index in range(1, 21)
    ]
    first = select_split(
        records,
        size=5,
        seed=123,
        max_coverage=10,
    )
    second = select_split(
        records,
        size=5,
        seed=123,
        max_coverage=10,
    )
    assert [row["cluster_index"] for row in first] == [
        row["cluster_index"] for row in second
    ]

    excluded = {int(row["cluster_index"]) for row in first}
    held_out = select_split(
        records,
        size=5,
        seed=456,
        max_coverage=10,
        excluded_indices=excluded,
    )
    assert excluded.isdisjoint(
        {int(row["cluster_index"]) for row in held_out}
    )

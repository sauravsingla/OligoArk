from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "benchmarks"
    / "run_external_lcrc_ont_benchmark.py"
)
SPEC = spec_from_file_location("run_external_lcrc_ont_benchmark", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

FORWARD_PRIMER = MODULE.FORWARD_PRIMER
REVERSE_PRIMER = MODULE.REVERSE_PRIMER
build_unique_kmer_index = MODULE.build_unique_kmer_index
orient_and_map_read = MODULE.orient_and_map_read
reverse_complement = MODULE.reverse_complement
select_split = MODULE.select_split
trim_universal_primers = MODULE.trim_universal_primers


def _reference(payload: str) -> str:
    assert len(payload) == 160
    return FORWARD_PRIMER + payload + REVERSE_PRIMER


def test_lcrc_mapping_is_orientation_invariant() -> None:
    references = [
        _reference("ACGTTGCA" * 20),
        _reference("GATTACAG" * 20),
    ]
    index = build_unique_kmer_index(references)
    forward = orient_and_map_read(references[1], index)
    reverse = orient_and_map_read(reverse_complement(references[1]), index)
    assert forward is not None
    assert reverse is not None
    assert forward[0] == 1
    assert reverse[0] == 1
    assert forward[1] == references[1]
    assert reverse[1] == references[1]


def test_lcrc_primer_trim_uses_only_universal_anchors() -> None:
    reference = _reference("ACGTTGCA" * 20)
    read = "GGG" + reference + "TTT"
    assert trim_universal_primers(read) == reference


def test_lcrc_select_split_is_deterministic_and_disjoint() -> None:
    references = [_reference("ACGTTGCA" * 20) for _ in range(20)]
    clusters = [[f"READ{index}_{copy}" for copy in range(12)] for index in range(20)]

    first = select_split(
        references,
        clusters,
        size=5,
        seed=123,
        max_coverage=10,
    )
    repeated = select_split(
        references,
        clusters,
        size=5,
        seed=123,
        max_coverage=10,
    )
    assert [row["cluster_index"] for row in first] == [
        row["cluster_index"] for row in repeated
    ]

    excluded = {int(row["cluster_index"]) - 1 for row in first}
    held_out = select_split(
        references,
        clusters,
        size=5,
        seed=456,
        max_coverage=10,
        excluded_indices=excluded,
    )
    assert {
        int(row["cluster_index"]) for row in first
    }.isdisjoint({int(row["cluster_index"]) for row in held_out})

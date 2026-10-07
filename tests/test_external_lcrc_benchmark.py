from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "benchmarks"
    / "run_external_lcrc_benchmark.py"
)
SPEC = spec_from_file_location("run_external_lcrc_benchmark", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

assign_read = MODULE.assign_read
build_unique_anchor_index = MODULE.build_unique_anchor_index
load_fasta = MODULE.load_fasta
reverse_complement = MODULE.reverse_complement
select_split = MODULE.select_split


def test_load_fasta_parses_multiline_sequences(tmp_path: Path) -> None:
    path = tmp_path / "refs.fa"
    path.write_text(
        ">1\nACGTACGT\nACGT\n>2\nTGCATGCA\nTGCA\n",
        encoding="utf-8",
    )
    assert load_fasta(path) == ["ACGTACGTACGT", "TGCATGCATGCA"]


def test_anchor_mapper_assigns_forward_and_reverse_reads() -> None:
    primer_left = "ATAATTGGCTCCTGCTTGCA"
    primer_right = "AATGTAGGCGGAAAGTGCAA"
    interior_a = (
        "ACACGACTGCATTAGGGGATCCGCTGACAATTCTTTGACATCATCGCGGGGTAAACCAATGCCC"
        "AATGTATTATTGGTCGGCATAATGACCCCTCCGCGTTACGATTAGGAATATTCCAAACAATTTG"
        "CCTGTTAACTCAGACTGGAGATGGTGGGCCCCAA"
    )
    interior_b = (
        "CTCCAGTAGGCTTACGGTGTTAAATTACCGTTGACGAGCTTTGAGTTTGCATAGCTGCCCCGGC"
        "AGATCTAGTCTTCAGGCACGCGCTCCTAAGCAGAAATGGAGTTCCGTGAAGTTCCCTGTGGTCT"
        "CACTCGACATGAGCCCCCCGCAACCTGGGGGGAA"
    )
    references = [
        primer_left + interior_a[:160] + primer_right,
        primer_left + interior_b[:160] + primer_right,
    ]
    anchors = build_unique_anchor_index(references)

    assigned = assign_read(references[0], references, anchors)
    assert assigned is not None
    assert assigned[0] == 0
    assert assigned[1] == references[0]

    reverse = assign_read(reverse_complement(references[1]), references, anchors)
    assert reverse is not None
    assert reverse[0] == 1
    assert reverse[1] == references[1]


def test_select_split_is_deterministic_and_disjoint() -> None:
    references = ["A" * 200 for _ in range(20)]
    clusters = [["A" * 200] * 12 for _ in range(20)]
    first = select_split(
        references,
        clusters,
        size=5,
        seed=123,
        max_coverage=10,
    )
    second = select_split(
        references,
        clusters,
        size=5,
        seed=123,
        max_coverage=10,
    )
    assert [row["cluster_index"] for row in first] == [
        row["cluster_index"] for row in second
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
    assert excluded.isdisjoint(
        {int(row["cluster_index"]) - 1 for row in held_out}
    )

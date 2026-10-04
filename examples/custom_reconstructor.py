"""Demonstrate the reconstruction extension point without ML dependencies."""

from oligoark.reconstruct import GraphConsensusReconstructor


class PrefixAwareScorer:
    """Simple custom scorer used only to demonstrate the plug-in contract."""

    def score(self, left: str, right: str) -> float:
        prefix = 0
        for left_base, right_base in zip(left, right, strict=False):
            if left_base != right_base:
                break
            prefix += 1
        return prefix / max(1, len(left), len(right))


reads = ["ACGTACGT", "ACGTACGA", "TTTTGGGG"]
reconstructor = GraphConsensusReconstructor(
    threshold=0.5,
    scorer=PrefixAwareScorer(),
    use_qgram_prefilter=False,
)
result = reconstructor.reconstruct(reads)
print(result)

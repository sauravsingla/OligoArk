"""Physical strand-length profiles for OligoArk experiments."""

from __future__ import annotations

from dataclasses import dataclass, replace

from .archive import ArchiveConfig
from .framing import frame_overhead_bytes


@dataclass(frozen=True)
class PhysicalStrandProfile:
    name: str
    target_nucleotides: int
    rs_nsym: int = 8
    parity_group_size: int = 8
    redundancy_scheme: str = "xor"
    fountain_redundancy: float = 0.25
    fountain_max_degree: int = 4
    fountain_layout: str = "random"
    adaptive_masks: bool = True
    mask_search_limit: int = 64
    compact_framing: bool = False
    compact_index_bytes: int = 3
    inline_mask_framing: bool = False
    indel_rescue: bool = False

    @property
    def packed_bytes(self) -> int:
        return self.target_nucleotides // 4

    @property
    def chunk_size(self) -> int:
        overhead = frame_overhead_bytes(
            self.rs_nsym,
            compact_framing=self.compact_framing,
            compact_index_bytes=self.compact_index_bytes,
            inline_mask_framing=self.inline_mask_framing,
        )
        maximum_frame_payload = 255 - (overhead - 1)
        maximum = min(
            maximum_frame_payload,
            self.packed_bytes - overhead,
        )
        if maximum < 8:
            raise ValueError(
                f"target_nucleotides={self.target_nucleotides} is too short for "
                f"the current OligoArk frame with rs_nsym={self.rs_nsym}"
            )
        return maximum

    @property
    def actual_nucleotides(self) -> int:
        overhead = frame_overhead_bytes(
            self.rs_nsym,
            compact_framing=self.compact_framing,
            compact_index_bytes=self.compact_index_bytes,
            inline_mask_framing=self.inline_mask_framing,
        )
        return 4 * (overhead + self.chunk_size)

    def to_archive_config(self, **overrides: object) -> ArchiveConfig:
        values: dict[str, object] = {
            "chunk_size": self.chunk_size,
            "rs_nsym": self.rs_nsym,
            "parity_group_size": self.parity_group_size,
            "adaptive_masks": self.adaptive_masks,
            "redundancy_scheme": self.redundancy_scheme,
            "fountain_redundancy": self.fountain_redundancy,
            "fountain_max_degree": self.fountain_max_degree,
            "fountain_layout": self.fountain_layout,
            "mask_search_limit": self.mask_search_limit,
            "compact_framing": self.compact_framing,
            "compact_index_bytes": self.compact_index_bytes,
            "inline_mask_framing": self.inline_mask_framing,
            "indel_rescue": self.indel_rescue,
        }
        values.update(overrides)
        return ArchiveConfig.from_mapping(values)

    def with_scheme(self, scheme: str) -> PhysicalStrandProfile:
        return replace(self, redundancy_scheme=scheme)


PHYSICAL_STRAND_PROFILES: dict[str, PhysicalStrandProfile] = {
    # Frozen legacy profiles: these retain the framing and RS settings used by earlier
    # OligoArk releases and published repository benchmark artifacts.
    "oligoark-152": PhysicalStrandProfile("oligoark-152", 152),
    "oligoark-200": PhysicalStrandProfile("oligoark-200", 200),
    "oligoark-248": PhysicalStrandProfile("oligoark-248", 248),
    # Research profiles with compact framing and lower per-strand RS overhead. These are
    # explicit names so the density improvement never changes the meaning of legacy profiles.
    "oligoark-152-compact": PhysicalStrandProfile(
        "oligoark-152-compact",
        152,
        rs_nsym=2,
        compact_framing=True,
        indel_rescue=True,
    ),
    "oligoark-152-compact-v3": PhysicalStrandProfile(
        "oligoark-152-compact-v3",
        152,
        rs_nsym=2,
        parity_group_size=24,
        fountain_redundancy=5.0 / 24.0,
        fountain_max_degree=24,
        fountain_layout="interleaved",
        # Large archives encounter rare frames with no valid legacy 4-bit mask.
        # The extended selector preserves 152 nt and can search all 256 masks.
        mask_search_limit=256,
        compact_framing=True,
        inline_mask_framing=True,
        indel_rescue=True,
    ),
    "oligoark-200-compact": PhysicalStrandProfile(
        "oligoark-200-compact",
        200,
        rs_nsym=2,
        compact_framing=True,
        indel_rescue=True,
    ),
    "oligoark-248-compact": PhysicalStrandProfile(
        "oligoark-248-compact",
        248,
        rs_nsym=2,
        compact_framing=True,
        indel_rescue=True,
    ),
    "scale-1024": PhysicalStrandProfile(
        "scale-1024",
        1024,
        rs_nsym=0,
        adaptive_masks=False,
        mask_search_limit=1,
        compact_framing=False,
        indel_rescue=False,
    ),
}


def physical_strand_profile(name: str) -> PhysicalStrandProfile:
    try:
        return PHYSICAL_STRAND_PROFILES[name]
    except KeyError as exc:
        choices = ", ".join(sorted(PHYSICAL_STRAND_PROFILES))
        raise ValueError(
            f"Unknown physical strand profile {name!r}; choose from {choices}"
        ) from exc

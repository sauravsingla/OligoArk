"""Exact fountain droplet XOR parity against the original bytewise reference."""

import random

from oligoark.ecc import xor_bytes


def test_integer_droplet_xor_randomized_parity() -> None:
    rng = random.Random(20261101)
    for width in (1, 2, 7, 32, 64, 128):
        for _ in range(100):
            chunks = [rng.randbytes(rng.randrange(width + 1)) for _ in range(rng.randrange(1, 25))]
            value = 0
            for chunk in chunks:
                shift = 8 * (width - len(chunk))
                value ^= int.from_bytes(chunk, "big") << shift
            assert value.to_bytes(width, "big") == xor_bytes(chunks, width)

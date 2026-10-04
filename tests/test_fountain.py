from oligoark.fountain import make_symbols, peel_decode


def test_fountain_style_symbol_can_recover_missing_chunk() -> None:
    chunks = [bytes([i]) * 8 for i in range(5)]
    symbols = make_symbols(chunks, count=30, width=8, seed=10)
    known = {i: chunk for i, chunk in enumerate(chunks) if i != 2}
    recovered = peel_decode(known, symbols, total=5, width=8)
    assert recovered[2] == chunks[2]

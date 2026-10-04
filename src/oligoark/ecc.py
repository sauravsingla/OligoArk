"""Pure-Python Reed-Solomon and XOR erasure helpers.

The RS implementation uses GF(2^8) with primitive polynomial 0x11d. It is intentionally
small and self-contained so OligoArk's core codec is reproducible without a binary runtime.
"""

from __future__ import annotations

from dataclasses import dataclass

_PRIMITIVE = 0x11D
_GF_EXP = [0] * 512
_GF_LOG = [0] * 256
_x = 1
for _i in range(255):
    _GF_EXP[_i] = _x
    _GF_LOG[_x] = _i
    _x <<= 1
    if _x & 0x100:
        _x ^= _PRIMITIVE
for _i in range(255, 512):
    _GF_EXP[_i] = _GF_EXP[_i - 255]


class ECCDecodeError(ValueError):
    """Raised when Reed-Solomon recovery cannot correct the payload."""


def _gf_mul(x: int, y: int) -> int:
    if x == 0 or y == 0:
        return 0
    return _GF_EXP[_GF_LOG[x] + _GF_LOG[y]]


def _gf_div(x: int, y: int) -> int:
    if y == 0:
        raise ZeroDivisionError("GF division by zero")
    if x == 0:
        return 0
    return _GF_EXP[(_GF_LOG[x] + 255 - _GF_LOG[y]) % 255]


def _gf_pow(x: int, power: int) -> int:
    return _GF_EXP[(_GF_LOG[x] * power) % 255]


def _gf_inverse(x: int) -> int:
    return _GF_EXP[255 - _GF_LOG[x]]


def _poly_scale(poly: list[int], x: int) -> list[int]:
    return [_gf_mul(coef, x) for coef in poly]


def _poly_add(p: list[int], q: list[int]) -> list[int]:
    out = [0] * max(len(p), len(q))
    for i, value in enumerate(p):
        out[i + len(out) - len(p)] ^= value
    for i, value in enumerate(q):
        out[i + len(out) - len(q)] ^= value
    return out


def _poly_mul(p: list[int], q: list[int]) -> list[int]:
    out = [0] * (len(p) + len(q) - 1)
    for j, qj in enumerate(q):
        for i, pi in enumerate(p):
            out[i + j] ^= _gf_mul(pi, qj)
    return out


def _poly_eval(poly: list[int] | bytes | bytearray, x: int) -> int:
    y = poly[0]
    for coef in poly[1:]:
        y = _gf_mul(y, x) ^ coef
    return y


def _generator_poly(nsym: int) -> list[int]:
    gen = [1]
    for i in range(nsym):
        gen = _poly_mul(gen, [1, _gf_pow(2, i)])
    return gen


def rs_encode(data: bytes, nsym: int) -> bytes:
    """Append Reed-Solomon parity symbols to a message shorter than 255 symbols."""
    if nsym <= 0:
        return data
    if nsym >= 255 or len(data) + nsym > 255:
        raise ValueError("RS(255) frame requires len(data) + nsym <= 255")
    gen = _generator_poly(nsym)
    out = bytearray(data) + bytearray(nsym)
    for i in range(len(data)):
        coef = out[i]
        if coef:
            for j in range(1, len(gen)):
                out[i + j] ^= _gf_mul(gen[j], coef)
    out[: len(data)] = data
    return bytes(out)


def _syndromes(msg: bytes | bytearray, nsym: int) -> list[int]:
    return [0] + [_poly_eval(msg, _gf_pow(2, i)) for i in range(nsym)]


def _find_error_locator(synd: list[int], nsym: int) -> list[int]:
    err_loc = [1]
    old_loc = [1]
    for i in range(nsym):
        k = i + 1
        delta = synd[k]
        for j in range(1, len(err_loc)):
            delta ^= _gf_mul(err_loc[-(j + 1)], synd[k - j])
        old_loc.append(0)
        if delta:
            if len(old_loc) > len(err_loc):
                new_loc = _poly_scale(old_loc, delta)
                old_loc = _poly_scale(err_loc, _gf_inverse(delta))
                err_loc = new_loc
            err_loc = _poly_add(err_loc, _poly_scale(old_loc, delta))
    while len(err_loc) and err_loc[0] == 0:
        del err_loc[0]
    errs = len(err_loc) - 1
    if errs * 2 > nsym:
        raise ECCDecodeError("too many symbol errors for configured Reed-Solomon parity")
    return err_loc


def _find_errors(err_loc: list[int], message_len: int) -> list[int]:
    errs = len(err_loc) - 1
    positions: list[int] = []
    for i in range(message_len):
        if _poly_eval(err_loc, _gf_pow(2, i)) == 0:
            positions.append(message_len - 1 - i)
    if len(positions) != errs:
        raise ECCDecodeError("could not locate all Reed-Solomon symbol errors")
    return positions


def _find_errata_locator(coef_positions: list[int]) -> list[int]:
    loc = [1]
    for position in coef_positions:
        loc = _poly_mul(loc, [_gf_pow(2, position), 1])
    return loc


def _find_error_evaluator(synd: list[int], err_loc: list[int], nsym: int) -> list[int]:
    product = _poly_mul(synd, err_loc)
    return product[-(nsym + 1) :]


def _correct_errata(msg: bytearray, synd: list[int], err_positions: list[int]) -> bytearray:
    coef_positions = [len(msg) - 1 - p for p in err_positions]
    err_loc = _find_errata_locator(coef_positions)
    err_eval = _find_error_evaluator(list(reversed(synd)), err_loc, len(err_loc) - 1)
    err_eval = list(reversed(err_eval))
    x_values = [_gf_pow(2, -(255 - p)) for p in coef_positions]
    corrections = bytearray(len(msg))
    for i, xi in enumerate(x_values):
        xi_inv = _gf_inverse(xi)
        locator_prime = 1
        for j, xj in enumerate(x_values):
            if j != i:
                locator_prime = _gf_mul(locator_prime, 1 ^ _gf_mul(xi_inv, xj))
        y = _poly_eval(list(reversed(err_eval)), xi_inv)
        y = _gf_mul(_gf_pow(xi, 1), y)
        magnitude = _gf_div(y, locator_prime)
        corrections[err_positions[i]] = magnitude
    return bytearray(_poly_add(list(msg), list(corrections)))


def rs_decode(data: bytes, nsym: int) -> bytes:
    """Correct errors and return the original message bytes."""
    if nsym <= 0:
        return data
    if len(data) > 255 or len(data) <= nsym:
        raise ECCDecodeError("invalid Reed-Solomon codeword length")
    msg = bytearray(data)
    synd = _syndromes(msg, nsym)
    if max(synd) == 0:
        return bytes(msg[:-nsym])
    err_loc = _find_error_locator(synd, nsym)
    err_positions = _find_errors(list(reversed(err_loc)), len(msg))
    corrected = _correct_errata(msg, synd, err_positions)
    if max(_syndromes(corrected, nsym)) != 0:
        raise ECCDecodeError("Reed-Solomon correction did not converge")
    return bytes(corrected[:-nsym])


@dataclass(frozen=True)
class ParityBlock:
    group_index: int
    payload: bytes


def xor_bytes(items: list[bytes], width: int) -> bytes:
    acc = bytearray(width)
    for item in items:
        padded = item.ljust(width, b"\x00")
        for i, value in enumerate(padded[:width]):
            acc[i] ^= value
    return bytes(acc)


def build_xor_parity(chunks: list[bytes], group_size: int, width: int) -> list[ParityBlock]:
    if group_size <= 1:
        return []
    parity: list[ParityBlock] = []
    for group_index, start in enumerate(range(0, len(chunks), group_size)):
        payload = xor_bytes(chunks[start : start + group_size], width)
        parity.append(ParityBlock(group_index, payload))
    return parity


def recover_one_missing(
    known: dict[int, bytes],
    parity: dict[int, bytes],
    total_chunks: int,
    group_size: int,
    width: int,
) -> dict[int, bytes]:
    """Recover at most one missing data chunk per XOR parity group."""
    recovered = dict(known)
    if group_size <= 1:
        return recovered
    for group_index, start in enumerate(range(0, total_chunks, group_size)):
        end = min(start + group_size, total_chunks)
        missing = [idx for idx in range(start, end) if idx not in recovered]
        if len(missing) != 1 or group_index not in parity:
            continue
        pieces = [recovered[idx] for idx in range(start, end) if idx in recovered]
        pieces.append(parity[group_index])
        recovered[missing[0]] = xor_bytes(pieces, width)
    return recovered

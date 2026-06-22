"""HAM: Integer encoding utilities for DNA sequences.

All DNA sequences are encoded as 2-bit-per-base integers:
  A=00, C=01, G=10, T=11

This module provides the shared encoding/decoding functions used across
matcher, dedup, and hash_builder modules.
"""

# ── Base-to-bits mapping ──
_BASE2BITS = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
_BITS2BASE = 'ACGT'
BASES = ('A', 'C', 'G', 'T')

# ── Byte-level lookup table (256 entries) ──
_BYTE2BITS = [-1] * 256
_BYTE2BITS[ord('A')] = 0
_BYTE2BITS[ord('C')] = 1
_BYTE2BITS[ord('G')] = 2
_BYTE2BITS[ord('T')] = 3

# ── 5bp chunk LUT (P0 optimisation) ──
# Index: b0<<8 | b1<<6 | b2<<4 | b3<<2 | b4.  Identity function since
# the index IS the correctly-ordered 10-bit encoding.
_5MER_ENC = list(range(1024))


def encode_barcode(bc: str) -> int:
    """Encode a 16bp DNA barcode into a 32-bit integer (2 bits/base)."""
    val = 0
    for ch in bc:
        val = (val << 2) | _BASE2BITS[ch]
    return val


def encode_seq(seq: str) -> int:
    """Encode a DNA sequence string as uint64 (2 bits/base)."""
    val = 0
    for ch in seq:
        val = (val << 2) | _BASE2BITS[ch]
    return val


def encode_seq_bytes(b: bytes, start: int, length: int) -> int:
    """Encode a DNA subsequence from bytes using pre-computed LUT."""
    val = 0
    for i in range(start, start + length):
        val = (val << 2) | _BYTE2BITS[b[i]]
    return val


def encode_umi_bytes(b: bytes, start: int) -> int:
    """Encode 12bp UMI from bytes as uint32 (24 bits), loop-unrolled."""
    v = _BYTE2BITS
    return (v[b[start]]   << 22 | v[b[start+1]]  << 20 |
            v[b[start+2]] << 18 | v[b[start+3]]  << 16 |
            v[b[start+4]] << 14 | v[b[start+5]]  << 12 |
            v[b[start+6]] << 10 | v[b[start+7]]  << 8  |
            v[b[start+8]] << 6  | v[b[start+9]]  << 4  |
            v[b[start+10]] << 2 | v[b[start+11]])


def decode_umi(val: int) -> str:
    """Decode uint32 UMI back to 12bp string."""
    chars = []
    for _ in range(12):
        chars.append(_BITS2BASE[val & 0x3])
        val >>= 2
    return ''.join(reversed(chars))


def encode_window_bigint(window_bytes: bytes) -> int:
    """Encode 26bp window as a 52-bit integer (P0 optimised).

    Uses 5x 5bp chunk LUT lookups + 1 trailing base, loop-unrolled.
    """
    v = _BYTE2BITS
    b = window_bytes
    i0 = (v[b[0]]  << 8 | v[b[1]]  << 6 | v[b[2]]  << 4 |
          v[b[3]]  << 2 | v[b[4]])
    i1 = (v[b[5]]  << 8 | v[b[6]]  << 6 | v[b[7]]  << 4 |
          v[b[8]]  << 2 | v[b[9]])
    i2 = (v[b[10]] << 8 | v[b[11]] << 6 | v[b[12]] << 4 |
          v[b[13]] << 2 | v[b[14]])
    i3 = (v[b[15]] << 8 | v[b[16]] << 6 | v[b[17]] << 4 |
          v[b[18]] << 2 | v[b[19]])
    i4 = (v[b[20]] << 8 | v[b[21]] << 6 | v[b[22]] << 4 |
          v[b[23]] << 2 | v[b[24]])
    last = v[b[25]]
    return (_5MER_ENC[i0] << 42 |
            _5MER_ENC[i1] << 32 |
            _5MER_ENC[i2] << 22 |
            _5MER_ENC[i3] << 12 |
            _5MER_ENC[i4] << 2  |
            last)


def extract_guides_from_bigint(big_int: int) -> list:
    """Extract 7x 20bp guide codes from a 52-bit window int via shift+mask."""
    MASK_40 = (1 << 40) - 1
    guides = []
    for offset in range(7):
        shift = (26 - offset - 20) * 2
        guides.append((big_int >> shift) & MASK_40)
    return guides


def generate_cb_variants(bc: str, max_hamming: int = 2):
    """Generate all Hamming <= max_hamming variants of a 16bp barcode.

    Yields (variant_str, hamming_distance).
    """
    yield bc, 0
    n = len(bc)
    # Hamming = 1
    for p1 in range(n):
        o1 = bc[p1]
        for a1 in BASES:
            if a1 == o1:
                continue
            yield bc[:p1] + a1 + bc[p1+1:], 1
    # Hamming = 2
    if max_hamming >= 2:
        for p1 in range(n):
            o1 = bc[p1]
            for a1 in BASES:
                if a1 == o1:
                    continue
                v1 = bc[:p1] + a1 + bc[p1+1:]
                for p2 in range(p1+1, n):
                    o2 = bc[p2]
                    for a2 in BASES:
                        if a2 == o2:
                            continue
                        yield v1[:p2] + a2 + v1[p2+1:], 2

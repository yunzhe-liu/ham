"""Regression tests for ham.encoding — window/guide bit-packing.

The extract_guides_from_bigint shift formula regressed in a past refactor
(commit 00913be) for guide_len==20 (the default 10xv3 chemistry): every
sliding offset extracted a truncated, zero-padded fragment instead of a
real guide_len-base substring, so exact/Hamming=1 matching against a real
guide hash silently produced ~0 hits. These tests pin down the correct,
generic behaviour so any future refactor of the window/guide encoding path
is caught immediately instead of shipping silently.
"""

from ham.encoding import (
    encode_seq, encode_window_bigint, extract_guides_from_bigint,
)


def _guide_at_offset(guide: str, offset: int, window_len: int, filler: str = "A") -> bytes:
    """Build a window of window_len bases with `guide` placed at `offset`."""
    pad_left = filler * offset
    pad_right = filler * (window_len - offset - len(guide))
    window = pad_left + guide + pad_right
    assert len(window) == window_len
    return window.encode()


def test_encode_window_bigint_matches_encode_seq_for_full_window():
    # When the window IS the guide (no slack), the window encoding must be
    # bit-identical to encoding the guide sequence directly.
    seq = "ACGTACGTACGTACGTACG"  # 19bp
    assert encode_window_bigint(seq.encode()) == encode_seq(seq)


def test_extract_guides_10xv3_shape_26bp_window_20bp_guide():
    """Regression test for the 00913be shift-formula bug."""
    guide = "ACGTACGTACGTACGTACGT"  # 20bp
    window_len, guide_len = 26, 20
    for offset in range(window_len - guide_len + 1):  # 0..6
        window_bytes = _guide_at_offset(guide, offset, window_len)
        big_int = encode_window_bigint(window_bytes)
        candidates = extract_guides_from_bigint(big_int, guide_len, window_len)
        assert len(candidates) == 7
        assert encode_seq(guide) in candidates, (
            f"guide at offset={offset} was not recovered by any sliding candidate"
        )


def test_extract_guides_5p_shape_19bp_window_19bp_guide_no_margin():
    guide = "ACGTACGTACGTACGTACG"[:19]
    window_len, guide_len = 19, 19
    big_int = encode_window_bigint(guide.encode())
    candidates = extract_guides_from_bigint(big_int, guide_len, window_len)
    assert candidates == [encode_seq(guide)]


def test_extract_guides_5p_shape_with_margin_recovers_drifted_guide():
    """New 25bp window (+/-3bp margin) around the 19bp 5' guide must recover
    the guide regardless of where within the margin it actually sits."""
    guide = "ACGTACGTACGTACGTACG"[:19]
    window_len, guide_len = 25, 19
    for offset in range(window_len - guide_len + 1):  # 0..6
        window_bytes = _guide_at_offset(guide, offset, window_len)
        big_int = encode_window_bigint(window_bytes)
        candidates = extract_guides_from_bigint(big_int, guide_len, window_len)
        assert encode_seq(guide) in candidates


def test_extract_guides_out_of_margin_does_not_falsely_match():
    """A guide sitting entirely outside the window must never appear among
    the extracted candidates (sanity check against false positives)."""
    guide = "ACGTACGTACGTACGTACG"[:19]
    window_len, guide_len = 19, 19  # no margin: only exact placement matches
    other = "TTTTTTTTTTTTTTTTTTT"
    big_int = encode_window_bigint(other.encode())
    candidates = extract_guides_from_bigint(big_int, guide_len, window_len)
    assert encode_seq(guide) not in candidates

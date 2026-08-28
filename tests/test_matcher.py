"""End-to-end regression tests: synthetic FASTQ pairs through ham match.

These pin down the fix for the guide_len==20 shift-formula regression
(commit 00913be broke 10xv3 exact/Hamming=1 matching entirely), and the 5'
chemistry window: guide anchor kept at the validated R2 position 16,
guide_len corrected to 20bp (standard SpCas9 protospacer length, was
19bp), with +/-3bp margin added on each side (mirroring 10xv3's
tolerance) -> window R2[13:39].
"""

import gzip

from ham.hash_builder import build_guide_hash
from ham.matcher import match_reads


def _write_fastq(path, seq):
    with open(path, "w") as f:
        f.write(f"@read1\n{seq}\n+\n{'I' * len(seq)}\n")


def _build_guide_hash(tmp_path, guide_seq, name="guides.pkl"):
    fasta = tmp_path / "guides.fasta"
    fasta.write_text(f">guide1\n{guide_seq}\n")
    return build_guide_hash(str(fasta), str(tmp_path / name))


CB = "AAAAAAAAAAAAAAAA"  # 16bp, arbitrary


def test_10xv3_exact_match_at_documented_position(tmp_path):
    """R2[28:48] holding the guide must match under 10xv3 (regression test
    for the commit-00913be shift-formula bug)."""
    guide = "ACGTACGTACGTACGTACGT"  # 20bp
    gh = _build_guide_hash(tmp_path, guide)

    r1 = str(tmp_path / "R1.fastq")
    r2 = str(tmp_path / "R2.fastq")
    _write_fastq(r1, CB + "T" * 12)  # 16bp CB + 12bp UMI
    _write_fastq(r2, "C" * 28 + guide + "G" * 6)  # window[28:54] = guide + 6bp

    wl = {CB}
    stats = match_reads(r1, r2, wl, gh, str(tmp_path / "hits.npz"), chemistry="10xv3")

    assert stats["matched"] == 1
    assert stats["exact_hit"] == 1


def test_5p_chemistry_matches_guide_at_validated_anchor(tmp_path):
    """10xv2-5p window is R2[13:39] (26bp), guide_len=20, guide anchored at
    the validated R2 position 16 (offset=3 within the window, i.e. the
    center candidate)."""
    guide = "ACGTACGTACGTACGTACGT"  # 20bp
    gh = _build_guide_hash(tmp_path, guide)

    r1 = str(tmp_path / "R1.fastq")
    r2 = str(tmp_path / "R2.fastq")
    _write_fastq(r1, CB + "T" * 10)  # 16bp CB + 10bp UMI (10xv2-5p)
    _write_fastq(r2, "C" * 16 + guide + "G" * 3)  # guide at anchor pos 16
    wl = {CB}

    stats = match_reads(r1, r2, wl, gh, str(tmp_path / "hits.npz"), chemistry="10xv2-5p")
    assert stats["matched"] == 1
    assert stats["exact_hit"] == 1


def test_5p_chemistry_recovers_guide_within_margin(tmp_path):
    """+/-3bp drift off the validated anchor (R2 position 16) must still
    match."""
    guide = "ACGTACGTACGTACGTACGT"
    gh = _build_guide_hash(tmp_path, guide)

    r1 = str(tmp_path / "R1.fastq")
    _write_fastq(r1, CB + "T" * 10)
    wl = {CB}

    for drift in (-3, 0, 3):
        prefix_len = 16 + drift
        r2 = str(tmp_path / f"R2_{drift}.fastq")
        r2_seq = "C" * prefix_len + guide + "G" * (45 - prefix_len - len(guide))
        _write_fastq(r2, r2_seq[:45])

        stats = match_reads(
            r1, r2, wl, gh, str(tmp_path / f"hits_{drift}.npz"),
            chemistry="10xv2-5p",
        )
        assert stats["matched"] == 1, f"expected a match at drift={drift}"


def test_5p_chemistry_does_not_match_beyond_margin(tmp_path):
    """A guide drifted 4bp beyond the anchor falls outside the +/-3bp
    margin and must NOT match."""
    guide = "ACGTACGTACGTACGTACGT"
    gh = _build_guide_hash(tmp_path, guide)

    r1 = str(tmp_path / "R1.fastq")
    _write_fastq(r1, CB + "T" * 10)
    wl = {CB}

    prefix_len = 16 + 4  # 4bp beyond the +/-3bp margin
    r2 = str(tmp_path / "R2_out.fastq")
    r2_seq = "C" * prefix_len + guide + "G" * (45 - prefix_len - len(guide))
    _write_fastq(r2, r2_seq[:45])

    stats = match_reads(
        r1, r2, wl, gh, str(tmp_path / "hits_out.npz"), chemistry="10xv2-5p",
    )
    assert stats["matched"] == 0

"""Named-chemistry field override: chem_cfg may be passed alongside a named
chemistry (not just 'custom') to override individual fields on top of that
chemistry's preset, e.g. a non-standard guide_len with an otherwise-standard
10xv2-5p layout."""

import pytest

from ham.hash_builder import build_guide_hash
from ham.matcher import CHEMISTRY_CONFIGS, match_reads


CB = "AAAAAAAAAAAAAAAA"


def _guide_hash(tmp_path, guide):
    fasta = tmp_path / "g.fasta"
    fasta.write_text(f">g1\n{guide}\n")
    return build_guide_hash(str(fasta), str(tmp_path / "g.pkl"))


def _write_fastq(path, seq):
    path.write_text(f"@r\n{seq}\n+\n{'I' * len(seq)}\n")


def test_named_chemistry_guide_len_override(tmp_path):
    """10xv2-5p defaults to guide_len=20; overriding to 19 (this dataset's
    actual protospacer length) must be honoured while cb/umi/window
    positions keep the named chemistry's defaults."""
    guide = "ACGTACGTACGTACGTACG"  # 19bp
    gh = _guide_hash(tmp_path, guide)

    r1 = tmp_path / "R1.fastq"
    r2 = tmp_path / "R2.fastq"
    _write_fastq(r1, CB + "T" * 10)  # unchanged 10xv2-5p UMI layout
    _write_fastq(r2, "C" * 16 + guide + "G" * 10)  # unchanged anchor pos 16
    wl = {CB}

    stats = match_reads(
        str(r1), str(r2), wl, gh, str(tmp_path / "hits.npz"),
        chemistry="10xv2-5p", chem_cfg={"guide_len": 19},
    )
    assert stats["matched"] == 1
    assert stats["exact_hit"] == 1


def test_named_chemistry_without_override_is_unaffected(tmp_path):
    """Passing no chem_cfg still uses the named chemistry's own defaults
    unchanged (guide_len=20 for 10xv2-5p)."""
    assert CHEMISTRY_CONFIGS["10xv2-5p"]["guide_len"] == 20
    guide = "ACGTACGTACGTACGTACGT"  # 20bp
    gh = _guide_hash(tmp_path, guide)

    r1 = tmp_path / "R1.fastq"
    r2 = tmp_path / "R2.fastq"
    _write_fastq(r1, CB + "T" * 10)
    _write_fastq(r2, "C" * 16 + guide + "G" * 10)
    wl = {CB}

    stats = match_reads(
        str(r1), str(r2), wl, gh, str(tmp_path / "hits.npz"),
        chemistry="10xv2-5p",
    )
    assert stats["matched"] == 1


def test_named_chemistry_rejects_unknown_override_key(tmp_path):
    gh = _guide_hash(tmp_path, "ACGTACGTACGTACGTACGT")
    r1 = tmp_path / "R1.fastq"
    r2 = tmp_path / "R2.fastq"
    _write_fastq(r1, CB + "T" * 10)
    _write_fastq(r2, "C" * 45)
    wl = {CB}

    with pytest.raises(ValueError, match="unknown keys"):
        match_reads(
            str(r1), str(r2), wl, gh, str(tmp_path / "hits.npz"),
            chemistry="10xv2-5p", chem_cfg={"not_a_real_field": 1},
        )


def test_custom_chemistry_still_requires_full_chem_cfg(tmp_path):
    gh = _guide_hash(tmp_path, "ACGTACGTACGTACGTACGT")
    r1 = tmp_path / "R1.fastq"
    r2 = tmp_path / "R2.fastq"
    _write_fastq(r1, CB + "T" * 10)
    _write_fastq(r2, "C" * 45)
    wl = {CB}

    with pytest.raises(ValueError, match="missing keys"):
        match_reads(
            str(r1), str(r2), wl, gh, str(tmp_path / "hits.npz"),
            chemistry="custom", chem_cfg={"guide_len": 19},  # only 1 of 7 keys
        )

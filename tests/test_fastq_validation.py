import pytest

from ham.hash_builder import build_guide_hash
from ham.matcher import match_reads


CB = "AAAAAAAAAAAAAAAA"


def _guide_hash(tmp_path):
    fasta = tmp_path / "g.fasta"
    fasta.write_text(">g1\nACGTACGTACGTACGTACGT\n")
    return build_guide_hash(str(fasta), str(tmp_path / "g.pkl"))


def test_missing_fastq_file_raises(tmp_path):
    gh = _guide_hash(tmp_path)
    with pytest.raises(FileNotFoundError):
        match_reads(
            str(tmp_path / "nope_R1.fastq"), str(tmp_path / "nope_R2.fastq"),
            {CB}, gh, str(tmp_path / "hits.npz"), chemistry="10xv3")


def test_r1_r2_record_count_mismatch_raises(tmp_path):
    gh = _guide_hash(tmp_path)
    r1 = tmp_path / "R1.fastq"
    r2 = tmp_path / "R2.fastq"
    # R1 has 2 records, R2 has only 1.
    r1.write_text("@r1\n" + "A" * 28 + "\n+\n" + "I" * 28 + "\n"
                   "@r2\n" + "A" * 28 + "\n+\n" + "I" * 28 + "\n")
    r2.write_text("@r1\n" + "C" * 54 + "\n+\n" + "I" * 54 + "\n")

    with pytest.raises(ValueError, match="record count mismatch"):
        match_reads(str(r1), str(r2), {CB}, gh, str(tmp_path / "hits.npz"),
                    chemistry="10xv3")


def test_malformed_fastq_header_raises(tmp_path):
    gh = _guide_hash(tmp_path)
    r1 = tmp_path / "R1.fastq"
    r2 = tmp_path / "R2.fastq"
    r1.write_text("NOT_A_HEADER\n" + "A" * 28 + "\n+\n" + "I" * 28 + "\n")
    r2.write_text("@r1\n" + "C" * 54 + "\n+\n" + "I" * 54 + "\n")

    with pytest.raises(ValueError, match="Malformed FASTQ"):
        match_reads(str(r1), str(r2), {CB}, gh, str(tmp_path / "hits.npz"),
                    chemistry="10xv3")

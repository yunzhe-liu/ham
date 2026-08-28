from ham.encoding import encode_seq
from ham.hash_builder import build_guide_hash


def test_build_guide_hash_exact_and_hamming1(tmp_path):
    fasta = tmp_path / "guides.fasta"
    fasta.write_text(">g1\nACGTACGTACGTACGTACGT\n>g2\nTTTTAAAACCCCGGGGTTTT\n")

    data = build_guide_hash(str(fasta), str(tmp_path / "guide_hash.pkl"))

    assert data["n_guides"] == 2
    assert data["idx_to_id"] == ["g1", "g2"]

    g1_idx = data["seq_to_idx"]["ACGTACGTACGTACGTACGT"]
    assert data["seq_to_idx_int"][encode_seq("ACGTACGTACGTACGTACGT")] == g1_idx

    # Single mismatch at position 0 must resolve to the same parent guide.
    variant = "CCGTACGTACGTACGTACGT"
    assert data["seq_to_idx"][variant] == g1_idx
    assert data["seq_to_idx_int"][encode_seq(variant)] == g1_idx

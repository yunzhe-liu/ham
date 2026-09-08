from ham.encoding import encode_barcode, generate_cb_variants
from ham.matcher import build_cb_hash, match_reads


WHITELIST_BARCODE = "AAAAAAAAAAAAAAAA"
HAMMING1_BARCODE = "CAAAAAAAAAAAAAAA"
HAMMING2_BARCODE = "CCAAAAAAAAAAAAAA"


def test_cell_barcode_defaults_to_hamming1():
    cb_hash = build_cb_hash({WHITELIST_BARCODE})

    assert encode_barcode(HAMMING1_BARCODE) in cb_hash["cb_map"]
    assert encode_barcode(HAMMING2_BARCODE) not in cb_hash["cb_map"]
    assert len(list(generate_cb_variants(WHITELIST_BARCODE))) == 49
    assert match_reads.__defaults__[3] == 1


def test_cell_barcode_hamming2_remains_an_explicit_option():
    cb_hash = build_cb_hash({WHITELIST_BARCODE}, cb_max_hamming=2)

    assert encode_barcode(HAMMING2_BARCODE) in cb_hash["cb_map"]

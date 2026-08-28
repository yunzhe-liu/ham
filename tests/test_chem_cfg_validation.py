import pytest

from ham.matcher import CHEMISTRY_CONFIGS, _validate_chem_cfg


@pytest.mark.parametrize("chemistry", ["10xv3", "10xv2-5p", "10xv2-5p-12umi"])
def test_named_chemistries_are_valid(chemistry):
    _validate_chem_cfg(CHEMISTRY_CONFIGS[chemistry], source=chemistry)


def test_rejects_missing_keys():
    with pytest.raises(ValueError, match="missing keys"):
        _validate_chem_cfg({"cb_start": 0, "cb_end": 16}, source="custom")


def test_rejects_window_shorter_than_guide():
    cfg = dict(CHEMISTRY_CONFIGS["10xv3"])
    cfg["guide_len"] = 100  # window is only 26bp
    with pytest.raises(ValueError, match="shorter than guide_len"):
        _validate_chem_cfg(cfg, source="custom")


def test_rejects_inverted_positions():
    cfg = dict(CHEMISTRY_CONFIGS["10xv3"])
    cfg["cb_end"] = cfg["cb_start"]
    with pytest.raises(ValueError, match="cb_end"):
        _validate_chem_cfg(cfg, source="custom")

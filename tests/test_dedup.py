from ham.dedup import dedup_umis_directional


def test_dedup_collapses_hamming1_neighbours_toward_dominant_umi():
    # 10 copies of the dominant UMI, 1 copy of a Hamming=1 neighbour that
    # should be absorbed into it.
    umis = ["AAAAAAAAAAAA"] * 10 + ["AAAAAAAAAAAC"]
    assert dedup_umis_directional(umis, threshold=1) == 1


def test_dedup_keeps_distant_umis_separate():
    umis = ["AAAAAAAAAAAA"] * 5 + ["TTTTTTTTTTTT"] * 5
    assert dedup_umis_directional(umis, threshold=1) == 2


def test_dedup_threshold_zero_is_exact_dedup_only():
    umis = ["AAAAAAAAAAAA", "AAAAAAAAAAAC", "AAAAAAAAAAAA"]
    assert dedup_umis_directional(umis, threshold=0) == 2


def test_dedup_empty_input():
    assert dedup_umis_directional([], threshold=1) == 0

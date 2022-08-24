#!/usr/bin/env ipython
from pathlib import Path
import bctools as bc
import pandas as pd
import os.path

HERE = Path(__file__).parent


def test_load():
    file1 = str(HERE / "data/SLX-22398.SITTA1.s_1.r_1_small.fq")
    file2 = str(HERE / "data/SLX-22398.SITTA1.s_1.r_2_small.fq")
    files = {"r1": file1, "r2": file2}
    bcu = bc.get_barcodes(files)


def test_debug():
    file1 = str(HERE / "data/SLX-22398.SITTA1.s_1.r_1_small.fq")
    file2 = str(HERE / "data/SLX-22398.SITTA1.s_1.r_2_small.fq")
    files = {"r1": file1, "r2": file2}
    bcu, tempDIR = bc.larry.process_larry(files, debug=True)

    assert isinstance(bcu, bc.cbu.CBUSeries)
    assert os.path.isfile(tempDIR.name + "/cbcumi1.fq")
    assert os.path.isfile(tempDIR.name + "/bar1.fq")
    assert os.path.isfile(tempDIR.name + "/cbcumi2.fq")
    assert os.path.isfile(tempDIR.name + "/bar2.fq")


def test_cutadapt():
    file1 = str(HERE / "data/mock_larry_cutadapt_r1.fq")
    file2 = str(HERE / "data/mock_larry_cutadapt_r2.fq")
    files = {"r1": file1, "r2": file2}
    cbu = bc.get_barcodes(files).sort_index()

    # fmt: off
    ind = pd.MultiIndex.from_tuples([
        ('AACCTGAAGTATTCCG', 'CACGTGCCGCCAAATCACTGCAGAAGCCGTTGGCAGCCAA', 'CGATCAGATGAA'),
        ('AACCTGAAGTATTCCG', 'CACGTGCCGCCAAATCACTGCAGAAGCCGTTGGCAGCCAA', 'CGATCAGATGAT'),
        ('AACCTGAAGTATTCCG', 'CACGTGCCGCCAAATCACTGCAGAAGCCGTTGGCAGCCAA', 'CGATCAGATGGT'),
        ('AACCTGAAGTATTCCG', 'CACGTGCCGCCAAATCACTGCAGAAGCCGTTGGCAGCCAA', 'CGATCAGATGCT'),
        ('TTTCTGAAGTATTCCG', 'GGGGTGCCGCCATTGGACTCGAGAAGCCGTTGGCAGGGTT', 'CGATCAGATGGT'),
        ('GGGGTGAAGTATTCCG', 'CCCCTGGGGGCATTGGACTCGAGAAGCCGTTGGCAGAAAA', 'CGATCAGATGGT')
    ])
    # fmt: on
    ref = pd.Series([1, 1, 1, 1, 1, 1], index=ind)
    ref = ref.sort_index()

    assert cbu.equals(ref)


def test_count():
    file1 = str(HERE / "data/mock_larry_counting_r1.fq")
    file2 = str(HERE / "data/mock_larry_counting_r2.fq")
    files = {"r1": file1, "r2": file2}
    cbu = bc.get_barcodes(files).sort_index()

    # fmt: off
    ind = pd.MultiIndex.from_tuples([
        ('AACCTGAAGTATTCCG', 'CACGTGCCGCCAAATCACTGCAGAAGCCGTTGGCAGCCAA', 'CGATCAGATAAA'),
        ('AACCTGAAGTATTCCG', 'CACGTGCCGCCAAATCACTGCAGAAGCCGTTGGCAGCCAA', 'CGATCAGATGGG'),
        ('AACCTGAAGTATTCCG', 'GGGGTGCCGCCATTGGACTCGAGAAGCCGTTGGCAGGGTT', 'CGATCAGATCCC'),
        ('TTGGAGAAGTATTCCG', 'GGGGTGCCGCCATTGGACTCGAGAAGCCGTTGGCAGGGTT', 'CGATCAGATCCC'),
        ('AAAAAGAAGTATTCCG', 'CCCCTGGGGGCATTGGACTCGAGAAGCCGTTGGCAGAAAA', 'CGATCAGATAAA'),
        ('AAAAAGAAGTATTCCG', 'CCCCTGGGGGCATTGGACTCGAGAAGCCGTTGGCAGAAAA', 'CGATCAGATTTT')
        # ('', '', ''),
    ])
    # fmt: on
    ref = pd.Series([2, 1, 1, 4, 5, 3], index=ind)
    ref = ref.sort_index()

    assert cbu.equals(ref)


def test_count_CBC():
    file1 = str(HERE / "data/mock_larry_counting_r1.fq")
    file2 = str(HERE / "data/mock_larry_counting_r2.fq")
    files = {"r1": file1, "r2": file2}
    cbu = bc.get_barcodes(
        files, valid_CBC=["TTGGAGAAGTATTCCG", "AAAAAGAAGTATTCCG"]
    ).sort_index()

    # fmt: off
    ind = pd.MultiIndex.from_tuples([
        ('TTGGAGAAGTATTCCG', 'GGGGTGCCGCCATTGGACTCGAGAAGCCGTTGGCAGGGTT', 'CGATCAGATCCC'),
        ('AAAAAGAAGTATTCCG', 'CCCCTGGGGGCATTGGACTCGAGAAGCCGTTGGCAGAAAA', 'CGATCAGATAAA'),
        ('AAAAAGAAGTATTCCG', 'CCCCTGGGGGCATTGGACTCGAGAAGCCGTTGGCAGAAAA', 'CGATCAGATTTT')
        # ('', '', ''),
    ])
    # fmt: on
    ref = pd.Series([4, 5, 3], index=ind)
    ref = ref.sort_index()

    assert cbu.equals(ref)

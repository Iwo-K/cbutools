#!/usr/bin/env ipython
from pathlib import Path
import bctools as bc


HERE = Path(__file__).parent
# assets_path = HERE / "assets/"

def test_load_larry():
    file1 = str(HERE / 'data/SLX-22398.SITTA1.s_1.r_1_small.fq')
    file2 = str(HERE / 'data/SLX-22398.SITTA1.s_1.r_2_small.fq')
    files = {'r1' : file1,
             'r2' : file2}
    bcu = bc.get_barcodes(files)

    # a = 1
    # b = 2
    # assert a == b

def test_BCU_filter_by_reads():
    pass

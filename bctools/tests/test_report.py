from pathlib import Path
import bctools as bc
import os.path

HERE = Path(__file__).parent


def test_larry1_report(tmp_path):
    """Generating and checking the report. Note: report content is not verified"""
    file1 = str(HERE / "data/SLX-22398.SITTA1.s_1.r_1_small.fq")
    file2 = str(HERE / "data/SLX-22398.SITTA1.s_1.r_2_small.fq")

    bc.get_barcodes_report(
        files={"r1": file1, "r2": file2},
        report_template="larry1",
        save_path=tmp_path,
        include_code=False,
    )

    assert os.path.isfile(tmp_path / "larry1.ipynb")
    assert os.path.isfile(tmp_path / "larry1.html")
    assert os.path.isfile(tmp_path / "larry1.csv")

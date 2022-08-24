from .larry import process_larry
from scipy.stats.contingency import crosstab
import numpy as np
import pandas as pd
from . import hamming
from . import cbu
from .cbu import load_barcodes
import papermill as pm
from pathlib import Path
from subprocess import run

HERE = Path(__file__).parent


def get_barcodes(files, type="larry", valid_CBC=None, save="barcodes.csv"):
    """Describe."""
    if type == "larry":
        counts = process_larry(files, valid_CBC=valid_CBC)
        return counts


def get_barcodes_report(
    files=None,
    report_template="larry1",
    save_path=".",
    prefix="larry1",
    make_html=True,
    include_code=True,
):
    if save_path is str:
        save_path = Path(save_path)
    ipynb_path = Path(save_path, prefix).with_suffix(".ipynb")
    csv_path = Path(save_path, prefix).with_suffix(".csv")

    if report_template == "larry1":
        template = str(HERE / "templates/larry1_template.ipynb")
        print("Converting notebook to HTML...")
        pm.execute_notebook(
            template,
            ipynb_path,
            parameters=dict(
                read1=files["r1"],
                read2=files["r2"],
                min_reads=2,
                min_hamming=2,
                min_umis=2,
                csv_path=str(csv_path),
            ),
        )

    print("Converting notebook to HTML...")
    if include_code:
        out = run(
            f"jupyter nbconvert --to html {ipynb_path}",
            capture_output=True,
            text=True,
            shell=True,
        )
    else:
        out = run(
            f"jupyter nbconvert --no-input --to html {ipynb_path}",
            capture_output=True,
            text=True,
            shell=True,
        )
    print(out.stderr)
    print(out.stdout)

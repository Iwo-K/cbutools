#!/usr/bin/env ipython

import numpy as np
import pandas as pd
from subprocess import run
import dnaio
from .cbu import CBUSeries
from tempfile import TemporaryDirectory


def process_larry(files, valid_CBC=None, debug=False):
    # Making a temporary directory
    temp = TemporaryDirectory(prefix=".intermediate_cutadapt_", dir="./")
    tempDIR = temp.name + "/"

    # Trimming the files with cutadapt to retain only larry matching barcodes
    command1 = (
        f"cutadapt -g TTGCTAGGAGAGACCATATG...ATGTCTGGATCCGATATCGC "
        f'-o {tempDIR}bar1.fq -p {tempDIR}cbcumi1.fq {files["r2"]} {files["r1"]} '
        f"--discard-untrimmed --pair-filter=both"
    )
    # Optionally save the report to json format with  --json=test.cutadapt.json
    out1 = run(command1, capture_output=True, text=True, shell=True)
    print(out1.stderr)
    print(out1.stdout)

    command2 = (
        f'cutadapt -g "NNNNTGNNNNCANNNNACNNNNGANNNNGTNNNNAGNNNN;min_overlap=36" '
        f"-o {tempDIR}bar2.fq -p {tempDIR}cbcumi2.fq {tempDIR}bar1.fq {tempDIR}cbcumi1.fq "
        f"--discard-untrimmed --action=retain --pair-filter=both"
    )
    out2 = run(command2, capture_output=True, text=True, shell=True)
    print(out2.stderr)
    print(out2.stdout)

    def read_barCBCUMI(x):
        bar = x[0].sequence
        cbcumi = x[1].sequence
        if (len(bar) == 40) and (
            len(cbcumi) == 28
        ):  # Only keeping LARRY barcodes with 40nt and CBC-UMI with 28
            cbc = cbcumi[:16]
            umi = cbcumi[16:28]
            k = (cbc, bar, umi)
            return k
        else:
            return "WRONG_LEN"

    # Accumulating everything in a dictionary with tuples as keys
    barfile = f"{tempDIR}bar2.fq"
    cbufile = f"{tempDIR}cbcumi2.fq"
    if valid_CBC is not None:
        counts = {"WRONG_LEN": 0, "CBC_ABSENT": 0}
        with dnaio.open(file1=barfile, file2=cbufile) as reader:
            for record in reader:
                k = read_barCBCUMI(record)
                if k[0] in valid_CBC:
                    if not k in counts:
                        counts[k] = 0
                    counts[k] += 1
                elif k == "WRONG_LEN":
                    counts["WRONG_LEN"] += 1
                else:
                    counts["CBC_ABSENT"] += 1
    else:
        counts = {"WRONG_LEN": 0, "CBC_ABSENT": 0}
        with dnaio.open(file1=barfile, file2=cbufile) as reader:
            for record in reader:
                k = read_barCBCUMI(record)
                if not k in counts:
                    counts[k] = 0
                counts[k] += 1

    wronglen_count = counts["WRONG_LEN"]
    cbcabsent_count = counts["CBC_ABSENT"]
    del counts["WRONG_LEN"]
    del counts["CBC_ABSENT"]

    index = pd.MultiIndex.from_tuples(
        [i for i in counts.keys()], names=("CBC", "Barcode", "UMI")
    )
    counts = CBUSeries([v for v in counts.values()], index=index)

    total_count = counts.sum() + wronglen_count + cbcabsent_count
    valid_count = counts.sum()
    print("============================================================")
    print("Counting LARRY barcodes and CBC UMIs")
    print(
        f"Reads with wrong length: {wronglen_count} ({wronglen_count/total_count*100:.2f}%)"
    )
    print(
        f"Reads with CBC absent in the list: {cbcabsent_count}"
        f"({cbcabsent_count/total_count*100:.2f}%)"
    )
    print(f"Valid reads: {valid_count} ({valid_count/total_count*100:.2f}%)")

    if debug:
        return counts, temp
    else:
        return counts

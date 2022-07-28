#!/usr/bin/env ipython

import numpy as np
from subprocess import run
import dnaio


def process_larry(files):
    # Trimming the files with cutadapt to retain only larry matching barcodes
    command1 = f'cutadapt -g TTGCTAGGAGAGACCATATG...ATGTCTGGATCCGATATCGC -o bar1.fq -p cbcumi1.fq {files["r2"]} {files["r1"]} --discard-untrimmed --pair-filter=both'
    # Optionally save the report to json format with  --json=test.cutadapt.json
    out1 = run(command1, capture_output=True, text=True, shell=True)
    print(out1.stderr)
    print(out1.stdout)

    command2 = f'cutadapt -g "NNNNTGNNNNCANNNNACNNNNGANNNNGTNNNNAGNNNN;min_overlap=36" -o bar2.fq -p cbcumi2.fq bar1.fq cbcumi1.fq --discard-untrimmed --action=retain --pair-filter=both'
    out2 = run(command2, capture_output=True, text=True, shell=True)
    print(out2.stderr)
    print(out2.stdout)

    # Reading the trimmed fastq files
    # A more efficient way will be to get the preallocate the sizes given the number of reads comming from cutadapt, for now just running a list
    # bars = np.array() #allocate the size
    # cbcumis = np.array()

    bars = []
    cbcumis = []

    with dnaio.open(file1="bar2.fq", file2="cbcumi2.fq") as reader:
        for record in reader:
            bars.append(record[0].sequence)
            cbcumis.append(record[1].sequence)
    print(len(bars))
    print(len(cbcumis))

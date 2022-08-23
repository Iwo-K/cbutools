#!/usr/bin/env ipython

from .larry import process_larry
from scipy.stats.contingency import crosstab
import numpy as np
import pandas as pd
from . import hamming
from . import cbu
from .cbu import load_barcodes


def get_barcodes(files, type="larry", valid_CBC=None, save="barcodes.csv"):
    """Describe."""
    if type == "larry":
        counts = process_larry(files, valid_CBC=valid_CBC)
        return counts

    # Save intermediate file (just a csv?)

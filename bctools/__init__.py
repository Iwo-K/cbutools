#!/usr/bin/env ipython

from .larry import process_larry
from scipy.stats.contingency import crosstab
import numpy as np
import pandas as pd
from . import hamming
from . import classes

def get_barcodes(files, type='larry', valid_CBC=None, save='barcodes.csv'):
    if type == 'larry':
        counts = process_larry(files)
        return counts

    #Save intermediate file (just a csv?)

    #Extract CBCs and UMIs - make a 3d array with counts?

def count_UMIs(counts):
    pass


def filter_byUMI(counts):
    pass

def assign_barcodes(counts):
    pass

#Can just use subsetting
# def filter_byreads(counts):
#     pass

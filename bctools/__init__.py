#!/usr/bin/env ipython

from .larry import process_larry
from scipy.stats.contingency import crosstab
import numpy as np
import pandas as pd

def get_barcodes(files, type='larry', valid_CBC=None, save='barcodes.csv'):
    if type == 'larry':
        bar, cbc, umi = process_larry(files)
        (barI, cbcI, umiI), values = crosstab(bar, cbc, umi)
        cbu = CBU(bar=barI, cbc=cbcI, umi=umiI, counts=values)
        return cbu

    #Save intermediate file (just a csv?)

    #Extract CBCs and UMIs - make a 3d array with counts?

class CBU:
    def __init__(self, bar, cbc, umi, counts):
        self.bar = bar
        self.cbc = cbc
        self.umi = umi
        self.counts = counts

    def __repr__(self):
        text = f"""CBU (CBC-Barcode-UMI) object with:
        {len(self.cbc)} CBCs
        {len(self.bar)} barcodes
        {len(self.umi)} UMIs
        Dimensions: {self.counts.shape}
        Non zero entires: {np.count_nonzero(self.counts)}"""
        return text

    def __print__(self):
        return self.__repr__()

    def summary(self):
        print('To be implemented')

    def filter_by_reads(self, min_reads=0):
        #1. reduce entires below threshold to 0
        #2. exclude redundant in each axes

        self.counts[self.counts < min_reads] = 0
        #TODO TO BE FINISHED

        #Removing redundant barcodes, CBCs, UMIs
        # sub = aa

    def __getitem__(self, arg):
        # TODO implement subsetting with boolean
        pass


def count_reads(bcfile):
    pass

def filter_byreads(counts):
    pass

def count_UMIs(counts):
    pass

def filter_byhamming(counts):
    pass

def filter_byUMI(counts):
    pass

def assign_barcodes(counts):
    pass

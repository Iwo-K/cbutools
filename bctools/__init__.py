#!/usr/bin/env ipython

from .larry import process_larry

def get_barcodes(files, type='larry', valid_CBC=None, save='barcodes.csv'):
    if type == 'larry':
        bcs = process_larry(files)

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

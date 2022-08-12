import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .hamming import hamming_filter
from math import ceil
from functools import wraps

# TODO how to add check that the indices are names correctly and in the right order? (maybe add an if clause and then the inex check?)
# maybe decorators? I have failed so far in the __init__ modification method
#
# TODO add simple summary print (no of barcodes, no of barcodes per cell with some deviations?)
# probably best if this summary is stored but then updated with every filtering step

def load_barcodes(file):
    """Loads a .csv files with barcodes and convert automatically to CBUSeries or CBseries"""

    data = pd.read_csv(file)
    cols = data.columns.tolist()
    if cols == ['CBC', 'Barcode', 'UMI', '0']:
        ser = pd.Series(data['0'])
        ind =  pd.MultiIndex.from_frame(data.loc[:, ['CBC', 'Barcode', 'UMI']])
        ser.index = ind
        return CBUSeries(ser)

    elif cols == ['CBC', 'Barcode', '0']:
        ser = pd.Series(data['0'])
        ind =  pd.MultiIndex.from_frame(data.loc[:, ['CBC', 'Barcode']])
        ser.index = ind
        return CBSeries(ser)

    else:
        raise Exception('This does not look like CBUSeries or CBSeries saved')


def check_CB_index(f):
    """Decorator to check the index before running some CBU-specific methods"""
    @wraps(f)
    def wrapper(*args, **kwargs):
        if args[0].index.names != ['CBC', 'Barcode']:
            raise Exception('The multiindex is not "CBC", "Barcode"')
        return f(*args, **kwargs)
    return wrapper


def filter_series(series, groupby=None, labels=["CBC", "Barcode", "UMI"], min_counts=0):
    """Filter a series (pd.Series, CBSeries or CBUSeries) based on counts,
    grouped by level or level combination of multiindex"""

    # If none or all three levels are supplied - just filter the values and return
    if (groupby == None) or sorted(groupby) == sorted(labels):
        return series[series >= min_counts]

    # Converting a string into a list
    if type(groupby) is str:
        groupby = [groupby]

    # Catching invalid groupby
    if (type(groupby) is list) and (0 < len(groupby) <= len(labels)):
        # Checking if all groupby elements are valid
        check = [i in labels for i in groupby]
        if not all(check):
            raise Exception(
                "Invalid groupby argument (should be None or combination of f{labels}"
            )
    else:
        raise Exception(
            "Invalid groupby argument (should be None or combination of f{labels}"
        )

    groupby = [x for x in labels if x in groupby]
    filtered = series.groupby(groupby).sum()
    filtered = filtered[filtered >= min_counts]

    # Getting the matching index in the original data (dropping levels that were grouped)
    todrop = [i for i in labels if i not in groupby]

    indexSUB = series.index.droplevel(level=todrop)
    return series[indexSUB.isin(filtered.index)]

def plot_groupby_hist(series, groupby, bins=50, vmax=None):
    """Plots histogram of reads, grouped as indicated in the groupby argument"""

    grouped = series.groupby(groupby).sum()
    plt.hist(np.log10(grouped), bins=bins)

    if vmax is None:
        vmax = grouped.max()
    ticks = ceil(np.log10(vmax)) + 1
    vmax = 10**(ticks-1)
    logpos = np.logspace(0, np.log10(vmax), ticks)
    plt.xticks(range(ticks), logpos)
    # plt.xscale('log')
    plt.yscale("log")
    plt.title(f"Number of reads per {groupby} combination")
    plt.show()

class CBSeries(pd.Series):
    def __init__(self, *args, **kwargs):
        super(CBSeries, self).__init__(*args, **kwargs)
        if not self.index.is_unique:
            raise Exception("Index is not unique")

    @property
    def _constructor(self):
        if not self.index.is_unique:
            raise Exception("Index is not unique")
        return CBSeries

    @check_CB_index
    def filter_by_UMI(self, groupby="Barcode", min_counts=0):
        labels = ["CBC", "Barcode"]
        return filter_series(self, groupby=groupby, labels=labels, min_counts=min_counts)

    def plot_hist(self, groupby='Barcode', *args, **kwargs):
        plot_groupby_hist(self, groupby=groupby, *args, **kwargs)

    @check_CB_index
    def assign_barcodes(self, dispr_filter=None):
        """Assigns barcode"""

        df = pd.DataFrame()

        for i in self.index.levels[0]:
            counts = self[i]
            if dispr_filter is not None:
                max_counts = max(counts)
                counts = counts[counts > (dispr_filter * max_counts)]
            barcodes = counts.index.tolist()
            row = pd.DataFrame({'Barcode_list' : [barcodes], 'Barcode_n' : len(barcodes)}, index=[i])
            df = pd.concat([df, (row)])

        def catl(x):
            ''' Function for concatenating strings '''
            return('-'.join(x))
        df['Barcode'] = df['Barcode_list'].apply(func = catl)
        return df

    def save_barcodes(self, *args, **kwargs):
        self.to_csv(*args, **kwargs)


def check_CBU_index(f):
    """Decorator to check the index before running some CBU-specific methods"""
    @wraps(f)
    def wrapper(*args, **kwargs):
        if args[0].index.names != ['CBC', 'Barcode', 'UMI']:
            raise Exception('The multiindex is not "CBC", "Barcode", "UMI"')
        return f(*args, **kwargs)
    return wrapper


class CBUSeries(pd.Series):
    def __init__(self, *args, **kwargs):
        super(CBUSeries, self).__init__(*args, **kwargs)
        if not self.index.is_unique:
            raise Exception("Index is not unique")

    @property
    def _constructor(self):
        return CBUSeries

    def plot_hist(self, groupby='Barcode', *args, **kwargs):
        plot_groupby_hist(self, groupby=groupby, *args, **kwargs)


    @check_CBU_index
    def filter_by_reads(self, groupby="Barcode", min_counts=0):
        """Filters by minimum expected number of reads. With the groupby column, level combinations are pertmitted"""
        labels = ["CBC", "Barcode", "UMI"]
        return filter_series(self, groupby=groupby, labels=labels, min_counts=min_counts)

    @check_CBU_index
    def filter_by_hamming(self, which="Barcode", min_distance=2):

        counts = self.groupby([which]).sum()
        tokeep, toreject, ties = hamming_filter(counts=counts, min_distance=min_distance)
        if len(ties) > 0:
            print(f"Ties detected between: {ties}")

        if which == 'CBC':
            return self.loc[tokeep, :, :]

        if which == 'Barcode':
            return self.loc[:, tokeep, :]

        if which == 'UMI':
            return self.loc[:, :, tokeep]

    @check_CBU_index
    def count_UMI(self):
        return CBSeries(self.groupby(["CBC", "Barcode"]).size())

    def save_barcodes(self, *args, **kwargs):
        self.to_csv(*args, **kwargs)

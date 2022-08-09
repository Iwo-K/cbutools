import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .hamming import hamming_filter
from math import ceil

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

    def filter_by_UMI(self, groupby="Barcode", min_counts=0):
        labels = ["CBC", "Barcode"]
        return filter_series(self, groupby=groupby, labels=labels, min_counts=min_counts)

    def assign_barcode(self):
        pass

    def plot_hist(self, groupby='Barcode', *args, **kwargs):
        plot_groupby_hist(self, groupby=groupby, *args, **kwargs)


class CBUSeries(pd.Series):
    def __init__(self, *args, **kwargs):
        super(CBUSeries, self).__init__(*args, **kwargs)
        if not self.index.is_unique:
            raise Exception("Index is not unique")
        # Where should these checks be made? Should I make a constructor function?
        # if self.index.nlevels != 3:
        #     raise Exception('Wrong number of index levels')
        # if self.index.names != ['CBC', 'Barcode', 'UMI']:
        #     raise Exception('Index names are not: "CBC", "Barcode, "UMI"')

    # setting this will mean that whenever a new Series is generated it will
    # also be CBUSeries
    @property
    def _constructor(self):
        if not self.index.is_unique:
            raise Exception("Index is not unique")
        if self.index.names == ["CBC", "Barcode", "UMI"]:
            return CBUSeries
        elif self.index.names == ["CBC", "Barcode"]:
            return CBSeries
        else:
            return pd.Series

    def plot_hist(self, groupby='Barcode', *args, **kwargs):
        plot_groupby_hist(self, groupby=groupby, *args, **kwargs)


    def filter_by_reads(self, groupby="Barcode", min_counts=0):
        """Filters by minimum expected number of reads. With the groupby column, level combinations are pertmitted"""
        labels = ["CBC", "Barcode", "UMI"]
        return filter_series(self, groupby=groupby, labels=labels, min_counts=min_counts)

    def filter_by_hamming(self, which="Barcode", threshold=2, collapse=True):

        counts = self.groupby([which]).sum()
        tokeep, torejet, ties = hamming_filter(counts=counts, threshold=threshold)
        if len(ties) > 0:
            print(f"Ties detected between: {ties}")
        return self.loc[:, tokeep, :]

    def count_UMI(self):
        return CBSeries(self.groupby(["CBC", "Barcode"]).size())

    # def filter_by_reads(self, groupby="Barcode", min_counts=0):
    #     """Filters by minimum expected number of reads. With the groupby column, level combinations are pertmitted"""
    #     # Add filtering min counts, min counts per cell, speify the levels for summary
    #     labels = ["CBC", "Barcode", "UMI"]

    #     # If none or all three levels are supplied - quick return
    #     if (groupby == None) or sorted(groupby) == sorted(labels):
    #         return self[self >= min_counts]

    #     # Converting single-element list to a string
    #     if (type(groupby) is list) and len(groupby) == 1:
    #         groupby = groupby[0]

    #     # Getting indices for a single index level (which can be easily extracted with get_level_values())
    #     if type(groupby) is str:
    #         if groupby in labels:
    #             indexSUB = tuple(self.index.get_level_values(groupby))
    #         else:
    #             raise Exception(
    #                 "Invalid groupby argument. Groupby can only be None, string or a list containing combinations of CBC, Barcode and UMI"
    #             )

    #     # Running if groupby is a list
    #     elif type(groupby) is list:
    #         if len(groupby) > 2:  # The 3-element list shoudl be caught above
    #             raise Exception(
    #                 "Invalid groupby argument. Groupby can only be None, string or a list containing combinations of CBC, Barcode and UMI"
    #             )

    #         # Checking if all groupby elements are valid
    #         check = [i in labels for i in groupby]
    #         if not all(check):
    #             raise Exception(
    #                 "Some of groupby arguments are not in the CBC,Barcoce, UMI set. Groupby argument can only be None, string or a list containing combinations of CBC, Barcode and UMI"
    #             )

    #         # Getting indices for >1 index level (as tuples)
    #         groupby = [x for x in labels if x in groupby]
    #         indexSUB = []
    #         for i in groupby:
    #             indexSUB.append(self.index.get_level_values(i))
    #         indexSUB = tuple(zip(*indexSUB))

    #     else:
    #         raise Exception(
    #             "Problem with the groupby argument. Groupby argument can only be None, string or a list containing combinations of CBC, Barcode and UMI"
    #         )

    #     # Filtering and subsetting index for those present in the filtered data
    #     filtered = self.groupby(groupby).sum()
    #     filtered = filtered[filtered >= min_counts]

    #     tokeep = [i in filtered.index for i in indexSUB]
    #     return self[tokeep]

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .hamming import hamming_filter

class CBUSeries(pd.Series):
    def __init__(self, *args, **kwargs):
        super(CBUSeries, self).__init__(*args, **kwargs)
        if not self.index.is_unique:
            raise Exception('Index is not unique)')
        # These checks should be moved to a constructor function maybe?
        # if self.index.nlevels != 3:
        #     raise Exception('Wrong number of index levels')
        # if self.index.names != ['CBC', 'Barcode', 'UMI']:
        #     raise Exception('Index names are not: "CBC", "Barcode, "UMI"')

    # setting this will mean that whenever a new Series is generated it will
    # also be CBUSeries
    @property
    def _constructor(self):
        if not self.index.is_unique:
            raise Exception('Index is not unique)')
        return CBUSeries

    def plot_reads_histogram(self, groupby='Barcode'):
        """ Plots histogra of reads, grouped as indicated in the groupby argument"""
        labels = ['CBC', 'Barcode', 'UMI']
        def plot_hist(x, title=''):
            plt.hist(np.log(x)/np.log(10), bins=50)
            plt.xticks(range(5), np.logspace(0,4,5))
            # plt.xscale('log')
            plt.yscale('log')
            plt.title(title + '\nMax reads plotted = 10000')
            plt.show()

        if (groupby == None) or sorted(groupby) == sorted(labels):
            plot_hist(self, title='Number of reads per CBC-Barcode-UMI combination')
            plt.show()

        else:
            grouped = self.groupby(groupby).sum()
            plot_hist(grouped, title=f'Number of reads per {groupby} combination')
            plt.show()


    def filter_by_reads(self, groupby='Barcode', min_counts=0):
        """ Filters by minimum expected number of reads. With the groupby column, level combinations are pertmitted"""
        # Add filtering min counts, min counts per cell, speify the levels for summary
        labels = ['CBC', 'Barcode', 'UMI']

        # If none or all three levels are supplied - quick return
        if (groupby == None) or sorted(groupby) == sorted(labels):
            return self[self >= min_counts]

        # Converting single-element list to a string
        if (type(groupby) is list) and len(groupby) == 1:
            groupby = groupby[0]

        # Getting indices for a single index level (which can be easily extracted with get_level_values())
        if type(groupby) is str:
            if groupby in labels:
                indexSUB = tuple(self.index.get_level_values(groupby))
            else:
                raise Exception('Invalid groupby argument. Groupby can only be None, string or a list containing combinations of CBC, Barcode and UMI')

        # Running if groupby is a list
        elif type(groupby) is list:
            if len(groupby) > 2: #The 3-element list shoudl be caught above
                raise Exception('Invalid groupby argument. Groupby can only be None, string or a list containing combinations of CBC, Barcode and UMI')

            # Checking if all groupby elements are valid
            check = [i in labels for i in groupby]
            if not all(check):
                raise Exception('Some of groupby arguments are not in the CBC,Barcoce, UMI set. Groupby argument can only be None, string or a list containing combinations of CBC, Barcode and UMI')

            # Getting indices for >1 index level (as tuples)
            groupby = [x for x in labels if x in groupby]
            indexSUB = []
            for i in groupby:
                indexSUB.append(self.index.get_level_values(i))
            indexSUB = tuple(zip(*indexSUB))

        else:
            raise Exception('Problem with the groupby argument. Groupby argument can only be None, string or a list containing combinations of CBC, Barcode and UMI')

        # Filtering and subsetting index for those present in the filtered data
        filtered = self.groupby(groupby).sum()
        filtered = filtered[filtered >= min_counts]

        tokeep = [i in filtered.index for i in indexSUB]
        return self[tokeep]

    def filter_by_hamming(self, which='Barcode', threshold=2, collapse=True):

        counts = self.groupby([which]).sum()
        tokeep, torejet, ties =  hamming_filter(counts=counts,
                                                threshold=threshold)
        if len(ties) > 0:
            print(f"Ties detected between: {ties}")
        return self.loc[:,tokeep,:]

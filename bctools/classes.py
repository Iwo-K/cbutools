import pandas as pd
import matplotlib.pyplot as plt

class CBUSeries(pd.Series):
    def __init__(self, *args, **kwargs):
        super(CBUSeries, self).__init__(*args, **kwargs)
        # These checks should be moved to a constructor function maybe?
        # if self.index.nlevels != 3:
        #     raise Exception('Wrong number of index levels')
        # if self.index.names != ['CBC', 'Barcode', 'UMI']:
        #     raise Exception('Index names are not: "CBC", "Barcode, "UMI"')

    # setting this will mean that whenever a new Series is generated it will
    # also be CBUSeries
    @property
    def _constructor(self):
        return CBUSeries

    def plot_read_histogram(self):
        plt.hist(self, bins=30)
        plt.yscale('log')
        plt.xscale('log')
        plt.title('Number of reads per CBC-Barcode-UMI combination')
        plt.show()

    def filter_by_reads(self, levels='Barcode', min_counts=0, min_counts_per_cell=0):
        #Add filtering min counts, min counts per cell, speify the levels for summary
        return self[self >= min_counts]

    def filter_by_hamming(self, which='barcode', hamming_distance=2, collapse=True):
        """
        Looks at pairwise hamming distances between indicated barcodes and
        if barcodes differ by <= hamming_distance, either keep the barcode with more
        counts (collapse=True) or sums the counts for both
        """
        def hamming(a,b):
            return np.sum([x1 != x2 for x1,x2 in zip(a, b)])

        filtered = {}
        for i in self.index.get_level_values(which):
            found = False
            for j in filtered:
                #Add to the filtered list

                #if found then decide what to do with and
                #Found=True
                break

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.core.debugger import set_trace

def map2d(x, y, func, symmetry=False):
    """Calls function elementwise on iterable 1d objects x and y and returns
    a matrix with the output values.
    By setting symmetry to True, function assumes symmetry in the output matrix
    thus only calculates the upper triangle and copying results into the lower
    triangle (thus cutting the computation time in half). Note: this involves
    computing the transpose so may consume more memory"""
    mat = np.zeros((len(x), len(y)))

    if not symmetry:
        for col in range(0, len(x)):
            for row in range(0, len(y)):
                 mat[row, col] = func(x[row], y[col])
    else:
        n = 0
        for col in range(0, len(x)):
            n += 1
            for row in range(n):
                mat[row, col] = func(x[row], y[col])
        mat = mat + np.triu(mat, k=1).T
    return mat


def hamming_distance(s1, s2):
    """Returns the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))

def filter_by_hamming(cbu, which='Barcode', threshold = 5):
    """Filters out barcodes based on Hamming distance.

    Calculates hamming distances between each pair of barcodes, then removes
    from each pair with distance below the threshold the less abundant barcode.

    Parameters
    ----------
    cbu : CBUSeries object
        child of pd.Series with multiindex: CBC, Barcode and UMI
    threshold : int
        Minimum hamming distance threshold, barcodes with distance < threshold are discarded
    Returns
    -------
    CBUSeries
        Filtered CBUSeries

    """
    seqs = cbu.index.get_level_values(which).unique().values
    seqs = seqs.astype(str)
    counts = cbu.groupby([which]).sum()

    # Computing pairwise hamming distances
    hdist = map2d(seqs, seqs, hamming_distance, symmetry=True)

    # Printing histogram for to see hamming distances
    n, bins, patches = plt.hist(hdist[np.triu_indices_from(hdist, k=1)], bins=50)
    plt.title('Pairwise hamming distance')
    plt.show()

    belowtr = np.nonzero(hdist < threshold)
    belowtr = np.array(belowtr).T

    toreject = []
    ties=[]
    for x1, x2 in belowtr:
        if(x1 != x2):
            x1 = seqs[x1]
            x1val = counts[x1]
            x2 = seqs[x2]
            x2val = counts[x2]
            if x1val < x2val: toreject.append(x1)
            elif x2val < x1val: toreject.append(x2)
            else:
                #How to resolve ties?
                ties.append([x1, x2])

    if len(toreject) == 0:
        print(f'No barcodes rejected\nNumber of ties: {ties}')
        return cbu
    else:
        toreject = pd.Series(toreject).unique()
        print(f"Rejected: {len(toreject)} barcodes\nNumber of ties: {ties}")
        return cbu[~cbu.index.get_level_values(which).isin(toreject)]

# solution from the LARRY github repo
# Sorts all barcodes and loop through barcodes
# The first caught barcode goes into the good_gfp_bcs dictionary
# and then each next barcode is first compared against good_gfp_bcs and
# if not matching then added as well
# This is a very performant solution as it does not do all the pairwise comaprisons
# but it depends on the sorted order, which means
# that it does not provide the actual true barcode sequence (gets collapsed to alphabeticaly first)
# may also cause some problems in edge cases with e.g. three barcodes a,b,c with a-b and b-c differing by 2 letters
# but a-c differing by let's say 3 letters
# Also does not store distances for inspection
#
# N_HAMMING = 2
# def hamming(bc1,bc2): return np.sum([x1 != x2 for x1,x2 in zip(bc1,bc2)])

# all_gfp_bcs = sorted(set([k[1] for k in gz.index.to_list()]))
# good_gfp_bcs = []
# bc_map = {}
# for i,bc1 in enumerate(all_gfp_bcs):
#     if i > 0 and i % 500 == 0: print('Mapped '+repr(i)+' out of '+repr(len(all_gfp_bcs))+' barcodes')
#     mapped = False
#     for bc2 in good_gfp_bcs:
#         if hamming(bc1,bc2) <= N_HAMMING:
#             mapped = True
#             bc_map[bc1] = bc2
#             break
#     if not mapped:
#         good_gfp_bcs.append(bc1)

# print('\nCollapsed '+repr(len(bc_map))+' barcodes')
# for bar in good_gfp_bcs: bc_map[bar] = bar

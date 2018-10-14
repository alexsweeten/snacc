import pandas as pd
import numpy as np


def read_dist(csv_file):
    df = pd.read_csv(csv_file, index_col=0, sep=None)
    return df.values


def read_dist_dataframe(csv_file):
    df = pd.read_csv(csv_file, index_col=0, sep=None)
    return df


def read_dist_values_names(csv_file):
    df = pd.read_csv(csv_file, index_col=0, sep=None)
    return df.index.values, df.values

    
def metrify(D):
    # Converts a raw distance matrix to one which
    # is symmetric and has zero diagonals
    D_sym = 0.5*(D + D.T)
    np.fill_diagonal(D_sym, 0.)
    return D_sym

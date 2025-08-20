import h5py
import numpy as np
from scipy.sparse import coo_matrix

with h5py.File("hamiltonians/example.h5", "r") as f:
    vals = f["ham_dp"][:]     # or combine ham_dc_real/imag
    row  = f["rc"][:,0] - 1   # adjust to 0-based
    col  = f["rc"][:,1] - 1
    H = coo_matrix((vals, (row, col)))

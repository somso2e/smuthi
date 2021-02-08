import numpy as np
from numba import complex64, complex128, float64, prange, jit


@jit(['complex64[:, :](complex64[:, :, :], float64[:])', 'complex128[:, :](complex128[:, :, :], float64[:])'],
     nopython=True, cache=True, nogil=True, parallel=True, fastmath=True)
def numba_trapz_3dim_array(y, x):
    shape = y.shape
    result = np.zeros((shape[0], shape[1]), dtype=type(y[0, 0, 0]))
    for i in prange(shape[0]):
        for j in range(shape[1]):
            out = .0 + .0j
            for k in range(len(y[i, j]) - 1):
                out += (x[k+1]-x[k]) * (y[i, j, k+1] + y[i, j, k])/2.0
                # here is still a room for optimization, I guess
            result[i, j] = out

    return result
import numpy as np
from numba import complex64, complex128, float64, prange, jit


@jit(['complex64[:, :](complex64[:, :, :], float64[:])', 
        'complex128[:, :](complex128[:, :, :], float64[:])'],
     nopython=True, cache=True, nogil=True, parallel=True, fastmath=True)
def numba_trapz_3dim_array(y, x):
    """
    This function can replace snippet

    'foo = np.trapz(y, x)'

    by

    'foo = numba_trapz_3dim_array(y, x)'
    """
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

@jit(['complex64[:, :, :] (float32[:], float32[:], float32[:], \
                                        complex64[:, :], complex64[:, :], complex64[:, :])',
      'complex128[:, :, :] (float64[:], float64[:], float64[:], \
                                        complex128[:, :], complex128[:, :], complex128[:, :])'],
     nopython=True, cache=True, nogil=True, parallel=True, fastmath=True)
def numba_3tensordots_1dim_times_2dim(x_float_1dim, y_float_1dim, z_float_1dim, \
                        x_complex_2dim, y_complex_2dim, z_complex_2dim):
    """
    This function can replace snippet

    'foo = np.tensordot(x_float_1dim, x_complex_2dim, axes=0) \n
    foo += np.tensordot(y_float_1dim, y_complex_2dim, axes=0) \n
    foo += np.tensordot(z_float_1dim, z_complex_2dim, axes=0)'

    by

    'foo = get_3_tensordots(x_float_1dim, y_float_1dim, z_float_1dim,
                                    x_complex_2dim, y_complex_2dim, z_complex_2dim)'
    """
    first_shape = x_float_1dim.shape
    second_shape = x_complex_2dim.shape
    result = np.zeros((first_shape[0], second_shape[0], second_shape[1]), dtype=type(x_complex_2dim[0, 0]))
    for i in prange(first_shape[0]):
        for j in range(second_shape[0]):
            for k in range(second_shape[1]):
                result[i, j, k] = x_float_1dim[i] * x_complex_2dim[j, k] + \
                                    y_float_1dim[i] * y_complex_2dim[j, k] + \
                                        z_float_1dim[i] * z_complex_2dim[j, k]

    return result
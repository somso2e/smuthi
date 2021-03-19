import numpy as np
from numba import complex64, complex128, float64, prange, jit, int32, types


@jit(['complex64[:](complex64[:, :, :], float64[:], int32)',
        'complex128[:](complex128[:, :, :], float64[:], int32)'],
     nopython=True, cache=True, nogil=True, fastmath=True, parallel=False)
def __numba_trapz_3dim_array_linear_part(y, x, i):
    one_dim_result = np.zeros((y.shape[1]), dtype=type(y[0, 0, 0]))
    for j in range(y.shape[1]):
        out = .0 + .0j
        for k in range(len(y[i, j]) - 1):
            out += (x[k+1]-x[k]) * (y[i, j, k+1] + y[i, j, k])/2.0
        one_dim_result[j] = out
    
    return one_dim_result


@jit(['complex64[:, :] (float32[:], float32[:], float32[:], \
                                        complex64[:, :], complex64[:, :], complex64[:, :], int32)',
      'complex128[:, :] (float64[:], float64[:], float64[:], \
                                        complex128[:, :], complex128[:, :], complex128[:, :], int32)'],
     nopython=True, cache=True, nogil=True, fastmath=True, parallel=False)
def __numba_3tensordots_1dim_times_2dim_linear_part(x_float_1dim, y_float_1dim, z_float_1dim, \
                        x_complex_2dim, y_complex_2dim, z_complex_2dim, i):
    second_shape = x_complex_2dim.shape
    two_dim_result = np.zeros((second_shape[0], second_shape[1]), dtype=type(x_complex_2dim[0, 0]))

    for j in range(second_shape[0]):
            for k in range(second_shape[1]):
                two_dim_result[j, k] = x_float_1dim[i] * x_complex_2dim[j, k] + \
                                    y_float_1dim[i] * y_complex_2dim[j, k] + \
                                        z_float_1dim[i] * z_complex_2dim[j, k]

    return two_dim_result


@jit(['complex64[:, :](complex64[:, :, :], float64[:])',
        'complex128[:, :](complex128[:, :, :], float64[:])'],
     nopython=True, cache=True, nogil=True, fastmath=True, parallel=False)
def numba_trapz_3dim_array(y, x):
    """
    This function can replace snippet

    'foo = np.trapz(y, x)'

    by

    'foo = numba_trapz_3dim_array(y, x)'
    """
    result = np.zeros((y.shape[0], y.shape[1]), dtype=type(y[0, 0, 0]))
    for i in range(y.shape[0]):
        result[i] =__numba_trapz_3dim_array_linear_part(y, x, i)

    return result


@jit(['complex64[:, :](complex64[:, :, :], float64[:])',
        'complex128[:, :](complex128[:, :, :], float64[:])'],
     nopython=True, cache=True, nogil=True, fastmath=True, parallel=True)
def numba_trapz_3dim_array_parallel(y, x):
    result = np.zeros((y.shape[0], y.shape[1]), dtype=type(y[0, 0, 0]))
    for i in prange(y.shape[0]):
        result[i] =__numba_trapz_3dim_array_linear_part(y, x, i)

    return result



@jit(['complex64[:, :, :] (float32[:], float32[:], float32[:], \
                                        complex64[:, :], complex64[:, :], complex64[:, :])',
      'complex128[:, :, :] (float64[:], float64[:], float64[:], \
                                        complex128[:, :], complex128[:, :], complex128[:, :])'],
     nopython=True, cache=True, nogil=True, fastmath=True, parallel=False)
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
    for i in range(first_shape[0]):
        result[i] = __numba_3tensordots_1dim_times_2dim_linear_part(x_float_1dim, y_float_1dim, z_float_1dim, \
                        x_complex_2dim, y_complex_2dim, z_complex_2dim, i)

    return result


@jit(['complex64[:, :, :] (float32[:], float32[:], float32[:], \
                                        complex64[:, :], complex64[:, :], complex64[:, :])',
      'complex128[:, :, :] (float64[:], float64[:], float64[:], \
                                        complex128[:, :], complex128[:, :], complex128[:, :])'],
     nopython=True, cache=True, nogil=True, fastmath=True, parallel=True)
def numba_3tensordots_1dim_times_2dim_parallel(x_float_1dim, y_float_1dim, z_float_1dim, \
                        x_complex_2dim, y_complex_2dim, z_complex_2dim):
    first_shape = x_float_1dim.shape
    second_shape = x_complex_2dim.shape
    result = np.zeros((first_shape[0], second_shape[0], second_shape[1]), dtype=type(x_complex_2dim[0, 0]))
    for i in prange(first_shape[0]):
        result[i] = __numba_3tensordots_1dim_times_2dim_linear_part(x_float_1dim, y_float_1dim, z_float_1dim, \
                        x_complex_2dim, y_complex_2dim, z_complex_2dim, i)

    return result


@jit(['types.UniTuple(complex64[:, :, :], 3)(complex64[:, :, :], complex64[:, :, :], complex64[:, :, :], complex64[:, :, :])',
        'types.UniTuple(complex128[:, :, :], 3)(complex128[:, :, :], complex128[:, :, :], complex128[:, :, :], complex128[:, :, :])'],
     nopython=True, cache=True, nogil=True, fastmath=True, parallel=True)
def evaluate_r_times_eikr(foo_x, foo_y, foo_z, exp_feed):
    '''
    Attention! Sometimes this function can decrease speed on 1 core mode.
    Here foo_x, foo_y, foo_z are supposed to be 2dim arrays with [None, :, :].
    This function can replace snippet

    exp_j = np.exp(1j * exp_feed)
    foo_x_eikr = foo_x * exp_j
    foo_y_eikr = foo_y * exp_j
    foo_z_eikr = foo_z * exp_j

    by

    foo_x_eikr, foo_y_eikr, foo_z_eikr = numba_multiple_on_exp(foo_x, foo_y, foo_z, kr).
    '''
    shape = exp_feed.shape
    result_x = np.zeros((shape[0], shape[1], shape[2]), dtype=type(exp_feed[0, 0, 0]))
    result_y = np.zeros((shape[0], shape[1], shape[2]), dtype=type(exp_feed[0, 0, 0]))
    result_z = np.zeros((shape[0], shape[1], shape[2]), dtype=type(exp_feed[0, 0, 0]))
    
    for j in prange(shape[1]):
        for k in range(shape[2]):
            foo_x_slice = foo_x[0, j, k]
            foo_y_slice = foo_y[0, j, k]
            foo_z_slice = foo_z[0, j, k]
            for i in range(shape[0]):
                exp_j = (np.cos(exp_feed[i, j, k].real) + 1.0j*np.sin(exp_feed[i, j, k].real)) * \
                                                            np.exp(-exp_feed[i, j, k].imag)
                result_x[i, j, k] = foo_x_slice * exp_j
                result_y[i, j, k] = foo_y_slice * exp_j
                result_z[i, j, k] = foo_z_slice * exp_j
                            
    return result_x, result_y, result_z
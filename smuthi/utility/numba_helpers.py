import numpy as np
from numba import complex64, complex128, float64, prange, jit, int32, int64, types


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


@jit(['void(complex64[:, :, :], complex64[:, :, :], complex64[:, :, :], complex64[:, :, :], \
            complex64[:, :, :], complex64[:, :, :], complex64[:, :, :], int32, int32)',
        'void(complex128[:, :, :], complex128[:, :, :], complex128[:, :, :], complex128[:, :, :], \
            complex128[:, :, :], complex128[:, :, :], complex128[:, :, :], int32, int32)'],
     nopython=True, cache=True, nogil=True, fastmath=True, parallel=False)
def __evaluate_r_times_eikr_linear_part(foo_x, foo_y, foo_z, exp_feed, result_x, result_y, result_z, j, i_shape):
    foo_x_slice = foo_x[0, j]
    foo_y_slice = foo_y[0, j]
    foo_z_slice = foo_z[0, j]
    for i in range(i_shape):
        exp_j = (np.cos(exp_feed[i, j].real) + 1.0j*np.sin(exp_feed[i, j].real)) * \
                                                    np.exp(-exp_feed[i, j].imag)
        result_x[i, j] = foo_x_slice * exp_j
        result_y[i, j] = foo_y_slice * exp_j
        result_z[i, j] = foo_z_slice * exp_j

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
        result[i] = x_float_1dim[i] * x_complex_2dim + \
                        y_float_1dim[i] * y_complex_2dim + \
                            z_float_1dim[i] * z_complex_2dim

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
        result[i] = x_float_1dim[i] * x_complex_2dim + \
                        y_float_1dim[i] * y_complex_2dim + \
                            z_float_1dim[i] * z_complex_2dim

    return result


@jit(['types.UniTuple(complex64[:, :, :], 3)(complex64[:, :, :], complex64[:, :, :], complex64[:, :, :], complex64[:, :, :])',
        'types.UniTuple(complex128[:, :, :], 3)(complex128[:, :, :], complex128[:, :, :], complex128[:, :, :], complex128[:, :, :])'],
     nopython=True, cache=True, nogil=True, fastmath=True, parallel=False)
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

    for j in range(shape[1]):
        __evaluate_r_times_eikr_linear_part(foo_x, foo_y, foo_z, exp_feed, result_x, result_y, result_z, j, shape[0])

    return result_x, result_y, result_z


@jit(['types.UniTuple(complex64[:, :, :], 3)(complex64[:, :, :], complex64[:, :, :], complex64[:, :, :], complex64[:, :, :])',
        'types.UniTuple(complex128[:, :, :], 3)(complex128[:, :, :], complex128[:, :, :], complex128[:, :, :], complex128[:, :, :])'],
     nopython=True, cache=True, nogil=True, fastmath=True, parallel=True)
def evaluate_r_times_eikr_parallel(foo_x, foo_y, foo_z, exp_feed):
    shape = exp_feed.shape
    result_x = np.zeros((shape[0], shape[1], shape[2]), dtype=type(exp_feed[0, 0, 0]))
    result_y = np.zeros((shape[0], shape[1], shape[2]), dtype=type(exp_feed[0, 0, 0]))
    result_z = np.zeros((shape[0], shape[1], shape[2]), dtype=type(exp_feed[0, 0, 0]))

    for j in prange(shape[1]):
        __evaluate_r_times_eikr_linear_part(foo_x, foo_y, foo_z, exp_feed, result_x, result_y, result_z, j, shape[0])

    return result_x, result_y, result_z


@jit(float64(int64), nopython=True, cache=True, nogil=True)
def jitted_prefactor(m):
    """Returns the prefactor :math:`\sqrt(\frac{(2*m+1)}{2(2m)!} (2m-1)!!`
    without using factorials nor bignum numbers, which makes it jittable.

    Args:
        m (int64): Argument (non-negative)

    Returns:
        :math:`\sqrt(\frac{(2*m+1)!!}{2(2m)!!}`
    """
    res = 1.
    for t in range(2,2*m+2,2):
        res += res/t # equivalent to res *= (t+1)/t, but retains all significant digits
    return (res/2)**0.5


@jit(float64(int64), nopython=True, cache=True, nogil=True)
def prefactor_expansion(m):
    """Expansion of :math:`\sqrt(\frac{(2*m+1)}{2(2m)!} (2m-1)!!` in the limit
    for large :math:`m`.
    The expansion converges very rapidly, so that taking the first 14 terms is
    enough to get every term with double precision down to :math:`m > 10`.
    Further terms are available, but they involve bignum integers.
    Therefore, for terms :math:`m < 10`, it's easier to perform the calculation
    directly.

    Args:
        m (int64): Argument (non-negative)

    Returns:
        :math:`\sqrt(\frac{(2*m+1)!!}{2(2m)!!}`
    """
    if m <= 10:
        return jitted_prefactor(m)
    x = (1/m)**(1/4) / 2
    res = 1/2/x \
        + 3/2*x**3 \
        - 23/4*x**7 \
        + 105/4*x**11 \
        - 1317/16*x**15 \
        + 1053/16*x**19 \
        - 132995/32*x**23 \
        + 3300753/32*x**27 \
        + 24189523/256*x**31 \
        - 7404427407/256*x**35 \
        - 45203760489/512*x**39 \
        + 8818244857071/512*x**43 \
        + 99932439090703/2048*x**47 \
        - 29944926937328991/2048*x**51 \
        # - 173768350561954907/4096*x**55 \
        # + 70443346589090375073/4096*x**59 \
        # + 3299174696912539072131/65536*x**63 \
        # - 1750098951789039471641607/65536*x**67 \
        # - 10306911268683450183973389/131072*x**71 \
        # + 6934593111025074608358220131/131072*x**75 \
        # + 82009190686072900341958975797/524288*x**79 \
        # - 68275536676414175755318233771549/524288*x**83 \
        # - 404866365213072360301912455070149/1048576*x**87 \
        # + 408755654050262791361092019679647223/1048576*x**91 \
        # + 9715746540636352989024348199533962007/8388608*x**95 \
        # - 11697824487409665503661150266687972482035/8388608*x**99
    return res / (np.pi)**(1/4)

"""This module contains helper functions for the evaluation of Ewald lattice sums."""

from numba import njit
import numpy as np


@njit()
def reciprocal_lattice_vec(a1, a2):
    """ 
    Args:
        a1 (numpy.ndarray):     lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):     lattice vector 2 in Carthesian coordinates
    Returns:
        2D lattice vectors in reciprocal space, https://en.wikipedia.org/wiki/Reciprocal_lattice
    """ 
    a3 = np.array([0, 0, 1])
    b1 = 2 * np.pi * np.cross(a2, a3) / np.linalg.det(np.concatenate((a1, a2, a3)).reshape(3, 3).astype(np.float64))
    b2 = 2 * np.pi * np.cross(a3, a1) / np.linalg.det(np.concatenate((a2, a3, a1)).reshape(3, 3).astype(np.float64))    
    return b1[:2], b2[:2]


@njit()
def n1n2_indices(x):
    """ Indices n1, n2 that specify lattice point R = n1*a1 + n2*a2 of unit cells in a distance x around
        the arbitrary central unit cell.   
    Args:
        x (int):                distance x from the arbitrary center unit cell
    Returns:
        n1 (numpy.ndarray):     index of lattice vector a1
        n2 (numpy.ndarray):     index of lattice vector a2
    """
    if x == 0:
        return np.array([0], np.float64), np.array([0], np.float64)
    
    array = np.zeros((2, 8*x), np.float64)   
    idx = 0
    for n1 in [-x, x]:
        for n2 in np.arange(-x, x+1): 
            array[:, idx] = [n1, n2]
            idx += 1
    for n1 in np.arange(-x+1, x):
        for n2 in [-x, x]:
            array[:, idx] = [n1, n2]
            idx += 1
    return array[0, :], array[1, :]


# factorials up to 20!
factorial_table = np.array([1, 1, 2, 6, 24, 120, 720, 5040, 40320,
                   362880, 3628800, 39916800, 479001600,
                   6227020800, 87178291200, 1307674368000,
                   20922789888000, 355687428096000, 6402373705728000,
                   121645100408832000, 2432902008176640000], np.int64)
# x! / 20! up to x = 30
factorial_table_2 = np.array([21, 462, 10626, 255024, 6375600, 165765600,
                              4475671200, 125318793600, 3634245014400,
                              109027350432000], np.int64)

@njit()
def numba_factorial(x):
    """ Returns factorials for x <= 30.
    Args:
        x (numpy.ndarray):      argument array of type np.int64
    Returns:
        res (numpy.ndarray):    2 dimensional array that splits the factorial into
                                two factors x! = x2! * x1! to allow storage as np.int64
                                First, for arguments 20 < x2 <= 30, 21*...*x2
                                Second, for arguments x1 <= 20, 1*2*...*x1
    """
    res = np.zeros((x.shape[0], 2), np.int64)
    for idx, xx in enumerate(x):
        if xx > 30:
            raise ValueError('Only factorials of x <= 30 available!')
        elif xx > 20:
            res[idx, 0] = factorial_table_2[xx - 21]
            res[idx, 1] = factorial_table[20]
        else:
            res[idx, 0] = 1
            res[idx, 1] = factorial_table[xx]
    return res


double_factorial_table = np.array([1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840,
                                   10395, 46080, 135135, 645120, 2027025,
                                   10321920, 34459425, 185794560, 654729075,
                                   3715891200, 13749310575, 81749606400,
                                   316234143225, 1961990553600, 7905853580625,
                                   51011754393600, 213458046676875,
                                   1428329123020800, 6190283353629375, 
                                   42849873690624000], np.int64)

@njit()
def numba_double_factorial(x):
    """ Returns double factorials for x <= 30.
    Args:
        x (numpy.ndarray):      argument array of type np.int64
    """
    if x > 30:
        raise ValueError
    return double_factorial_table[x]


###############################################################################
from numba.extending import get_cython_function_address
import ctypes
_PTR = ctypes.POINTER
_dble = ctypes.c_double
_ptr_dble = _PTR(_dble)


def load_wofz(name):
    try:
        addr = get_cython_function_address(name, "wofz")
    except: return None
    else: return addr


import_wofz_possible_file_paths = ["cython_speedups",
         "smuthi.utility.cython.cython_speedups",
         "utility.cython.cython_speedups"]


def raise_wofz_exception():
    raise Exception('''
                    In smuthi/periodic/ewald_helper.py Numba needs to know the location of the function "wofz".\n
                    This function is imported from the Cython file in smuthi/utility/cython/cython_speedups_api.h.\n
                    We tried three main methods to import this file, but all three failed. :( \n
                    This likely means that you are trying to run Smuthi for a non-standard directory and \n
                    Smuthi is unable to deterine where to find the cython_speedups file. If you need to run \n
                    Smuthi in a non-standard directory, you can add the proper file paths to the import_wofz_possible_file_paths\n
                    variable on line 118 in the ewald_helper file. 
                    
                    Alternatively, this error may occur because cython_speedups_api.h does not exist and needs to be\n
                    recompiled from cython_speedups.pyx. If this is the case, please follow the instructions in the file\n
                    smuthi/utility/cython/compile_speedups_from_c.py or, if compile_speedups.c does not exist \n
                    follow the instructions in compile_speedups_from_pyx.py. 
                    ''')


try: 
    for name in import_wofz_possible_file_paths:
        addr = load_wofz(name)
        if addr != None:
            break    
except:
    raise_wofz_exception()

if addr is None:
    raise_wofz_exception()

functype = ctypes.CFUNCTYPE(None, _dble, _dble, _ptr_dble, _ptr_dble)
wofz_fn = functype(addr)


@njit()
def numba_wofz_complex(x):
    """ Numba compatible version of the Faddevva function - replacing scipy.special.wofz().
        The cwrapper is copied from: https://github.com/numba/numba/issues/3086
    """
    out_real = np.empty(1,dtype=np.float64)
    out_imag = np.empty(1,dtype=np.float64)
    wofz_fn(np.real(x), np.imag(x), out_real.ctypes, out_imag.ctypes)    
    return complex(out_real[0] + 1j * out_imag[0])
###############################################################################


@njit()
def upper_incomplete_gamma_fun(n, x_vec):
    """ Recursion formula of the upper incomplete Gamma function.  
        Kambe, Z. Naturforschg. 22a, 322-330, 1967. (42) und Appendix 2 (A8), (A9)
    Args:
        n (int):                    order
        x_vec (numpy.ndarray):      argument array
    Returns:
        gamma_vec (numpy.ndarray):  values of the upper incomplete gamma function for 0 ... n 
                                    dimension: len(x_vec) times (n+1)          
    """   
    gamma_vec = np.zeros((x_vec.shape[0], n+1), np.complex128)
    for idx_x, x in enumerate(x_vec):         
        arg_x = np.arctan2(x.imag, x.real)
        if abs(arg_x - np.pi) < 1e-08:
            arg_x = -arg_x
            val = np.sqrt(np.abs(x)) * np.exp(0.5j * (arg_x + np.pi))
            w = numba_wofz_complex(val)
            gamma_vec[idx_x, 0] = np.sqrt(np.pi) * np.exp(-x) * w 
            for idx_n in range(1, n+1):
                b = 0.5 - idx_n
                gamma_vec[idx_x, idx_n] = (gamma_vec[idx_x, idx_n-1] + x ** b * np.exp(-x)) / b
        else:
            val = np.sqrt(np.abs(x)) * np.exp(0.5j * (arg_x + np.pi))
            w = numba_wofz_complex(val)
            gamma_vec[idx_x, 0] = np.sqrt(np.pi) * np.exp(-x) * w
            for idx_n in range(1, n+1):
                b = 0.5 - idx_n
                gamma_vec[idx_x, idx_n] = (gamma_vec[idx_x, idx_n-1] - x ** b * np.exp(-x)) / b
    return gamma_vec


@njit()
def int_recursion(L, eta, k, R_vec):
    """ Integral recursion for the Ewald sum's real space summand 
        Kambe, Z. Naturforschg. 22a, 322-330, 1967. Appendix 2, (A11) - (A14)
    Args:
        L (int):                multipole degree
        eta (float):            Ewald sum separation parameter 
        k (complex):            initial field's wavenumber
        R_vec (numpy.ndarray):  real space unit cell displacement vector     
    Returns:
        Ewald sum's real space summand integral 
    """
    shape = R_vec.shape[0]
    alpha = k ** 2 / (4 * eta ** 2)
    I = np.zeros((shape, L+2), np.complex128)
    w = np.zeros(shape, np.complex128)
    for idx in range(shape):
        w[idx] = numba_wofz_complex(np.sqrt(alpha) + 1j * k * R_vec[idx] / (2 * np.sqrt(alpha)))
    exp_term =  np.exp(alpha - (k * R_vec) ** 2 / (4 * alpha))
    
    I[:, 0] = np.sqrt(np.pi) * exp_term * w.imag 
    I[:, 1] = np.sqrt(np.pi) * 2 / (k * R_vec) * exp_term * w.real
    
    for idx in range(2, L+2):
        I[:, idx] = (2 / (k * R_vec)) ** 2 * ((2 * (idx-2) + 1) / 2 * I[:, idx-1] - I[:, idx-2] + alpha ** (-(idx-2) - 0.5) * exp_term)        
    return I[:, -1]


@njit()
def delta_n(n_max, gamma, cz, eta):
    """ Integral recursion for the Ewald sum's reciprocal space summand to account for the
        coupling between one particle (S_0i) and a different particle's (S_0j) periodic arrangement.
        Kambe, Z. Naturforschg. 23a, 1280-1294, 1968. Appendix 3, (A3.1) - (A.3.7)
    Args:
        n_max (int):            maximal order
        gamma (numpy.ndarray):  sqrt(k ** 2 - kgt ** 2)
        cz (float):             particle displacement in along z-axis
        eta (float):            Ewald sum separation parameter 
    Returns:
        Ewald sum's reciprocal space summand intergal
    """
    shape = gamma.shape[0]
    I = np.zeros((shape, n_max + 1), np.complex128)
    w_mn = np.zeros(shape, np.complex128)
    w_pl = np.zeros(shape, np.complex128)
    
    x = -gamma ** 2 / (4 * eta ** 2)
    rootx = np.zeros(gamma.shape[0], np.complex128)
    z = np.zeros(gamma.shape[0], np.complex128)     
    
    sign = np.where(x.real < 0, 1.0, 0.0)
    rootx = -1j * (sign * np.abs(x)) ** 0.5 + ((1 - sign) * x) ** 0.5
    z = sign * (gamma * cz) + (1 - sign) * 1j * np.abs(gamma * cz)
    
    zsqr = (gamma * cz) ** 2
    exp_term = np.exp(-x + zsqr / (4 * x))
    for idx in range(shape):
        w_mn[idx] = numba_wofz_complex(-z[idx] / (2 * rootx[idx]) + 1j * rootx[idx])
        w_pl[idx] = numba_wofz_complex(z[idx] / (2 * rootx[idx]) + 1j * rootx[idx])
            
    I[:, 0] = np.pi ** 0.5 / 2 * exp_term * (w_mn + w_pl)
    if n_max == 0:
        return I
    I[:, 1] = 1j * np.pi ** 0.5 / z * exp_term * (w_mn - w_pl)
    
    for idx in range(2, n_max + 1):
        I[:, idx] = 4 / zsqr * ((1.5 - idx) * I[:, idx-1] - I[:, idx-2] + rootx * x ** (1 - idx) * exp_term) 
    return I
        
        
@njit()
def normalization_spherical_harmonics(M):  
    """ Normalization factor between the definition of spherical harmonics in Smuthi
        (following Doicu, Light Scattering by Systems of Particles, Berlin, 2006.) and the definitions
        in Kambe and Beutel (following Tsang, Theory of Microwave Remote Sensing, 1985.)
        for details see: https://doi.org/10.5445/IR/1000141626 (Appendix C, p.135).
    """
    if M >= 0:
        return np.sqrt(2 * np.pi) * (-1) ** (-M)
    else:
        return np.sqrt(2 * np.pi)
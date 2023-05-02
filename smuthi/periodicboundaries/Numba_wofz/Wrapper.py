from numba.extending import get_cython_function_address
from numba import vectorize, njit
import ctypes
import numpy as np
import numba as nb

_PTR = ctypes.POINTER
_dble = ctypes.c_double
_ptr_dble = _PTR(_dble)

addr = get_cython_function_address("special", "wofz")
functype = ctypes.CFUNCTYPE(None, _dble, _dble, _ptr_dble, _ptr_dble)
wofz_fn = functype(addr)


@nb.njit
def numba_wofz_complex(x):
    out_real = np.empty(1,dtype=np.float64)
    out_imag = np.empty(1,dtype=np.float64)
    wofz_fn(np.real(x), np.imag(x), out_real.ctypes, out_imag.ctypes)
    
    return complex(out_real[0] + 1j * out_imag[0])


from scipy.special import wofz

val=complex(1+2j)
expected = wofz(val)
got = numba_wofz_complex(val)

np.testing.assert_allclose(got, expected)
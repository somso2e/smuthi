#!python
#cython: language_level=3
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

"""
Created on Tue Apr  5 10:27:28 2022

This code uses Cython to speed up multiple functions commonly used in Smuthi.
The code also has functions for calculating the explicit direct particle-particle
coupling matrix within a single layer. A hash table method is used to speed 
up the calculation of the explicit direct coupling matrix in both the 2D and 3D 
case. T

NOTE: The code above this comment section is important. 
cython: language_level=3 -> This enforces python3 standards.
Most importantly 1/3 = 0.33333. language_level=2 -> 1/3 = 0

@author: parkerwray
"""

import scipy.special
import numpy as np

# python imports before cython imports, so cython overwrites.
cimport scipy.special.cython_special
cimport numpy as np
cimport cython
from cython.parallel import prange

## Note: All C libraries can be found for Cython in https://github.com/cython/cython/tree/master/Cython/Includes
cdef extern from "complex.h" nogil:
    double complex cexp(double complex)
    double complex cpow(double complex x, double complex y)
    double complex csqrt(double complex z)
    double cabs(double complex z)
    double cimag(double complex z)
    double creal(double complex z)    
    
cdef extern from "<math.h>" nogil:
    double fabs(double x)
    double sqrt(double x)
    double atan2(double y, double x)
    double pow(double x, double y)
    
cdef extern from "<stdlib.h>" nogil:
    long abs (long number)
    
    
    
    


##############################################################################
# V1: Implemented all methods to do 2D coupling (monolayer particle coupling)
# using C, both directly and with look up tables. Also C version of particle 
# outgoing-to-outgoing translations
#
# V2: Implemented methods to do 3D coupling. Differentiated between a 2D and 3D
# Legendre and direct coupling functions. Look up tables are NOT available for
# 3D and therefore, explicitly setting the name to 2D was not done. 
#
# V3: Minor changes that I no longer remember... :(
#
# V4: Major changes to all functions such that the Gil is released (at least 
# partially) in all cdef functions. This is good manners, especially for 
# programs like Dask, that check the status of a python program and does not
# know how to interpret when a Cython function holds the Gil then disapears. 
# https://stackoverflow.com/questions/54769782/distributed-worker-warning-heartbeat-to-scheduler-failed
# The release of the Gil also allows multithreading, which is implemented using
# Cython's prange function.
# https://cython.readthedocs.io/en/latest/src/userguide/parallelism.html
# https://docs.dask.org/en/stable/best-practices.html
#
# V5: Minor bug fixes comparing python to cython. Note: AB coeffs are now 
# enforced real/imag. All int -> long
# 
# V6: Proved the ability to call Cython functions from Numba using the "api" command.
# This was shown on the wofz function. 
#
# Cython vs. Python implementations were tested for speed and accuracy. 
# Values output by the Cython funcitons sub-functions match Python sub-functions
# up to the 16 decimal point.
# Cython can give a speed boost from 80-800x, dependent on what functions 
# are being compared and the state of the machine at runtime. 
##############################################################################



'''
##############################################################################
#****************************************************************************#
#****************************************************************************#
#******************** Sub-functions accelerated with C **********************#
#****************************************************************************#
#****************************************************************************#
##############################################################################

This section details all the subfunctions.
'''

##############################################################################
# Calculate wofz special function in a numba-compatible fashion
##############################################################################

# Numba  requires a compatible version of the Faddevva function - replacing scipy.special.wofz().
# The cwrapper is used in smuthi/periodicboundaries/ewald_helper.py and 
# is copied from: https://github.com/numba/numba/issues/3086
# The cython portion of this solution is found below. The purpose is to have a 
# singular file housing all cython functions because this is easier to track/manage.
# This file is automatically compiled for the user and the "api" generates 
# a header file that numba can find. 


cdef api wofz(double in1_real,double in1_imag,double *out_real,double *out_imag):
  cdef double complex z
  z.real=in1_real
  z.imag=in1_imag
  cdef double complex out = scipy.special.cython_special.wofz(z)
  out_real[0]=out.real
  out_imag[0]=out.imag
  return



##############################################################################
# Recursively calculate factorial (n!) - Speed ~2x when called from Python wrap
##############################################################################
cdef double c_factorial(double n) nogil:
    if n == 0:
        return 1
    elif n < 0:
        return 0 #should be NAN!
    else:
        return n * c_factorial(n-1)

def factorial(double n):
    cdef double val
    with nogil:
         val = c_factorial(n)
    return val

##############################################################################
# Recursively calculate double factorial (n!!) - Speed ~1x when called from Python wrap
##############################################################################  
cdef double c_double_factorial(double n) nogil:
    if n in (0, 1):
        return 1
    elif n < 0:
        return 0 #should be NAN!
    else:
        return n * c_double_factorial(n - 2)   

def double_factorial(double n):
    cdef double val
    with nogil:
        val = c_double_factorial(n)
    return val

##############################################################################
# Determine matrix size needed for direct coupling between 2 particles
############################################################################## 
cdef long c_blocksize(long l_max, long m_max) nogil:
    return c_multi_to_single_index(1, l_max, m_max,l_max,m_max) + 1

def blocksize(long l_max, long m_max):
    cdef long val
    with nogil:
        val = c_blocksize(l_max, m_max)
    return val

##############################################################################
# Map tau, l, m to a single linear index - Speed ~0.3x when called from Python wrap
############################################################################## 
cdef long c_multi_to_single_index(long tau,
                                 long l,
                                 long m,
                                 long l_max,
                                 long m_max) nogil:
    cdef long tau_blocksize
    cdef long n
    
    tau_blocksize = m_max * (m_max + 2) + (l_max - m_max) * (2 * m_max + 1)
    n = tau * tau_blocksize
    if (l - 1) <= m_max:
        n += (l - 1) * (l - 1 + 2)
    else:
        n += m_max * (m_max + 2) + (l - 1 - m_max) * (2 * m_max + 1)
    n += m + min(l, m_max)
    return n

def multi_to_single_index(long tau, long l, long m, long l_max, long m_max):
    cdef long val
    with nogil:
        val = c_multi_to_single_index(tau, l, m, l_max, m_max)
    return val

##############################################################################
#  Spherical Hankel function - Speed ~1x when called from Python wrap
############################################################################## 
cdef double complex c_spherical_hankel(long n, double complex x) nogil:
    cdef double complex spherej
    cdef double complex spherey
    cdef double complex sphereh
    
    spherej = scipy.special.cython_special.spherical_jn(n, x)
    if x == 0:
        spherey = -1e24
    else:
        spherey = scipy.special.cython_special.spherical_yn(n, x)
        
    #if spherej < 1e-15:
    #    spherej = 0
    #if spherey < 1e-15:
    #    spherey = 0        
    
    sphereh = spherej + 1j * spherey  
    return sphereh  

def spherical_hankel(long n,double complex x):
    cdef double complex val
    with nogil:
        val = c_spherical_hankel(n,x)
    return val

##############################################################################
#  Spherical Bessel function for translation matrix
############################################################################## 
cdef double complex c_spherical_bessel(long n, double complex x) nogil:
    return scipy.special.cython_special.spherical_jn(n, x)


def spherical_bessel(long n,double complex x):
    cdef double complex val
    with nogil:
        val = c_spherical_bessel(n,x)
    return val

##############################################################################
#  2D Legendre function for direct coupling block - Speed ~153x when called from Python wrap
############################################################################## 
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:,:] c_legendre_normalized_2D(long lmax) nogil:
    
    cdef double[:,:] plm
    cdef long l
    cdef long m

    with gil:
        plm = np.zeros((lmax+1,lmax+1),dtype = np.float64)    
    
    plm[0][0] = sqrt(2)/2
    plm[1][0] = 0

    for l in range(1, lmax):
        plm[l + 1][0] = (-l / (l + 1)) * sqrt((2 * l + 3) / (2 * l - 1)) * plm[l-1][0]

    for m in range(1, lmax + 1):
        plm[m][m] = sqrt((2 * m + 1) / (2 * c_factorial(2 * m))) * c_double_factorial(2 * m - 1) 

        for l in range(m, lmax):
            plm[l + 1][m] = -sqrt(((2 * l + 3) * (l - m) * (l + m)) / ((2 * l - 1) * (l + 1 - m) * (l + 1 + m))) * plm[l - 1][m]

    return plm  

def legendre_normalized_2D(long lmax):
    cdef double[:,:] val
    with nogil:
         val = c_legendre_normalized_2D(lmax)
    return np.asarray(val, dtype = np.float64)

##############################################################################
#  3D Legendre function for direct coupling block
############################################################################## 
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:,:] c_legendre_normalized_3D(double ct, double st, long lmax) nogil:

    cdef double[:,:] plm
    cdef long l
    cdef long m

    with gil:
        plm = np.zeros((lmax+1,lmax+1),dtype = np.float64)    
    

    plm[0][0] = sqrt(2)/2
    plm[1][0] = sqrt(3/2) * ct

    for l in range(1, lmax):
        plm[l + 1][0] = (1 / (l + 1) * sqrt((2 * l + 1) * (2 * l + 3)) * ct * plm[l][0] -
                         l / (l + 1) * sqrt((2 * l + 3) / (2 * l - 1)) * plm[l-1][0])
    for m in range(1, lmax + 1):
        plm[m][m] = sqrt((2 * m + 1) / (2 * c_factorial(2 * m))) * c_double_factorial(2 * m - 1) * st**m
        for l in range(m, lmax):
            plm[l + 1][m] = (sqrt((2 * l + 1) * (2 * l + 3) / ((l + 1 - m) * (l + 1 + m))) * ct * plm[l][m] -
                             sqrt((2 * l + 3) * (l - m) * (l + m) / ((2 * l - 1) * (l + 1 - m) * (l + 1 + m))) *
                             plm[l - 1][m])
    return plm  

def legendre_normalized_3D(double ct, double st, long lmax):
    cdef double[:,:] val
    with nogil:
        val = c_legendre_normalized_3D(ct, st, lmax)
    return np.asarray(val)


##############################################################################
#  AB5 coefficients - Speed ~2x when called from Python wrap
############################################################################## 
from numba.core.typing import cffi_utils
from pywigxjpf_ffi import ffi, lib
import pywigxjpf_ffi
cffi_utils.register_module(pywigxjpf_ffi)
nb_wig3jj = pywigxjpf_ffi.lib.wig3jj

lib.wig_table_init(100,9)
lib.wig_temp_init(100)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef (double complex, double complex) c_ab5_coefficients(long l1,
                        long m1,
                        long l2, 
                        long m2, 
                        long p) nogil:
    
    cdef double complex jfac = 0+1j*0
    cdef double complex fac1 = 0+1j*0
    cdef double complex fac2a = 0+1j*0
    cdef double complex fac2b = 0+1j*0
    cdef double complex a ,b
    cdef double wig1, wig2a, wig2b,

    jfac = cpow(1.0j,(abs(m1 - m2) - abs(m1) - abs(m2) + l2 - l1 + p)) * (cpow((-1), abs((m1 - m2))))
    fac1 = sqrt(((2 * l1 + 1) * (2 * l2 + 1)) / (2 * l1 * (l1 + 1) * l2 * (l2 + 1)))
    fac2a = (l1 * (l1 + 1) + l2 * (l2 + 1) - p * (p + 1)) * sqrt(2 * p + 1)
    fac2b = sqrt((l1 + l2 + 1 + p) * (l1 + l2 + 1 - p) * (p + l1 - l2) * (p - l1 + l2) * (2 * p + 1))
    
    # Note that arguments are in two_j = 2*j
    with gil:
        wig1 = nb_wig3jj(2*l1, 2*l2, 2*p, 2*m1, -m2*2, -(m1 - m2)*2)
        wig2a = nb_wig3jj(2*l1, 2*l2, 2*p, 0, 0, 0)
        wig2b = nb_wig3jj(2*l1, 2*l2, 2*(p - 1), 0, 0, 0)

    # Instead of 0 the result gives numbers <e-15. Force these to 0 so it matches python.
    # The problem is probably in the jfac 
    a = creal(jfac * fac1 * fac2a * wig1 * wig2a)+1j*0  # As are always real.
    b = 0+1j*cimag(jfac * fac1 * fac2b * wig1 * wig2b)  # Bs are always imag.
    
    return a, b

def cab5_coefficients(long l1, long m1, long l2, long m2, long p):
    cdef double complex a
    cdef double complex b
    with nogil:
        a,b = c_ab5_coefficients(l1, m1, l2, m2, p)
    return complex(a), complex(b)


'''
##############################################################################
#****************************************************************************#
#****************************************************************************#
#**************** 3D direct coupling block and translation ******************#
#********************* Includes multithread options *************************#
#****************************************************************************#
##############################################################################

This section details the direct coupling block and translation between 2 
particles that can be anywhere in the same layer.

Since the direct coupling block is actually just a translation using a different
spherical bessel funciton, use a "kind" flag to chose the bessel function. Then, 
implement python wrappers to call this flag appropriately based on implementation.
This prevents coding error and makes the code more readable.

Options for multithreading are added. 
In case these options cause errors for different implementations/configurations, 
multithreading compatible functions are given a different name so a fall 
back option is accessable.
'''

##############################################################################
# Translate block between 2 particles in 3D
##############################################################################
#@cython.cdivision(True) <- DO NOT ENABLE THIS. IT WILL DO INTERGER DIVISION!
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double complex [:,:] c_svwf_translate_3D(double complex k,
                                                  double dx,
                                                  double dy,
                                                  double dz,
                                                  long lmax1,
                                                  long lmax2,
                                                  long mmax1,
                                                  long mmax2, 
                                                  long kind) nogil:

    # Type things and calculate simple 1 line formulas
    # Note: Math functions are optimized by standard c libraries
    cdef long blocksize1 = c_blocksize(lmax1, mmax1)
    cdef long blocksize2 = c_blocksize(lmax2, mmax2)
    cdef double complex d = csqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))
    cdef double d_real = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))
    cdef double cos_theta = dz / d_real
    cdef double sin_theta = sqrt(dx**2 + dy**2) / d_real
    cdef double complex kd = k*d
    cdef double[:,:] legendre
    cdef double complex[:,::1] w 
    cdef double complex [:] sph
    cdef int[:] m1_array, m2_array
    cdef long im1, im2, l1, l2, ld, n1, n2, tau1, tau2, i, m1, m2
    cdef double complex A, B, a5, b5, eimph
    cdef double phi 
    
    # All numpy initialize functions must be done with the Gil. 
    # Other methods such as malloc or pointers can be used for pure C
    # Though, I see no reason. 
    with gil:
        phi = np.arctan2(dy, dx)
        w = np.zeros((blocksize1, blocksize2), dtype=np.complex128, order = 'C')   
        # If the particle is itself then coupling matrix = 0 matrix.
        if d == 0:
            return w
        sph = np.zeros((lmax1+lmax2+1), dtype=np.complex128)
        m1_array = np.arange(-mmax1, mmax1+1, 1, dtype = np.intc)
        m2_array = np.arange(-mmax2, mmax2+1, 1, dtype = np.intc)

    if kind == 0: # 0 = incoming wave, 1 = outgoing wave. 
    #NOTE: an long flag is used instead of a string because Cython 
    # has troublesome conventions about datatyping strings from python v2 to v3
    # See more on this at https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_hankel(n, kd) 
    else:
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_bessel(n, kd)       
        
    legendre = c_legendre_normalized_3D(cos_theta, sin_theta, lmax1 + lmax2) 

    
    for im1 in range(0, 2*mmax1 + 1):
        m1 = m1_array[im1]
        for im2 in range(0, 2*mmax2 + 1):
            m2 = m2_array[im2]
            eimph = cexp(1j * (m2 - m1) * phi)
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    A = 0
                    B = 0
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1):  
                        a5, b5 = c_ab5_coefficients(l2, m2, l1, m1, ld)
                        A += a5 * sph[ld] * legendre[ld][abs(m1 - m2)]
                        B += b5 * sph[ld] * legendre[ld][abs(m1 - m2)]
                    A, B = eimph * A, eimph * B
                    for tau1 in range(2):
                        n1 = c_multi_to_single_index(tau1, l1, m1, lmax1, mmax1)
                        for tau2 in range(2):
                            n2 = c_multi_to_single_index(tau2, l2, m2, lmax2, mmax2)
                            if tau1 == tau2:
                                w[n1, n2] = A
                            else:
                                w[n1, n2] = B

    return w

##############################################################################
# Translate block between 2 particles in 3D with multithreading enabled
##############################################################################
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double complex [:,:] c_svwf_translate_3D_threaded(double complex k,
                                                  double dx,
                                                  double dy,
                                                  double dz,
                                                  long lmax1,
                                                  long lmax2,
                                                  long mmax1,
                                                  long mmax2, 
                                                  long kind) nogil:

    # Type things and calculate simple 1 line formulas
    # Note: Math functions are optimized by standard c libraries
    cdef long blocksize1 = c_blocksize(lmax1, mmax1)
    cdef long blocksize2 = c_blocksize(lmax2, mmax2)
    cdef double complex d = csqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))
    cdef double d_real = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))
    cdef double cos_theta = dz / d_real
    cdef double sin_theta = sqrt(dx**2 + dy**2) / d_real
    cdef double complex kd = k*d
    cdef double[:,:] legendre
    cdef double complex[:,::1] w 
    cdef double complex [:] sph
    cdef int[:] m1_array, m2_array
    cdef long im1, im2, l1, l2, ld, n1, n2, tau1, tau2, i, m1, m2
    cdef double complex A, B, a5, b5, eimph
    cdef double phi 
    
    # All numpy initialize functions must be done with the Gil. 
    # Other methods such as malloc or pointers can be used for pure C
    # Though, I see no reason. 
    with gil:
        phi = np.arctan2(dy, dx)
        w = np.zeros((blocksize1, blocksize2), dtype=np.complex128, order = 'C')   
        # If the particle is itself then coupling matrix = 0 matrix.
        if d == 0:
            return w
        sph = np.zeros((lmax1+lmax2+1), dtype=np.complex128)
        m1_array = np.arange(-mmax1, mmax1+1, 1, dtype = np.intc)
        m2_array = np.arange(-mmax2, mmax2+1, 1, dtype = np.intc)

    if kind == 0: # 0 = incoming wave, 1 = outgoing wave. 
    #NOTE: an long flag is used instead of a string because Cython 
    # has troublesome conventions about datatyping strings from python v2 to v3
    # See more on this at https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_hankel(n, kd) 
    else:
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_bessel(n, kd)       
        
    legendre = c_legendre_normalized_3D(cos_theta, sin_theta, lmax1 + lmax2) 

    
    for im1 in prange(0, 2*mmax1 + 1): #<- only change to thread since gil already released!
        m1 = m1_array[im1]
        for im2 in range(0, 2*mmax2 + 1):
            m2 = m2_array[im2]
            eimph = cexp(1j * (m2 - m1) * phi)
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    A = 0
                    B = 0
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1):  
                        a5, b5 = c_ab5_coefficients(l2, m2, l1, m1, ld)
                        A += a5 * sph[ld] * legendre[ld][abs(m1 - m2)]
                        B += b5 * sph[ld] * legendre[ld][abs(m1 - m2)]
                    A, B = eimph * A, eimph * B
                    for tau1 in range(2):
                        n1 = c_multi_to_single_index(tau1, l1, m1, lmax1, mmax1)
                        for tau2 in range(2):
                            n2 = c_multi_to_single_index(tau2, l2, m2, lmax2, mmax2)
                            if tau1 == tau2:
                                w[n1, n2] = A
                            else:
                                w[n1, n2] = B

    return w


# Have option to thread build into the translate function for easier reuse.
def svwf_translate_3D(double complex k, double dx, double dy, double dz,
                      long lmax1, long lmax2, long mmax1, long mmax2,
                      long kind, threaded = False):
    cdef double complex [:,:] val
    
    if threaded == True:
        with nogil:
            val = c_svwf_translate_3D_threaded(k, dx, dy, dz, lmax1, lmax2, mmax1, mmax2, kind)
    else:
        with nogil:
            val = c_svwf_translate_3D(k, dx, dy, dz, lmax1, lmax2, mmax1, mmax2, kind)
    return np.asarray(val, dtype = np.complex128)

# Special case of translate. Make wrapper for readability in smuthi.
def direct_coupling_block_3D(double complex k, double dx, double dy, double dz,
                              long lmax1, long lmax2, long mmax1, long mmax2, threaded = False):
    cdef kind = 0
    return svwf_translate_3D(k, dx, dy,  dz,
                             lmax1,lmax2,mmax1,mmax2,
                             kind, threaded)


'''
##############################################################################
#****************************************************************************#
#****************************************************************************#
#********** 3D direct coupling block and translation hash   *****************#
#********************* Includes multithread options *************************#
#****************************************************************************#
##############################################################################

The same AB5 coefficients are often shared between particle pairs. In Python
these are memoized. In Cython, we memoize a hash table of particle order 
pairs (e.g., 2->2 , 3->2, 2->3, ...). This table is then used for fast C loops.

This section details the direct coupling block and translation between 2 
particles that can be anywhere in the same layer. Options for multithreading 
are added. In case these options cause errors for different implementations, 
multithreading compatible functions are given a different name so a fall 
back option is accessable.
'''

##############################################################################
#  Make  hash table of ab5 for a particle order pair (WITH GIL)
##############################################################################
# """
# This hash table is particularly good when all particles have the same order. 
# """
@cython.boundscheck(False)
@cython.wraparound(False)
cdef c_hard_ab5_coefficient_hash_table(long lmax1, 
                                 long lmax2, 
                                 long mmax1, 
                                 long mmax2):   #<-- Note: This function holds the gil!
    
    cdef long ab5counter = 0
    cdef long m1, m2, l1, l2, ld
    for m1 in range(-mmax1, mmax1 + 1):
        for m2 in range(-mmax2, mmax2 + 1):
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1): 
                        ab5counter += 1
    
   
    cdef double complex[:] a5_array, b5_array
    cdef double complex a5, b5
    
    a5_array = np.ones((ab5counter,), dtype=np.complex128)
    b5_array = np.ones((ab5counter,), dtype=np.complex128)
   
    ab5counter = 0
    for m1 in range(-mmax1, mmax1 + 1):
        for m2 in range(-mmax2, mmax2 + 1):
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1): 
                        a5, b5 = c_ab5_coefficients(l2, m2, l1, m1, ld)
                        a5_array[ab5counter] = a5
                        b5_array[ab5counter] = b5
                        ab5counter+=1
                        
    return a5_array, b5_array

def _ab5_coefficient_hash_table(long lmax1, long lmax2, long mmax1, long mmax2):
    cdef double complex[:] a5_array, b5_array
    a5_array, b5_array = c_hard_ab5_coefficient_hash_table(lmax1, lmax2, mmax1, mmax2)
    return np.asarray(a5_array, dtype = np.complex128),  np.asarray(b5_array, dtype = np.complex128)



##############################################################################
# Translation 3D with  Hash
##############################################################################
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double complex [:,:]  c_svwf_translate_3D_from_hash_table(long blocksize1,long blocksize2,
                            double complex[:,:] w,
                            double complex[:] sph,
                            double complex k,
                            double dx, double dy, double dz,
                            long lmax1, long lmax2, long mmax1, long mmax2,
                            double complex[:] a5_array, 
                            double complex[:] b5_array, 
                            long kind) nogil:



    # calculate distance vector
    cdef double complex d = csqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))
    cdef double d_real = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))
    cdef double cos_theta = dz / d_real
    cdef double sin_theta = sqrt(dx**2 + dy**2) / d_real
    cdef double complex kd = k*d
    cdef double[:,:] legendre
    cdef int[:] m1_array, m2_array
    
    # calculate coupling matrix
    cdef long im1, im2, l1, l2, ld, n1, n2, tau1, tau2, i, m1, m2
    cdef double complex A, B, a5, b5, eimph
    cdef double phi = atan2(dy, dx)
    
    with gil:
        # If the particle is itself then coupling matrix = 0 matrix.
        m1_array = np.arange(-mmax1, mmax1+1, 1, dtype = np.intc)
        m2_array = np.arange(-mmax2, mmax2+1, 1, dtype = np.intc)
        
        
    if d == 0:
        return w    

    if kind == 0: # 0 = incoming wave, 1 = outgoing wave. 
    #NOTE: an long flag is used instead of a string because Cython 
    # has troublesome conventions about datatyping strings from python v2 to v3
    # See more on this at https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_hankel(n, kd) 
    else:
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_bessel(n, kd)      
    legendre = c_legendre_normalized_3D(cos_theta, sin_theta, lmax1 + lmax2) 

    cdef long ab5counter = 0
    for im1 in range(0, 2*mmax1 + 1): # <- Only necessary change for threading since gil is already released 
        m1 = m1_array[im1]
        for im2 in range(0, 2*mmax2 + 1):
            m2 = m2_array[im2]
            eimph = cexp(1j * (m2 - m1) * phi)
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    A = 0
                    B = 0
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1):  
                        A += a5_array[ab5counter] * sph[ld] * legendre[ld][abs(m1 - m2)]
                        B += b5_array[ab5counter] * sph[ld] * legendre[ld][abs(m1 - m2)]
                        ab5counter+=1
                    A, B = eimph * A, eimph * B
                    for tau1 in range(2):
                        n1 = c_multi_to_single_index(tau1, l1, m1, lmax1, mmax1)
                        for tau2 in range(2):
                            n2 = c_multi_to_single_index(tau2, l2, m2, lmax2, mmax2)
                            if tau1 == tau2:
                                w[n1, n2] = A
                            else:
                                w[n1, n2] = B

    return w


##############################################################################
# Translation 3D with  Hash and Multithreading
##############################################################################
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double complex [:,:]  c_svwf_translate_3D_from_hash_table_threaded(long blocksize1,long blocksize2,
                            double complex[:,:] w,
                            double complex[:] sph,
                            double complex k,
                            double dx, double dy, double dz,
                            long lmax1, long lmax2, long mmax1, long mmax2,
                            double complex[:] a5_array, 
                            double complex[:] b5_array,
                            long kind) nogil:



    # calculate distance vector
    cdef double complex d = csqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))
    cdef double d_real = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))
    cdef double cos_theta = dz / d_real
    cdef double sin_theta = sqrt(dx**2 + dy**2) / d_real
    cdef double complex kd = k*d
    cdef double[:,:] legendre
    cdef int[:] m1_array, m2_array
    
    # calculate coupling matrix
    cdef long im1, im2, l1, l2, ld, n1, n2, tau1, tau2, i, m1, m2
    cdef double complex A, B, a5, b5, eimph
    cdef double phi = atan2(dy, dx)
    
    with gil:
        # If the particle is itself then coupling matrix = 0 matrix.
        m1_array = np.arange(-mmax1, mmax1+1, 1, dtype = np.intc)
        m2_array = np.arange(-mmax2, mmax2+1, 1, dtype = np.intc)
        
        
    if d == 0:
        return w    

    if kind == 0: # 0 = incoming wave, 1 = outgoing wave. 
    #NOTE: an long flag is used instead of a string because Cython 
    # has troublesome conventions about datatyping strings from python v2 to v3
    # See more on this at https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_hankel(n, kd) 
    else:
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_bessel(n, kd)      
    legendre = c_legendre_normalized_3D(cos_theta, sin_theta, lmax1 + lmax2) 

    cdef long ab5counter = 0
    for im1 in prange(0, 2*mmax1 + 1): # <- Only necessary change for threading since gil is already released 
        m1 = m1_array[im1]
        for im2 in range(0, 2*mmax2 + 1):
            m2 = m2_array[im2]
            eimph = cexp(1j * (m2 - m1) * phi)
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    A = 0
                    B = 0
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1):  
                        A += a5_array[ab5counter] * sph[ld] * legendre[ld][abs(m1 - m2)]
                        B += b5_array[ab5counter] * sph[ld] * legendre[ld][abs(m1 - m2)]
                        ab5counter+=1
                    A, B = eimph * A, eimph * B
                    for tau1 in range(2):
                        n1 = c_multi_to_single_index(tau1, l1, m1, lmax1, mmax1)
                        for tau2 in range(2):
                            n2 = c_multi_to_single_index(tau2, l2, m2, lmax2, mmax2)
                            if tau1 == tau2:
                                w[n1, n2] = A
                            else:
                                w[n1, n2] = B

    return w


# Have option to thread build into the translate function for easier reuse.
def svwf_translate_3D_from_hash_table(long blocksize1,long blocksize2,
                            double complex[:,:] w,
                            double complex[:] sph,
                            double complex k,
                            double dx, double dy, double dz,
                            long lmax1, long lmax2, long mmax1, long mmax2,
                            double complex[:] a5_array, 
                            double complex[:] b5_array,
                            long kind, threaded = False):
    
    cdef double complex [:,:] val
    if threaded == True:
        with nogil:
            val = c_svwf_translate_3D_from_hash_table_threaded(blocksize1,blocksize2,
                                            w,sph,k, dx, dy,  dz,
                                            lmax1,lmax2,mmax1,mmax2,
                                            a5_array, b5_array,
                                            kind)
    else:
        with nogil:
            val = c_svwf_translate_3D_from_hash_table(blocksize1,blocksize2,
                                            w,sph,k, dx, dy,  dz,
                                            lmax1,lmax2,mmax1,mmax2,
                                            a5_array, b5_array,
                                            kind)
            
    return np.asarray(val, dtype = np.complex128)

# Special case of translate. Make wrapper for readability in smuthi.
def direct_coupling_block_3D_from_hash_table(long blocksize1,long blocksize2,
                                double complex[:,:] w,
                                double complex[:] sph,
                                double complex k,
                                double dx, double dy, double dz,
                                long lmax1, long lmax2, long mmax1, long mmax2,
                                double complex[:] a5_array, 
                                double complex[:] b5_array,                             
                                threaded = False):
    cdef kind = 0
    return svwf_translate_3D_from_hash_table(blocksize1,blocksize2,
                                    w,sph,k, dx, dy,  dz,
                                    lmax1,lmax2,mmax1,mmax2,
                                    a5_array, b5_array,
                                    kind, threaded)


'''
##############################################################################
#****************************************************************************#
#****************************************************************************#
#**************** 2D direct coupling block and translation ******************#
#********************* Includes multithread options *************************#
#****************************************************************************#
##############################################################################

This section details the direct coupling block and translation between 2 
particles that can be anywhere in the same x-y plane. Options for multithreading 
are added. In case these options cause errors for different implementations, 
multithreading compatible functions are given a different name so a fall 
back option is accessable.
'''


##############################################################################
# Translate block between 2 particles in 2D (Legendre = 90 deg)
##############################################################################


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double complex [:,:] c_svwf_translate_2D(double complex k,
                                                  double dx,
                                                  double dy,
                                                  double dz,
                                                  long lmax1,
                                                  long lmax2,
                                                  long mmax1,
                                                  long mmax2,
                                                  long kind) nogil:

    # Get index
    cdef long blocksize1 = c_blocksize(lmax1, mmax1)
    cdef long blocksize2 = c_blocksize(lmax2, mmax2)
    
    # calculate coupling matrix
    cdef long im1, im2, l1, l2, ld, n1, n2, tau1, tau2, i, m1, m2
    cdef double complex A, B, a5, b5, eimph, d, phi, kd
    cdef double complex[:,::1] w
    cdef double complex [:] sph
    cdef int[:] m1_array, m2_array 
    cdef double[:,:] legendre
    
    with gil:
        # calculate distance vector
        d = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))
        phi = np.arctan2(dy, dx)

        #initialize coupling matrix
        w = np.zeros((blocksize1, blocksize2), dtype=np.complex128, order = 'C')
        # If the particle is itself then coupling matrix = 0 matrix.
        if d == 0:
            return w
        sph = np.zeros((lmax1+lmax2+1), dtype=np.complex128)
        m1_array = np.arange(-mmax1, mmax1+1, 1, dtype = np.intc)
        m2_array = np.arange(-mmax2, mmax2+1, 1, dtype = np.intc)



    # spherical functions  
    kd = k*d
    if kind == 0: # 0 = incoming wave, 1 = outgoing wave. 
    #NOTE: an long flag is used instead of a string because Cython 
    # has troublesome conventions about datatyping strings from python v2 to v3
    # See more on this at https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_hankel(n, kd) 
    else:
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_bessel(n, kd)  
    legendre = c_legendre_normalized_2D(lmax1 + lmax2) 

    
    for im1 in range(0, 2*mmax1 + 1):
        m1 = m1_array[im1]
        for im2 in range(0, 2*mmax2 + 1):
            m2 = m2_array[im2]
            eimph = cexp(1j * (m2 - m1) * phi)
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    A = 0
                    B = 0
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1):  
                        a5, b5 = c_ab5_coefficients(l2, m2, l1, m1, ld)
                        A += a5 * sph[ld] * legendre[ld][abs(m1 - m2)]
                        B += b5 * sph[ld] * legendre[ld][abs(m1 - m2)]
                    A, B = eimph * A, eimph * B
                    for tau1 in range(2):
                        n1 = c_multi_to_single_index(tau1, l1, m1, lmax1, mmax1)
                        for tau2 in range(2):
                            n2 = c_multi_to_single_index(tau2, l2, m2, lmax2, mmax2)
                            if tau1 == tau2:
                                w[n1, n2] = A
                            else:
                                w[n1, n2] = B

    return w



##############################################################################
# Translate block between 2 particles in 2D (Legendre = 90 deg) with threading
##############################################################################
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double complex [:,:] c_svwf_translate_2D_threaded(double complex k,
                                                  double dx,
                                                  double dy,
                                                  double dz,
                                                  long lmax1,
                                                  long lmax2,
                                                  long mmax1,
                                                  long mmax2,
                                                  long kind) nogil:

    # Get index
    cdef long blocksize1 = c_blocksize(lmax1, mmax1)
    cdef long blocksize2 = c_blocksize(lmax2, mmax2)
    
    # calculate coupling matrix
    cdef long im1, im2, l1, l2, ld, n1, n2, tau1, tau2, i, m1, m2
    cdef double complex A, B, a5, b5, eimph, d, phi, kd
    cdef double complex[:,::1] w
    cdef double complex [:] sph
    cdef int[:] m1_array, m2_array 
    cdef double[:,:] legendre
    
    with gil:
        # calculate distance vector
        d = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))
        phi = np.arctan2(dy, dx)

        #initialize coupling matrix
        w = np.zeros((blocksize1, blocksize2), dtype=np.complex128, order = 'C')
        # If the particle is itself then coupling matrix = 0 matrix.
        if d == 0:
            return w
        sph = np.zeros((lmax1+lmax2+1), dtype=np.complex128)
        m1_array = np.arange(-mmax1, mmax1+1, 1, dtype = np.intc)
        m2_array = np.arange(-mmax2, mmax2+1, 1, dtype = np.intc)



    # spherical functions  
    kd = k*d
    if kind == 0: # 0 = incoming wave, 1 = outgoing wave. 
    #NOTE: an long flag is used instead of a string because Cython 
    # has troublesome conventions about datatyping strings from python v2 to v3
    # See more on this at https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_hankel(n, kd) 
    else:
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_bessel(n, kd)  
    legendre = c_legendre_normalized_2D(lmax1 + lmax2) 

    
    for im1 in prange(0, 2*mmax1 + 1):
        m1 = m1_array[im1]
        for im2 in range(0, 2*mmax2 + 1):
            m2 = m2_array[im2]
            eimph = cexp(1j * (m2 - m1) * phi)
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    A = 0
                    B = 0
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1):  
                        a5, b5 = c_ab5_coefficients(l2, m2, l1, m1, ld)
                        A += a5 * sph[ld] * legendre[ld][abs(m1 - m2)]
                        B += b5 * sph[ld] * legendre[ld][abs(m1 - m2)]
                    A, B = eimph * A, eimph * B
                    for tau1 in range(2):
                        n1 = c_multi_to_single_index(tau1, l1, m1, lmax1, mmax1)
                        for tau2 in range(2):
                            n2 = c_multi_to_single_index(tau2, l2, m2, lmax2, mmax2)
                            if tau1 == tau2:
                                w[n1, n2] = A
                            else:
                                w[n1, n2] = B

    return w


# Have option to thread build into the translate function for easier reuse.
def svwf_translate_2D(double complex k, double dx, double dy, double dz,
                      long lmax1, long lmax2, long mmax1, long mmax2,
                      long kind, threaded = False):
    
    cdef double complex [:,:] val

    if threaded == True:
        with nogil:
            val = c_svwf_translate_2D_threaded(k, dx, dy, dz, 
                                              lmax1, lmax2, mmax1, mmax2,
                                              kind)
    else:
        with nogil:
            val = c_svwf_translate_2D(k, dx, dy, dz, 
                                     lmax1, lmax2, mmax1, mmax2,
                                     kind)
            
    return np.asarray(val, dtype = np.complex128)

# Special case of translate. Make wrapper for readability in smuthi.
def direct_coupling_block_2D(double complex k, double dx, double dy, double dz,
                              long lmax1, long lmax2, long mmax1, long mmax2,                          
                              threaded = False):
    cdef kind = 0
    return svwf_translate_2D(k, dx, dy, dz,
                             lmax1, lmax2, mmax1, mmax2, 
                             kind, threaded)

'''
##############################################################################
#****************************************************************************#
#****************************************************************************#
#********** 2D direct coupling block and translation hash   *****************#
#********************* Includes multithread options *************************#
#****************************************************************************#
##############################################################################

The same AB5 coefficients are often shared between particle pairs. In Python
these are memoized. In Cython, we memoize a hash table of particle order 
pairs (e.g., 2->2 , 3->2, 2->3, ...). This table is then used for fast C loops.

This section details the direct coupling block and translation between 2 
particles that can be anywhere in the same layer. Options for multithreading 
are added. In case these options cause errors for different implementations, 
multithreading compatible functions are given a different name so a fall 
back option is accessable.
'''



##############################################################################
#  Make  hash table of ab5 and legendre for a particle pair
##############################################################################
# """
#     This hash table is considered "hard" because it is dependent on the order 
#     of which particle is the emitter and which is the reciever. This hash 
#     table is particularly good when all particles have the same order. 
#     If orders are different, a hash of just a5 and b5 as multidimensional 
#     arrays may be a better option because you dont have to run the hash for 
#     e.g., 2->3 and 3->2. But, no tests have been run to verify there is truly 
#     a speed increase... so, unless this version of a hash becomes the 
#     bottleneck I DO NOT recomend writing and debugging a new more "generalized"
#     hash..
# """
@cython.boundscheck(False)
@cython.wraparound(False)
cdef c_hard_ab5_coefficient_and_legendre_hash_table(long lmax1, 
                                             long lmax2, 
                                             long mmax1, 
                                             long mmax2):
    
    cdef long ab5counter = 0
    cdef long m1, m2, l1, l2, ld
    for m1 in range(-mmax1, mmax1 + 1):
        for m2 in range(-mmax2, mmax2 + 1):
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1): 
                        ab5counter += 1
    
    cdef double complex[:] a5leg_array, b5leg_array
    cdef double complex a5, b5, leg
    cdef double [:,:] legendre

    a5leg_array = np.ones((ab5counter,), dtype=np.complex128)
    b5leg_array = np.ones((ab5counter,), dtype=np.complex128)

    ab5counter = 0
    legendre = c_legendre_normalized_2D(lmax1+lmax2+1)
    for m1 in range(-mmax1, mmax1 + 1):
        for m2 in range(-mmax2, mmax2 + 1):
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1): 
                        a5, b5 = c_ab5_coefficients(l2, m2, l1, m1, ld)
                        leg =  complex(legendre[ld][abs(m1 - m2)])
                        a5leg_array[ab5counter] = leg*a5
                        b5leg_array[ab5counter] = leg*b5
                        ab5counter+=1
                        
    return a5leg_array, b5leg_array

def _ab5_coefficient_and_legendre_hash_table(long lmax1, long lmax2, long mmax1, long mmax2):
    cdef double complex[:] a5leg_array, b5leg_array 
    a5leg_array, b5leg_array = c_hard_ab5_coefficient_and_legendre_hash_table(lmax1, lmax2, mmax1, mmax2)
    return np.asarray(a5leg_array, dtype = np.complex128),  np.asarray(b5leg_array, dtype = np.complex128)

##############################################################################
#  Direct coupling block using  hash table
##############################################################################
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double complex[:,:] c_svwf_translate_2D_from_hash_table(long blocksize1,
                                                            long blocksize2,
                                                            double complex[:,:] w,
                                                            double complex[:] sph,
                                                            double complex k,
                                                            double dx,
                                                            double dy,
                                                            double dz,
                                                            long lmax1,
                                                            long lmax2,
                                                            long mmax1,
                                                            long mmax2,
                                                            double complex[:] a5leg_hash_table,
                                                            double complex[:] b5leg_hash_table,
                                                            long kind) nogil:      

    cdef double phi = atan2(dy, dx)
    cdef double complex d = csqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))

    # If the particle is itself then coupling matrix = 0 matrix.
    if d == 0:
        return w
    
    # spherical functions
    cdef double complex kd = k*d
    
    if kind == 0: # 0 = incoming wave, 1 = outgoing wave. 
    #NOTE: an long flag is used instead of a string because Cython 
    # has troublesome conventions about datatyping strings from python v2 to v3
    # See more on this at https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_hankel(n, kd) 
    else:
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_bessel(n, kd)  
  

    # calculate coupling matrix
    cdef long m1, m2, l1, l2, ld, n1, n2, tau1, tau2, idx_ab5leg, idx_eim 
    cdef double complex A, B, eimph
                                                        
    # the particle coupling operator is the transpose of the SVWF translation operator
    # therefore, (l1,m1) and (l2,m2) are interchanged:
    idx_ab5leg = 0       
    for m1 in range(-mmax1, mmax1 + 1):
        for m2 in range(-mmax2, mmax2 + 1):
            eimph = cexp(1j * (m2 - m1) * phi) 
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    A = 0
                    B = 0
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1):  # if ld<abs(m1-m2) then P=0
                       
                        A += a5leg_hash_table[idx_ab5leg]*sph[ld] 
                        B += b5leg_hash_table[idx_ab5leg]*sph[ld] 
                        idx_ab5leg += 1
                                                        
                    A, B = eimph * A, eimph * B
                        
                    for tau1 in range(2):
                        n1 = c_multi_to_single_index(tau1, l1, m1, lmax1, mmax1)
                        for tau2 in range(2):
                            n2 = c_multi_to_single_index(tau2, l2, m2, lmax2, mmax2)
                            if tau1 == tau2:
                                w[n1, n2] = A
                            else:
                                w[n1, n2] = B

    return w
##############################################################################
#  Direct coupling block using  hash table
##############################################################################
@cython.boundscheck(False)
@cython.wraparound(False)
cdef double complex[:,:] c_svwf_translate_2D_from_hash_table_threaded(long blocksize1,
                                                            long blocksize2,
                                                            double complex[:,:] w,
                                                            double complex[:] sph,
                                                            double complex k,
                                                            double dx,
                                                            double dy,
                                                            double dz,
                                                            long lmax1,
                                                            long lmax2,
                                                            long mmax1,
                                                            long mmax2,
                                                            double complex[:] a5leg_hash_table,
                                                            double complex[:] b5leg_hash_table,
                                                            long kind) nogil:      

                                 


    cdef double phi = atan2(dy, dx)
    cdef double complex d = csqrt(pow(dx,2) + pow(dy,2) + pow(dz,2))

    # If the particle is itself then coupling matrix = 0 matrix.
    if d == 0:
        return w
    
    # spherical functions
    cdef double complex kd = k*d
    
    if kind == 0: # 0 = incoming wave, 1 = outgoing wave. 
    #NOTE: an long flag is used instead of a string because Cython 
    # has troublesome conventions about datatyping strings from python v2 to v3
    # See more on this at https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_hankel(n, kd) 
    else:
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = c_spherical_bessel(n, kd)  
  

    # calculate coupling matrix
    cdef long m1, m2, l1, l2, ld, n1, n2, tau1, tau2, idx_ab5leg, idx_eim 
    cdef double complex A, B, eimph
                                                        
    # the particle coupling operator is the transpose of the SVWF translation operator
    # therefore, (l1,m1) and (l2,m2) are interchanged:
    idx_ab5leg = 0       
    for m1 in prange(-mmax1, mmax1 + 1):
        for m2 in range(-mmax2, mmax2 + 1):
            eimph = cexp(1j * (m2 - m1) * phi) 
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    A = 0
                    B = 0
                    for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1):  # if ld<abs(m1-m2) then P=0
                       
                        A += a5leg_hash_table[idx_ab5leg]*sph[ld] 
                        B += b5leg_hash_table[idx_ab5leg]*sph[ld] 
                        idx_ab5leg += 1
                                                        
                    A, B = eimph * A, eimph * B
                        
                    for tau1 in range(2):
                        n1 = c_multi_to_single_index(tau1, l1, m1, lmax1, mmax1)
                        for tau2 in range(2):
                            n2 = c_multi_to_single_index(tau2, l2, m2, lmax2, mmax2)
                            if tau1 == tau2:
                                w[n1, n2] = A
                            else:
                                w[n1, n2] = B

    return w


# Have option to thread build into the translate function for easier reuse.
def svwf_translate_2D_from_hash_table(long blocksize1,long blocksize2,
                            double complex[:,:] w,
                            double complex[:] sph,
                            double complex k,
                            double dx, double dy, double dz,
                            long lmax1, long lmax2, long mmax1, long mmax2,
                            double complex[:] a5leg_array, 
                            double complex[:] b5leg_array,
                            long kind, threaded = False):    

    cdef double complex [:,:] val
    if threaded == True:
        with nogil:
            val = c_svwf_translate_2D_from_hash_table_threaded(blocksize1,blocksize2,
                                                    w,sph,k,dx, dy, dz,
                                                    lmax1, lmax2, mmax1, mmax2,
                                                    a5leg_array,b5leg_array,
                                                    kind)
    else:
        with nogil:
            val = c_svwf_translate_2D_from_hash_table(blocksize1,blocksize2,
                                            w,sph,k,dx, dy, dz,
                                            lmax1, lmax2, mmax1, mmax2,
                                            a5leg_array,b5leg_array,
                                            kind)
            
    return np.asarray(val, dtype = np.complex128)

# Special case of translate. Make wrapper for readability in smuthi.
def direct_coupling_block_2D_from_hash_table(long blocksize1,long blocksize2,
                            double complex[:,:] w,
                            double complex[:] sph,
                            double complex k,
                            double dx, double dy, double dz,
                            long lmax1, long lmax2, long mmax1, long mmax2,
                            double complex[:] a5leg_array, 
                            double complex[:] b5leg_array,
                            threaded = False):
    cdef long kind = 0
    return svwf_translate_2D_from_hash_table(blocksize1,blocksize2,
                                    w,sph,k,dx, dy, dz,
                                    lmax1, lmax2, mmax1, mmax2,
                                    a5leg_array,b5leg_array,
                                    kind, threaded)



'''
##############################################################################
#****************************************************************************#
#****************************************************************************#
#*******************  Class for managing translations  **********************#
#********************* Includes multithread options *************************#
#****************************************************************************#
##############################################################################

This section defines a class that can more easily abstract all the translation 
options defined above. This could be used as a more readable method to implement
the different options. 

Note: This class is based on the idea that the hash options are always faster
compared to the direct calculation. Also, the initial cost in setup and memory 
of the hash is minimal. Therefore, it should always be used as the default 
option. 

Futhermore, the class is a python class so that getter and setter methods do not
need to be defined in order to pickle. 
'''


# class SWE_DirectTranslation():
    
#     def __init__(self, lmax1, mmax1, lmax2, mmax2):
#         self.lmax1 = int(lmax1)
#         self.lmax2 = int(lmax2)
#         self.mmax1 = int(mmax1)
#         self.mmax2 = int(mmax2)
#         self.a5leg_array = None
#         self.b5leg_array = None
#         self.a5_array = None
#         self.b5_array = None  
        
#         return 
        
#     def is_same(self, lmax1, mmax1, lmax2, mmax2):
#         return (self.lmax1 == lmax1 and self.lmax2 == lmax2 and self.mmax1 == mmax1 and self.mmax2 == mmax2)
        
#     def direct_coupling_block(self, k, dx, dy, dz, threaded = False):
#         harmonic_kind = 'outgoing to regular'
#         return self.translate(self, k, dx, dy, dz,
#                   harmonic_kind, threaded)
        
#     def translate(self, k, dx, dy, dz,
#                   harmonic_kind, threaded = False):
        
#         if harmonic_kind == 'outgoing to regular':
#             kind = int(0)
#         elif harmonic_kind == 'outgoing to outgoing':
#             kind = int(1)
#         else:
#             raise Exception(f'''The svwf translation matrix transforms a particle into a basis
#         of either "incoming" or "outgoing" spherical vector harmonics.
#         \n You did not properly specify if the harmonic_kind should be "incoming" or "outgoing"
#         \n Instead you specified harmonic_kind = {harmonic_kind}''') 
        
#         if dz == 0:
#             return self._translate_2d(k, dx, dy, dz, kind, threaded)
#         else: 
#             return self._translate_3d(k, dx, dy, dz, kind, threaded)
        

#     def __translate_2d(self, k, dx, dy, dz, kind, threaded):
        
#         if (self.a5leg_array == None or self.b5leg_array == None):
#             self.a5leg_array, self.b5leg_array = c_hard_ab5_coefficient_and_legendre_hash_table(self.lmax1, 
#                                                                                             self.lmax2, 
#                                                                                             self.mmax1, 
#                                                                                             self.mmax2)
#         return csvwf_translate_2D_from_hash_table(complex(k),float(dx),float(dy),float(dz),
#                          self.lmax1,self.lmax2,self.mmax1,self.mmax2,
#                          self.a5leg_hash_table, self.b5leg_hash_table,
#                          kind, threaded)  
    
#     def __translate_3d(self, k, dx, dy, dz, kind, threaded):
        
#         if (self.a5_array == None or self.b5_array == None):
#             self.a5_array, self.b5_array = chard_ab5_coefficient_hash_table(self.lmax1, 
#                                                                         self.lmax2, 
#                                                                         self.mmax1, 
#                                                                         self.mmax2)
            
#         return csvwf_translate_3D_from_hash_table(complex(k),float(dx),float(dy),float(dz),
#                          self.lmax1,self.lmax2,self.mmax1,self.mmax2,
#                          self.a5_hash_table, self.b5_hash_table,
#                          kind, threaded)  
        





















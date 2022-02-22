"""This module contains helper functions to evaluate the coupling between a 
single particle and a periodic particle arrangement."""

import smuthi.fields as flds
from numba import njit
import numpy as np


def a5b5_lookup(lmax1, lmax2, mmax1, mmax2):
    """ Prepare lookup table of smuthi.vector_wave_functions.ab5_coefficients().
        Allows Numba support for the evaluation of the coupling matrix of periodic particle arrangements.
    Args:
        lmax1 (int):        l=1,...lmax: Original wave's SVWF multipole degree
        mmax1 (int):        m=-mmax,...,mmax: Original wave's SVWF multipole order
        lmax2 (int):        l=1,...lmax: Partial wave's SVWF multipole degree
        mmax2 (int):        m=-mmax,...,mmax: Partial wave's SVWF multipole order
    Returns: 
        Lookup table for a5 and b5 coefficients. 
        The respective coefficients can be called via [a=0 b=1, l1 - 1, l2 - 1, m1 + mmax2, m2 + mmax2, p].
    """

    a5b5mat = np.zeros([2, lmax1, lmax2, 2*mmax1+1, 2*mmax2+1, lmax1+lmax2+1], np.complex128)
    for l1 in range(1, lmax1 + 1):
        for l2 in range(1, lmax2 + 1):
            for m1 in range(-min(mmax1, l1), min(mmax1, l1) + 1):
                for m2 in range(-min(mmax2, l2), min(mmax2, l2) + 1):
                    M = m2 - m1
                    for L in range(max(abs(l1 - l2), abs(M)), l1 + l2 + 1): # if L < abs(M) then p=0   
                        a5b5mat[:, l1 - 1, l2 - 1, m1 + mmax1, m2 + mmax2, L] = flds.transformations.ab5_coefficients(l2, m2, l1, m1, L)
    return a5b5mat


@njit() 
def transformation_coefficients_vwf(tau, l, m, pol, kp=None, kz=None, pilm=None, taulm=None, dagger=False):
    r""" Numba compatible copy of smuthi.vector_wave_functions.transformation_coefficients_vwf()
    Transformation coefficients B to expand SVWF in PVWF and vice versa:

    .. math::
        B_{\tau l m,j}(x) = -\frac{1}{\mathrm{i}^{l+1}} \frac{1}{\sqrt{2l(l+1)}} (\mathrm{i} \delta_{j1} + \delta_{j2})
        (\delta_{\tau j} \tau_l^{|m|}(x) + (1-\delta_{\tau j} m \pi_l^{|m|}(x))

    For the definition of the :math:`\tau_l^m` and :math:`\pi_l^m` functions, see
    `A. Doicu, T. Wriedt, and Y. A. Eremin: "Light Scattering by Systems of Particles", Springer-Verlag, 2006
    <https://doi.org/10.1007/978-3-540-33697-6>`_

    Args:
        tau (int):          SVWF polarization, 0 for spherical TE, 1 for spherical TM
        l (int):            l=1,... SVWF multipole degree
        m (int):            m=-l,...,l SVWF multipole order
        pol (int):          PVWF polarization, 0 for TE, 1 for TM
        kp (numpy array):         PVWF in-plane wavenumbers
        kz (numpy array):         complex numpy-array: PVWF out-of-plane wavenumbers
        pilm (numpy array):       pilm[l, m, *ct.shape]
        taulm (numpy array):      tauilm[l, m, *ct.shape]
        dagger (bool):      switch on when expanding PVWF in SVWF and off when expanding SVWF in PVWF

    Returns:
        Transformation coefficient as array (size like kp).
    """
    # numba cannot handle these "is None" questions
    # if pilm is None: 
    #     k = np.sqrt(kp**2 + kz**2)
    #     ct = kz / k
    #     st = kp / k
    #     _, pilm, taulm = sf.legendre_normalized_numbed(ct, st, l)

    if tau == pol:
        sphfun = taulm[l, abs(m)]
    else:
        sphfun = m * pilm[l, abs(m)]

    if dagger:
        if pol == 0:
            prefac = -1 / (-1j) ** (l + 1) / np.sqrt(2 * l * (l + 1)) * (-1j)
        elif pol == 1:
            prefac = -1 / (-1j) ** (l + 1) / np.sqrt(2 * l * (l + 1)) * 1
        else:
            raise ValueError('pol must be 0 (TE) or 1 (TM)')
    else:
        if pol == 0:
            prefac = -1 / (1j) ** (l + 1) / np.sqrt(2 * l * (l + 1)) * (1j)
        elif pol ==1:
            prefac = -1 / (1j) ** (l + 1) / np.sqrt(2 * l * (l + 1)) * 1
        else:
            raise ValueError('pol must be 0 (TE) or 1 (TM)')

    B = prefac * sphfun
    return B


@njit()
def index_block(i, lmax_array, mmax_array):
    """ Numba compatible copy of smuthi.linearsystem.linear_system.System_matrix.index_block()
    Args:
        i (int):                        number of particle
        lmax_array (numpy.ndarray):     lmax of each particle 
        mmax_array (numpy.ndarray):     mmax of each particle
    Returns:
        Indices that correspond to the coefficients for that particle   
    """
    blocksizes = np.zeros((i + 1), np.int32)
    for idx in range(i + 1):
        blocksizes[idx] = flds.blocksize(lmax_array[idx], mmax_array[idx])
    return (np.sum(blocksizes[:i]), np.sum(blocksizes))


###############################################################################
#                  Numba compatible smuthi.layer functions                    #                               
###############################################################################
@njit()
def fresnel_r(pol, kz1, kz2, n1, n2):
    """ Numba compatible copy of smuthi.layers.fresnel_r()
        Fresnel reflection coefficient.
    Args:
        pol (int):              polarization (0=TE, 1=TM)
        kz1 (float or array):   incoming wave's z-wavenumber (k*cos(alpha1))
        kz2 (float or array):   transmitted wave's z-wavenumber (k*cos(alpha2))
        n1 (float or complex):  first medium's complex refractive index (n+ik)
        n2 (float or complex):  second medium's complex refractive index (n+ik)
    Returns:
        Complex Fresnel reflection coefficient (float or array)
    """
    if pol == 0:
        return (kz1 - kz2) / (kz1 + kz2)
    else:
        return (n2 ** 2 * kz1 - n1 ** 2 * kz2) / (n2 ** 2 * kz1 + n1 ** 2 * kz2)
    

@njit()    
def fresnel_t(pol, kz1, kz2, n1, n2):
    """ Numba compatible copy of smuthi.layers.fresnel_t()
        Fresnel transmission coefficient.
    Args:
        pol (int):              polarization (0=TE, 1=TM)
        kz1 (float or array):   incoming wave's z-wavenumber (k*cos(alpha1))
        kz2 (float or array):   transmitted wave's z-wavenumber (k*cos(alpha2))
        n1 (float or complex):  first medium's complex refractive index (n+ik)
        n2 (float or complex):  second medium's complex refractive index (n+ik)
    Returns:
        Complex Fresnel transmission coefficient (float or array)
    """
    if pol == 0:
        return 2 * kz1 / (kz1 + kz2)
    else:
        return 2 * n1 * n2 * kz1 / (n2 ** 2 * kz1 + n1 ** 2 * kz2)    


@njit()
def interface_transition_matrix(pol, kz1, kz2, n1, n2):
    """ Numba compatible copy of smuthi.layers.interface_transition_matrix()
        Interface transition matrix to be used in the Transfer matrix algorithm.
    Args:
        pol (int):              polarization (0=TE, 1=TM)
        kz1 (float or complex):   incoming wave's z-wavenumber (k*cos(alpha1))
        kz2 (float or complex):   transmitted wave's z-wavenumber (k*cos(alpha2))
        n1 (float or complex):  first medium's complex refractive index (n+ik)
        n2 (float or complex):  second medium's complex refractive index (n+ik)
    Returns:
        Interface transition matrix as 2x2 numpy array 
    """
    t = fresnel_t(pol, kz1, kz2, n1, n2)
    r = fresnel_r(pol, kz1, kz2, n1, n2)
    return 1 / t * np.array(((1, r), (r, 1)), np.complex128)


@njit()
def layer_propagation_matrix(kz, d):
    """ Numba compatible copy of smuthi.layers.layer_propagation_matrix()
        Layer propagation matrix to be used in the Transfer matrix algorithm.
    Args:
        kz (float or complex):  z-wavenumber (k*cos(alpha))
        d  (float):             thickness of layer
    Returns:
        Layer propagation matrix as 2x2 numpy array 
    """
    arr = np.array(((np.exp(-1j * kz * d), 0), (0, np.exp(1j * kz * d))), np.complex128)    
    # remove inf (and nan), inf appears for thick layers and evanescent waves
    if np.isinf(arr[0, 0]):
        arr[0, 0] = 1.79769313e+307
        arr[1, 1] = 1e-308       
    return arr
    

@njit()
def layersystem_transfer_matrix(pol, layer_d, layer_n, kpar, omega):
    """ Numba compatible copy of smuthi.layers.layersystem_transfer_matrix()
        Transfer matrix of a planarly layered medium.
    Args:
        pol (int):                  polarization(0=TE, 1=TM)
        layer_d (numpy.ndarray):    layer thicknesses (float64)
        layer_n (numpy.ndarray):    complex layer refractive indices
        kpar (complex):               in-plane wavenumber
        omega (float):              angular frequency in units of c=1: omega=2*pi/lambda
    Returns:
        Transfer matrix as 2x2 numpy array 
    """
    layer_kz = []
    for n in layer_n:
        kz = np.sqrt((omega * n) ** 2 - kpar ** 2 + 0j)
        if kz.imag < 0:
            kz = -kz
        layer_kz.append(kz)
    tmat = np.array(((1, 0), (0, 1)), np.complex128)
    for i in range(layer_d.shape[0] - 1):
        dmat = interface_transition_matrix(pol, layer_kz[i], layer_kz[i + 1], layer_n[i], layer_n[i + 1])
        pmat = layer_propagation_matrix(layer_kz[i], layer_d[i])
        tmat = np.dot(tmat, np.dot(pmat, dmat))
    return tmat


@njit()
def layersystem_scattering_matrix(pol, layer_d, layer_n, kpar, omega):
    """ Numba compatible copy of smuthi.layers.layersystem_scattering_matrix()
        Scattering matrix of a planarly layered medium.
    Args:
        pol (int):                  polarization(0=TE, 1=TM)
        layer_d (numpy.ndarray):    layer thicknesses (float64)
        layer_n (numpy.ndarray):    complex layer refractive indices
        kpar (complex):             in-plane wavenumber
        omega (float):              angular frequency in units of c=1: omega=2*pi/lambda
    Returns:
        Scattering matrix as 2x2 numpy array 
    """
    layer_kz = []
    for n in layer_n:
        kz = np.sqrt((omega * n) ** 2 - kpar ** 2 + 0j)
        if kz.imag < 0:
            kz = -kz
        layer_kz.append(kz)
    smat = np.array(((1, 0), (0, 1)), np.complex128)
    for i in range(layer_d.shape[0] - 1):
        dmat = interface_transition_matrix(pol, layer_kz[i], layer_kz[i + 1], layer_n[i], layer_n[i + 1])
        pmat = layer_propagation_matrix(layer_kz[i], layer_d[i])
        tmat = np.dot(pmat, dmat)
        s11 = smat[0, 0] / (tmat[0, 0] - smat[0, 1] * tmat[1, 0])
        s12 = (smat[0, 1] * tmat[1, 1] - tmat[0, 1]) / (tmat[0, 0] - smat[0, 1] * tmat[1, 0])
        s21 = smat[1, 1] * tmat[1, 0] * s11 + smat[1, 0]
        s22 = smat[1, 1] * tmat[1, 0] * s12 + smat[1, 1] * tmat[1, 1]
        smat = np.array(((s11, s12), (s21, s22)), np.complex128)
    return smat


@njit()
def layersystem_response_matrix(pol, layer_d, layer_n, kpar, omega, fromlayer, tolayer):
    """ Numba compatible copy of smuthi.layers.layersystem_response_matrix()
        Layer system response matrix of a planarly layered medium.
    Args:
        pol (int):                          polarization(0=TE, 1=TM)
        layer_d (numpy.ndarray):            layer thicknesses (float64)
        layer_n (numpy.ndarray):            complex layer refractive indices
        kpar (numpy.ndarray):               in-plane wavenumber
        omega (float):                      angular frequency in units of c=1: omega=2*pi/lambda
        fromlayer (int):                    number of layer where the excitation is located
        tolayer (int):                      number of layer where the response is evaluated
    Returns:
        Layer system response matrix as 2x2xN array with kpar is array of len = N.
    """
    
    # if kpar is an array of only one value, the function returns a (2x2x1) matrix /different to the original function
    # necessary since numba expects the function to always return arrays of the same number of dimensions
    
    shape = len(kpar)
    if not shape == 1: # is kpar an array? then use recursive call to fill an 2 x 2 x N ndarray
        result = np.zeros((2, 2, shape), np.complex128)
        for i in range(shape):
            kp = np.array([kpar[i]], np.complex128)
            result[:, :, i] = layersystem_response_matrix(pol, layer_d, layer_n, kp,
                                                          omega, fromlayer, tolayer).reshape(2, 2)
        return result

    layer_d_above = np.append(np.array([0], np.float64),  layer_d[fromlayer:])
    layer_n_above = np.append(np.array([layer_n[fromlayer]], np.complex128),  layer_n[fromlayer:])
    smat_above = layersystem_scattering_matrix(pol, layer_d_above, layer_n_above, kpar[0], omega)
    layer_d_below = np.append(layer_d[:fromlayer], np.array([0], np.float64))
    layer_n_below = layer_n[:fromlayer + 1]
    smat_below = layersystem_scattering_matrix(pol, layer_d_below, layer_n_below, kpar[0], omega)
    lmat = np.dot(np.linalg.inv(np.array(((1, -smat_below[0, 1]), (-smat_above[1, 0], 1)), np.complex128)),
                          np.array(((0, smat_below[0, 1]), (smat_above[1, 0], 0)), np.complex128))
    if tolayer > fromlayer:
        tmat_fromto = layersystem_transfer_matrix(pol, layer_d[fromlayer:tolayer + 1], layer_n[fromlayer:tolayer + 1],
                                                  kpar[0], omega)
        lmat = np.dot(np.linalg.inv(tmat_fromto), lmat + np.array(((1, 0), (0, 0)), np.complex128))
    elif tolayer < fromlayer:
        tmat_fromto = layersystem_transfer_matrix(pol, layer_d[tolayer:fromlayer + 1], layer_n[tolayer:fromlayer + 1],
                                                  kpar[0], omega)
        lmat = np.dot(tmat_fromto, lmat + np.array(((0, 0), (0, 1)), np.complex128))
    return lmat.reshape(2, 2, 1)
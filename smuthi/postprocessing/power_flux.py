"""Manage post processing steps to evaluate power flux"""

import smuthi.fields
import numpy as np

def power_flux_through_zplane(vacuum_wavelength, z, upgoing_pwe=None, downgoing_pwe=None):
    """Evaluate time averaged power flux though a plane of z=const.

    Args:
        vacuum_wavelength (float):  Vacuum wavelength in length units.
        z (float):                  plane height z
        upgoing_pwe (PlaneWaveExpansion):   of kind "upgoing"
        downgoing_pwe (PlaneWaveExpansion): of kind "downgoing"

    Returns:
        Time averaged energy flux.

    """
    
    if not upgoing_pwe and not downgoing_pwe:
        return 0     
    if upgoing_pwe:
        if not upgoing_pwe.valid(np.array([0]), np.array([0]), np.array([z])):
            raise ValueError('Upgoing PlaneWaveExpansion is not valid.')
    if downgoing_pwe:
        if not downgoing_pwe.valid(np.array([0]), np.array([0]), np.array([z])):
            raise ValueError('Downgoing PlaneWaveExpansion is not valid.')
    
    if upgoing_pwe and downgoing_pwe:
        if upgoing_pwe.kind == downgoing_pwe.kind:
            raise ValueError('Both PlaneWaveExpansions are of same kind.')
        if not np.isclose(upgoing_pwe.k, downgoing_pwe.k):
            raise ValueError('Wavenumbers are not compatible.')
        if not all(np.isclose(upgoing_pwe.k_parallel, downgoing_pwe.k_parallel)):
            raise ValueError('In-plane wavenumbers are not compatible.')
        if not all(np.isclose(upgoing_pwe.azimuthal_angles, downgoing_pwe.azimuthal_angles)):
            raise ValueError('Azimuthal angles are not compatible.')
            
    if upgoing_pwe:
        k = upgoing_pwe.k
        kpar = upgoing_pwe.k_parallel
        alpha = upgoing_pwe.azimuthal_angles
    else:
        k = downgoing_pwe.k
        kpar = downgoing_pwe.k_parallel
        alpha = downgoing_pwe.azimuthal_angles        
    omega = smuthi.fields.angular_frequency(vacuum_wavelength)
    kz = smuthi.fields.k_z(k_parallel=kpar, k=k)
           
    if z == float('-inf') or z == float('inf'): # no evanescent waves necessary
        if upgoing_pwe:
            kpar = upgoing_pwe.k_parallel[np.where(upgoing_pwe.k_parallel <= upgoing_pwe.k)[0]]
            gpl = upgoing_pwe.coefficients[:, np.where(upgoing_pwe.k_parallel <= upgoing_pwe.k)[0], :]
        if downgoing_pwe:
            kpar = downgoing_pwe.k_parallel[np.where(downgoing_pwe.k_parallel <= downgoing_pwe.k)[0]]
            gmn = downgoing_pwe.coefficients[:, np.where(downgoing_pwe.k_parallel <= downgoing_pwe.k)[0], :]
        kz = smuthi.fields.k_z(k_parallel=kpar, k=k)
    else: # translate to z
        if upgoing_pwe:
            ejkzdz = np.exp(1j * upgoing_pwe.k_z_grid() * (z - upgoing_pwe.reference_point[2]))
            gpl = upgoing_pwe.coefficients * ejkzdz
        if downgoing_pwe:
            ejkzdz = np.exp(1j * downgoing_pwe.k_z_grid() * (z - downgoing_pwe.reference_point[2]))
            gmn = downgoing_pwe.coefficients * ejkzdz
    
    if any(kpar.imag):
        raise ValueError('complex in-plane wavenumbers are not allowed')
    
    if not upgoing_pwe:
        gpl = np.zeros(gmn.shape, dtype=complex)
    if not downgoing_pwe:
        gmn = np.zeros(gpl.shape, dtype=complex)         
    
    prefac = 2 * np.pi ** 2 / omega            
    gpl2_mn_gmn2 = np.abs(gpl) ** 2 - np.abs(gmn) ** 2
    im_gpl_t_gmnconj = (gpl * np.conjugate(gmn)).imag
    integrand = kz.real[None, : , None] * gpl2_mn_gmn2 - 2 * kz.imag[None, :, None] * im_gpl_t_gmnconj           
    return prefac * np.sum(np.trapz(np.trapz(integrand, alpha, axis=2) * kpar, kpar, axis=1))


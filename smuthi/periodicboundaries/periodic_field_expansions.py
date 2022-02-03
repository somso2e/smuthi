"""This module contains post processing functions to evaluate the scattered
electromagnetic field of periodic particle arrangements."""

import smuthi.fields as flds
import smuthi.fields.expansions as fldex
from smuthi.fields.transformations import transformation_coefficients_vwf
from smuthi.postprocessing.far_field import FarField
import smuthi.periodicboundaries.periodic_particle_coupling as ppc
import numpy as np


def periodic_swe_to_pwe_conversion(initial_field, layer_system, particle, a1, a2):
    """ Converts the scattered field's SphericalWaveExpansion of the arbitrary central particle within
        a periodic extend of particles into a PlaneWaveExpansion object and translates it to
        the scattering layer's anchor point.
    Args:
        initial_field (smuthi.initial_field.PlaneWave):     initial plane wave object
        layer_system (smuthi.layer.LayerSystem):            stratified medium
        particle (smuthi.particles.Particle):               particle (periodic arrangement) emitting scattered field
        a1 (numpy.ndarray):                                 lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):                                 lattice vector 2 in Carthesian coordinates 
    Returns:
        pwe_up (smuthi.fields.expansions.PlaneWaveExpansion):       upwards propagating plane wave expansion of the
                                                                    scattered field of a periodic extend of a single particle 
        pwe_down (smuthi.fields.expansions.PlaneWaveExpansion):     downwards propagating plane wave expansion of the
                                                                    scattered field of a periodic extend of a single particle 
    """
    
    # prepare kparallel and alpha discretization of the PWE
    i_swe = layer_system.layer_number(particle.position[2]) 
    reference_point = [0, 0, layer_system.reference_z(i_swe)]
    if initial_field.polar_angle < np.pi:
        pfe = initial_field.piecewise_field_expansion(layer_system).expansion_list[2 * i_swe]
    else:
        pfe = initial_field.piecewise_field_expansion(layer_system).expansion_list[2 * i_swe + 1]
    k = pfe.k
    k0t = np.array([pfe.k_parallel[0] * np.cos(pfe.azimuthal_angles),
                    pfe.k_parallel[0] * np.sin(pfe.azimuthal_angles)])
    
    A = np.linalg.norm(np.cross(a1, a2))    
    k_parallel, azimuthal_angles, kpar_alpha_comb, gmax = ppc.discrete_scattering_angles(k, k0t, pfe.azimuthal_angles[0], a1, a2)
    
    # prepare pwe objects    
    lower_z_up = particle.position[2]
    upper_z_up = layer_system.upper_zlimit(i_swe)      
    pwe_up = fldex.PlaneWaveExpansion(k=k, k_parallel=k_parallel, azimuthal_angles=azimuthal_angles, kind='upgoing',
                                      reference_point=reference_point, lower_z=lower_z_up, upper_z=upper_z_up)
    lower_z_down = layer_system.lower_zlimit(i_swe)
    upper_z_down = particle.position[2]
    pwe_down = fldex.PlaneWaveExpansion(k=k, k_parallel=k_parallel, azimuthal_angles=azimuthal_angles, kind='downgoing',
                                        reference_point=reference_point, lower_z=lower_z_down, upper_z=upper_z_down)
    
    # prepare kpar and alpha lookups
    alpha_vec = np.zeros((2 * gmax + 1) ** 2)
    kpar_vec = np.zeros((2 * gmax + 1) ** 2)
    k_idx = np.zeros((2 * gmax + 1) ** 2, dtype=int)
    a_idx = np.zeros((2 * gmax + 1) ** 2, dtype=int)
    s = 0
    for idx in range(len(k_parallel)):
        alpha_vec[s:s+len(kpar_alpha_comb[idx])] = azimuthal_angles[kpar_alpha_comb[idx]]
        kpar_vec[s:s+len(kpar_alpha_comb[idx])] = k_parallel[idx]
        k_idx[s:s+len(kpar_alpha_comb[idx])] = idx
        a_idx[s:s+len(kpar_alpha_comb[idx])] = kpar_alpha_comb[idx]
        s += len(kpar_alpha_comb[idx])
        
    kz_up = flds.k_z(k_parallel=kpar_vec, k=pwe_up.k)
    
    kzk = kz_up * k 
    for l in range(1, particle.l_max + 1):
        for m in range(-min(particle.m_max, l), min(particle.m_max, l) + 1):
            for tau in range(2):
                n = flds.multi_to_single_index(tau, l, m, particle.l_max, particle.m_max)
                eima = np.exp(1j * m * alpha_vec)
                for pol in range(2):
                    B_up = transformation_coefficients_vwf(tau, l, m, pol=pol, kp=kpar_vec, kz=kz_up,
                                                           pilm_list=None, taulm_list=None, dagger=False) 
                    B_down = transformation_coefficients_vwf(tau, l, m, pol=pol, kp=kpar_vec, kz=-kz_up,
                                                             pilm_list=None, taulm_list=None, dagger=False) 
                    pwe_up.coefficients[pol, k_idx, a_idx] += particle.scattered_field.coefficients[n] * eima * B_up
                    pwe_down.coefficients[pol, k_idx, a_idx] += particle.scattered_field.coefficients[n] * eima * B_down
    pwe_up.coefficients[:, k_idx, a_idx] = 2 * np.pi / (A * kzk) * pwe_up.coefficients[:, k_idx, a_idx]
    pwe_down.coefficients[:, k_idx, a_idx] = 2 * np.pi / (A * kzk) * pwe_down.coefficients[:, k_idx, a_idx]
    
    # translate the pwe objects to the reference point
    rpwe_mn_rswe = np.array(reference_point) - np.array(particle.scattered_field.reference_point)
    
    agrid = pwe_up.azimuthal_angle_grid()
    kpgrid = pwe_up.k_parallel_grid()
    kx = kpgrid * np.cos(agrid)
    ky = kpgrid * np.sin(agrid)
    kz_up = pwe_up.k_z_grid()
    kz_down = pwe_down.k_z_grid()

    kvec_up = np.array([kx, ky, kz_up])
    kvec_down = np.array([kx, ky, kz_down])
    ejkrSiS_up = np.nan_to_num(np.exp(1j * np.tensordot(kvec_up, rpwe_mn_rswe, axes=([0], [0]))))
    ejkrSiS_down = np.exp(1j * np.tensordot(kvec_down, rpwe_mn_rswe, axes=([0], [0])))
    
    pwe_up.coefficients = pwe_up.coefficients * ejkrSiS_up[None, :, :]
    pwe_down.coefficients = pwe_down.coefficients * ejkrSiS_down[None, :, :]
        
    return pwe_up, pwe_down    

# Note: The FarField.signal corresponds to the plane wave's power per interface area (Theobald 2021 dissertation, eq.(6.29))
# It differs from the smuthi.postprocessing.far_field.FarField object that is constructed from a PlaneWaveExpansion
# with a continuous angular distribution.
# Not all methods of FarField objects are applicable! Might be a good idea to open up a new class.
def periodic_pwe_to_ff_conversion(initial_field, layer_system, plane_wave_expansion):
    """Compute the far field of a plane wave expansion object of discrete plane waves.

    Args:
        initial_field (smuthi.initial_field.PlaneWave):     initial plane wave object
        layer_system (smuthi.layer.LayerSystem):            stratified medium
        plane_wave_expansion (smuthi.fields.expansions.PlaneWaveExpansion):     
                                                            plane wave expansion of the scattered field of a
                                                            periodic extend of a single particle that is to be
                                                            convert into a periodic far field object.
    Returns:
        A smuthi.postprocessing.far_field.FarField object.
    """
    omega = flds.angular_frequency(initial_field.vacuum_wavelength)
    if 0 <= initial_field.polar_angle <= np.pi/2:
        n0 = layer_system.refractive_indices[0]
    elif np.pi / 2 <= initial_field.polar_angle < np.pi:
        n0 = layer_system.refractive_indices[-1]
    k0 = 2 * np.pi * n0 / initial_field.vacuum_wavelength
    k = plane_wave_expansion.k
    kp = plane_wave_expansion.k_parallel
    if plane_wave_expansion.kind == 'upgoing':
        polar_angles = np.arcsin(kp / k)
    elif plane_wave_expansion.kind == 'downgoing':
        polar_angles = np.pi - np.arcsin(kp / k)
    else:
        raise ValueError('PWE type not specified')
    if any(polar_angles.imag):
        raise ValueError('complex angles are not allowed')
    polar_angles = polar_angles[np.logical_not(np.isnan(polar_angles))]
    azimuthal_angles = plane_wave_expansion.azimuthal_angles
    
    intens = (k0 / (2 * omega) * np.cos(polar_angles)[None, :, None] \
              * abs(plane_wave_expansion.coefficients[:, :len(polar_angles), :]) ** 2)
        
    srt_idcs = np.argsort(polar_angles)  # reversing order in case of downgoing
    ff = FarField(polar_angles=polar_angles[srt_idcs], azimuthal_angles=azimuthal_angles)
    ff.signal = intens[:, srt_idcs, :]
    return ff

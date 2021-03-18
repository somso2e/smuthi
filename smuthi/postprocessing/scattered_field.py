"""Manage post processing steps to evaluate the scattered electric field"""

import sys
from tqdm import tqdm
import smuthi.fields as flds
import smuthi.fields.expansions as fldex
import smuthi.fields.transformations as trf
import numpy as np
from smuthi.layers import layersystem_response_matrix
from smuthi.fields.vector_wave_functions import plane_vector_wave_function
from smuthi.fields.transformations import transformation_coefficients_vwf

def scattered_field_piecewise_expansion(vacuum_wavelength, particle_list, layer_system, k_parallel='default',
                                        azimuthal_angles='default', angular_resolution=None, layer_numbers=None):
    """Compute a piecewise field expansion of the scattered field.

    Args:
        vacuum_wavelength (float):                  vacuum wavelength
        particle_list (list):                       list of smuthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem):   stratified medium
        k_parallel (numpy.ndarray or str):          in-plane wavenumbers array.
                                                    if 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
        azimuthal_angles (numpy.ndarray or str):    azimuthal angles array
                                                    if 'default', use smuthi.fields.default_azimuthal_angles
        angular_resolution (float):                 If provided, angular arrays are generated with this angular
                                                    resolution over the default angular range
        layer_numbers (list):                       if specified, append only plane wave expansions for these layers


    Returns:
        scattered field as smuthi.field_expansion.PiecewiseFieldExpansion object

    """
    if layer_numbers is None:
        layer_numbers = range(layer_system.number_of_layers())

    if type(k_parallel) == str and k_parallel == 'default':
        k_parallel = flds.default_Sommerfeld_k_parallel_array

    if angular_resolution is not None:
        azimuthal_angles, _ = flds.angular_arrays(angular_resolution)
    if type(azimuthal_angles) == str and azimuthal_angles == 'default':
        azimuthal_angles = flds.default_azimuthal_angles

    sfld = fldex.PiecewiseFieldExpansion()
    for i in layer_numbers:
        # layer mediated scattered field ---------------------------------------------------------------------------
        k = flds.angular_frequency(vacuum_wavelength) * layer_system.refractive_indices[i]
        ref = [0, 0, layer_system.reference_z(i)]
        vb = (layer_system.lower_zlimit(i), layer_system.upper_zlimit(i))
        pwe_up = fldex.PlaneWaveExpansion(k=k, k_parallel=k_parallel, azimuthal_angles=azimuthal_angles, kind='upgoing',
                                          reference_point=ref, lower_z=vb[0], upper_z=vb[1])
        pwe_down = fldex.PlaneWaveExpansion(k=k, k_parallel=k_parallel, azimuthal_angles=azimuthal_angles,
                                            kind='downgoing', reference_point=ref, lower_z=vb[0], upper_z=vb[1])
        
        for particle in tqdm(particle_list, desc='Scatt. field expansion (%i)'%i, file=sys.stdout,
                                        bar_format='{l_bar}{bar}| elapsed: {elapsed} ' 'remaining: {remaining}'):
            add_up, add_down = trf.swe_to_pwe_conversion(particle.scattered_field, k_parallel, azimuthal_angles,
                                                         layer_system, i, True)
            pwe_up = pwe_up + add_up
            pwe_down = pwe_down + add_down

            pwe_up.validity_conditions.append(particle.is_outside)
            pwe_down.validity_conditions.append(particle.is_outside)

        # in bottom_layer, suppress upgoing waves, and in top layer, suppress downgoing waves
        if i > 0:
            sfld.expansion_list.append(pwe_up)
        if i < layer_system.number_of_layers()-1:
            sfld.expansion_list.append(pwe_down)

    # direct field ---------------------------------------------------------------------------------------------
    for particle in particle_list:
        sfld.expansion_list.append(particle.scattered_field)
        sfld.validity_conditions.append(particle.is_outside)

    return sfld


def scattered_field_pwe(vacuum_wavelength, particle_list, layer_system, layer_number, k_parallel='default',
                        azimuthal_angles='default', angular_resolution=None, include_direct=True, include_layer_response=True,
                        only_l=None, only_m=None, only_pol=None, only_tau=None):
    """Calculate the plane wave expansion of the scattered field of a set of particles.

    Args:
        vacuum_wavelength (float):          Vacuum wavelength (length unit)
        particle_list (list):               List of Particle objects
        layer_system (smuthi.layers.LayerSystem):  Stratified medium
        layer_number (int):                 Layer number in which the plane wave expansion should be valid
        k_parallel (numpy.ndarray or str):      in-plane wavenumbers array.
                                                if 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
        azimuthal_angles (numpy.ndarray or str):azimuthal angles array
                                                if 'default', use smuthi.fields.default_azimuthal_angles
        angular_resolution (float):             If provided, angular arrays are generated with this angular resolution
                                                over the default angular range
        include_direct (bool):                  If True, include the direct scattered field
        include_layer_response (bool):          If True, include the layer system response
        only_pol (int):  if set to 0 or 1, only this plane wave polarization (0=TE, 1=TM) is considered
        only_tau (int):  if set to 0 or 1, only this spherical vector wave polarization (0 — magnetic, 1 — electric) is
                         considered
        only_l (int):    if set to positive number, only this multipole degree is considered
        only_m (int):    if set to non-negative number, only this multipole order is considered

    Returns:
        A tuple of PlaneWaveExpansion objects for upgoing and downgoing waves.
    """

    sys.stdout.write('Evaluating scattered field plane wave expansion in layer number %i ...\n'%layer_number)
    sys.stdout.flush()

    if type(k_parallel) == str and k_parallel == 'default':
        k_parallel = flds.default_Sommerfeld_k_parallel_array

    if angular_resolution is not None:
        azimuthal_angles, _ = flds.angular_arrays(angular_resolution)
    if type(azimuthal_angles) == str and azimuthal_angles == 'default':
        azimuthal_angles = flds.default_azimuthal_angles

    omega = flds.angular_frequency(vacuum_wavelength)
    k = omega * layer_system.refractive_indices[layer_number]
    z = layer_system.reference_z(layer_number)
    vb = (layer_system.lower_zlimit(layer_number), layer_system.upper_zlimit(layer_number))
    pwe_up = fldex.PlaneWaveExpansion(k=k, k_parallel=k_parallel, azimuthal_angles=azimuthal_angles, kind='upgoing',
                                      reference_point=[0, 0, z], lower_z=vb[0], upper_z=vb[1])
    pwe_down = fldex.PlaneWaveExpansion(k=k, k_parallel=k_parallel, azimuthal_angles=azimuthal_angles, kind='downgoing',
                                        reference_point=[0, 0, z], lower_z=vb[0], upper_z=vb[1])

    for iS, particle in enumerate(tqdm(particle_list, desc='Scatt. field pwe          ', file=sys.stdout,
                                        bar_format='{l_bar}{bar}| elapsed: {elapsed} ' 'remaining: {remaining}')):

        i_iS = layer_system.layer_number(particle.position[2])

        # direct contribution
        if i_iS == layer_number and include_direct:
            pu, pd = trf.swe_to_pwe_conversion(swe=particle.scattered_field, k_parallel=k_parallel,
                                               azimuthal_angles=azimuthal_angles, layer_system=layer_system,
                                               only_l=only_l, only_m=only_m, only_pol=only_pol, only_tau=only_tau)
            pwe_up = pwe_up + pu
            pwe_down = pwe_down + pd

        # layer mediated contribution
        if include_layer_response:
            pu, pd = trf.swe_to_pwe_conversion(swe=particle.scattered_field, k_parallel=k_parallel,
                                               azimuthal_angles=azimuthal_angles, layer_system=layer_system,
                                               layer_number=layer_number, layer_system_mediated=True,
                                               only_l=only_l, only_m=only_m, only_pol=only_pol, only_tau=only_tau)
            pwe_up = pwe_up + pu
            pwe_down = pwe_down + pd

    return pwe_up, pwe_down


def evaluate_scattered_field_stat_phase_approx(x, y, z, vacuum_wavelength, particle_list, layer_system):
    """Evaluate the scattered electric field for N particles on a substrate. The substrate reflection is evaluated
    by means of the stationary phase approximation, as presented in
    "A quick way to approximate a Sommerfeld-Weyl_type Sommerfeld integral" by W.C. Chew (1988).

    See also the technical note "Usage of the stationary phase approximation in SMUTHI" by A. Egel (2020)

    The stationary phase approximation is expected to yield good results for field points far away from the particles.

    ********************************************************************************************************************
    Note: This function assumes that the particles are located in the upper layer of a two-layer system (particles on
    substrate). For other cases, this function does not apply.
    ********************************************************************************************************************

    Args:
        x (float or numpy.ndarray):                                  x-coordinates of query points
        y (float or numpy.ndarray):                                  y-coordinates of query points
        z (float or numpy.ndarray):                                  z-coordinates of query points
        vacuum_wavelength (float):                          Vacuum wavelength :math:`\lambda` (length unit)
        particle_list (list):                               List of Particle objects
        layer_system (smuthi.layers.LayerSystem):           Stratified medium

    Returns:
        Tuple of (E_x, E_y, E_z) numpy.ndarray objects with the Cartesian coordinates of complex electric field.
    """
    x = np.array(x, ndmin=1)
    y = np.array(y, ndmin=1)
    z = np.array(z, ndmin=1)

    # assert a two layer system
    assert (len(layer_system.refractive_indices) == 2)

    # assert that all particles are in top layer (this function applies only to this special case)
    for particle in particle_list:
        assert (layer_system.layer_number(particle.position[2]) == 1)

    # assert that all query points are in top layer (this function applies only to this special case)
    for zi in z:
        assert (layer_system.layer_number(zi) == 1)

    ex = np.zeros(x.shape, dtype=complex)
    ey = np.zeros(x.shape, dtype=complex)
    ez = np.zeros(x.shape, dtype=complex)

    k = layer_system.wavenumber(1, vacuum_wavelength)

    for particle in tqdm(particle_list, desc='Scatt. field st. ph. appr.', file=sys.stdout,
                                        bar_format='{l_bar}{bar}| elapsed: {elapsed} ' 'remaining: {remaining}'):

        # direct scattered field
        exdir, eydir, ezdir = particle.scattered_field.electric_field(x, y, z)
        ex += exdir
        ey += eydir
        ez += ezdir

        # layer-reflected scattered field
        xtilde = x - particle.position[0]
        ytilde = y - particle.position[1]
        ztilde = z + particle.position[2]
        rtilde = np.sqrt(xtilde**2 + ytilde**2 + ztilde**2)

        # stationary phase wave-vector
        k0x = k * xtilde / rtilde
        k0y = k * ytilde / rtilde
        k0z = k * ztilde / rtilde
        k0par = np.hypot(k0x, k0y)
        alpha0 = np.arctan2(k0y, k0x)

        prefac = -1j * np.exp(1j * k * rtilde) / (k * rtilde)

        for pol in range(2):
            L = layersystem_response_matrix(pol=pol,
                                            layer_d=layer_system.thicknesses,
                                            layer_n=layer_system.refractive_indices,
                                            kpar=k0par,
                                            omega=flds.angular_frequency(vacuum_wavelength),
                                            fromlayer=1,
                                            tolayer=1)
            fresnel_coeff = L[0, 1]
            pvwf_direction = plane_vector_wave_function(0, 0, 0, k0par, alpha0, k0z, pol)
            print(pvwf_direction[2][-1])

            for tau in range(2):
                for m in range(-particle.scattered_field.m_max, particle.scattered_field.m_max + 1):
                    eima = np.exp(1j * m * alpha0)
                    for l in range(max(1, abs(m)), particle.scattered_field.l_max + 1):
                        b = particle.scattered_field.coefficients_tlm(tau, l, m)
                        Bnj = transformation_coefficients_vwf(tau=tau, l=l, m=m, pol=pol, kp=k0par, kz=-k0z)
                        scalar_fac = prefac * fresnel_coeff * b * eima * Bnj
                        ex += scalar_fac * pvwf_direction[0]
                        ey += scalar_fac * pvwf_direction[1]
                        ez += scalar_fac * pvwf_direction[2]

    return ex, ey, ez

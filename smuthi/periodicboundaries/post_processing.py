"""This module contains post processing functions to evaluate the scattered
electromagnetic field of periodic particle arrangements."""

import smuthi.fields as flds
import smuthi.fields.expansions as fldex
from smuthi.fields.transformations import transformation_coefficients_vwf
from smuthi.postprocessing.far_field import FarField
import smuthi.utility.cuda as cu
try:
    import pycuda.autoinit
    import pycuda.driver as drv
    from pycuda import gpuarray
    from pycuda.compiler import SourceModule
    import pycuda.cumath
except:
    pass
import smuthi.periodicboundaries.particle_coupling as pbcoup
import smuthi.periodicboundaries.expansions_cuda as cu_src
import numpy as np
import copy


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
    k_parallel, azimuthal_angles, kpar_alpha_comb, gmax = pbcoup.discrete_scattering_angles(k, k0t, pfe.azimuthal_angles[0], a1, a2)
    
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


def total_field_periodic_plane_wave_expansion(layer_system, initial_field, particle_list, a1, a2):
    """ Compute plane wave expansion of the total transmitted and the total reflected field.
    Args:
        layer_system (smuthi.layer.LayerSystem):            stratified medium
        initial_field (smuthi.initial_field.PlaneWave):     initial plane wave object
        particle_list (list):                               list of particles within one unit cell
        a1 (numpy.ndarray):                                 lattice vector 1 in carthesian coordinates
        a2 (numpy.ndarray):                                 lattice vector 2 in carthesian coordinates 
    Returns:
        Tuple of (transmitted, reflected) plane wave propagating through a planarly layered medium containing 
        periodic particle arrangements
    """

    # evaluate each periodic particle's scattered plane wave expansions   
    pfe_up, pfe_down = [], []
    for particle in particle_list:
        pwe_up, pwe_down = periodic_swe_to_pwe_conversion(initial_field=initial_field,
                                                              layer_system=layer_system,
                                                              particle=particle,
                                                              a1=a1, a2=a2)
        pfe_up.append(pwe_up)
        pfe_down.append(pwe_down)

    # prepare piecewise field expansion of the initial field
    pfe_init = initial_field.piecewise_field_expansion(layer_system)
    pfe_init_top = pfe_init.expansion_list[-2:]
    pfe_init_bottom = pfe_init.expansion_list[:2]
    
    # prepare piecewise field expansion of the scattered field 
    for idx in range(len(pfe_up)):
        pwe_up, pwe_down = pfe_up[idx], pfe_down[idx]
        i_sca = layer_system.layer_number(pwe_up.reference_point[2])
        i_max = layer_system.number_of_layers() - 1
        if i_sca == 0:
            pwe_up_top, __ = layer_system.response((pwe_up, pwe_down), i_sca, i_max)
            pwe_down_bot = pwe_down + layer_system.response((pwe_up, pwe_down), i_sca, 0)[1]
        elif i_sca == i_max:
            pwe_up_top = pwe_up + layer_system.response((pwe_up, pwe_down), i_sca, i_max)[0]
            __, pwe_down_bot = layer_system.response((pwe_up, pwe_down), i_sca, 0)
        else:
            pwe_up_top, __ = layer_system.response((pwe_up, pwe_down), i_sca, i_max)
            __, pwe_down_bot = layer_system.response((pwe_up, pwe_down), i_sca, 0)
            
        if idx == 0:
            pfe_up_top = pwe_up_top
            pfe_down_bottom = pwe_down_bot
        else:
            pfe_up_top += pwe_up_top
            pfe_down_bottom += pwe_down_bot
            
    # assign pfe to Transmittance and Reflectance ############################     
    if 0 <= initial_field.polar_angle < np.pi / 2:
        pfe_sca_T, pfe_init_T = pfe_up_top, pfe_init_top[0]
        pfe_sca_R, pfe_init_R = pfe_down_bottom, pfe_init_bottom[1]
    else:
        pfe_sca_T, pfe_init_T = pfe_down_bottom, pfe_init_bottom[1]
        pfe_sca_R, pfe_init_R = pfe_up_top, pfe_init_top[0]
    
    # total fields pfe #######################################################
    pfe_in_list = [pfe_init_T, pfe_init_R]
    pfe_sca_list = [pfe_sca_T, pfe_sca_R]
    pfe_total_list = [copy.deepcopy(pfe_sca_T), copy.deepcopy(pfe_sca_R)]
    for idx in range(2):
        if pfe_in_list[idx].reference_point == pfe_sca_list[idx].reference_point:
            idx_a = np.argwhere(pfe_sca_list[idx].azimuthal_angles == pfe_in_list[idx].azimuthal_angles[0])[0][0]
            idx_kpar = np.argwhere(pfe_sca_list[idx].k_parallel == pfe_in_list[idx].k_parallel[0])[0][0]
        else:
            raise ValueError('Plane wave expansions are not compatible!')
    
        pfe_total_list[idx].coefficients[:, idx_kpar, idx_a] += pfe_in_list[idx].coefficients.reshape(2)
    
        # prevent rounding error @ kpar == k #################################
        pfe_total_list[idx].k_parallel[np.isclose(pfe_total_list[idx].k_parallel, pfe_total_list[idx].k)] = pfe_total_list[idx].k
    
    pfe_total_T = pfe_total_list[0]
    pfe_total_R = pfe_total_list[1]
    return pfe_total_T, pfe_total_R


# Field evaluation follows the procedure of Smuthi's "old" electric field
# and magnetic field routines 
def electromagnetic_nearfield(pwe, x, y, z, chunksize=50, field_type='electric', vacuum_wavelength=None):
    """Evaluate electric or magnetic field.      
    Args:
        pwe (smuthi.field_expansion.PlaneWaveExpansion):    plane wave expansion emitted by a periodic particle
                                                            arrangement 
        x (numpy.ndarray):    x-coordinates of query points
        y (numpy.ndarray):    y-coordinates of query points
        z (numpy.ndarray):    z-coordinates of query points
        chunksize (int):      number of field points that are simultaneously evaluated when running in CPU mode
        field_type (str):     field type to be evaluated; either 'electric' or 'magnetic'
        vacuum_wavelength (float):  vacuum wavelength in length units (only necessary for magnetic fields).    
    Returns:
        Tuple of (E_x, E_y, E_z) or (M_x, M_y, M_z) numpy.ndarray objects with the Cartesian coordinates of
        complex electric field.
    """    
    emx = np.zeros(x.shape, dtype=complex)
    emy = np.zeros(x.shape, dtype=complex)
    emz = np.zeros(x.shape, dtype=complex)

    xr = x[pwe.valid(x, y, z)] - pwe.reference_point[0]
    yr = y[pwe.valid(x, y, z)] - pwe.reference_point[1]
    zr = z[pwe.valid(x, y, z)] - pwe.reference_point[2]
    
    if vacuum_wavelength:
        omega = flds.angular_frequency(vacuum_wavelength)
    
    if cu.use_gpu and xr.size and len(pwe.k_parallel) > 1:  # run calculations on gpu
            
        re_kp_d = gpuarray.to_gpu(pwe.k_parallel.real.astype(np.float32))
        im_kp_d = gpuarray.to_gpu(pwe.k_parallel.imag.astype(np.float32))
        
        re_kz_d = gpuarray.to_gpu(pwe.k_z().real.astype(np.float32))
        im_kz_d = gpuarray.to_gpu(pwe.k_z().imag.astype(np.float32))
        
        alpha_d = gpuarray.to_gpu(pwe.azimuthal_angles.astype(np.float32))
        
        xr_d = gpuarray.to_gpu(xr.astype(np.float32))
        yr_d = gpuarray.to_gpu(yr.astype(np.float32))
        zr_d = gpuarray.to_gpu(zr.astype(np.float32))
        
        re_g_te_d = gpuarray.to_gpu(pwe.coefficients[0, :, :].real.astype(np.float32))
        im_g_te_d = gpuarray.to_gpu(pwe.coefficients[0, :, :].imag.astype(np.float32))
        re_g_tm_d = gpuarray.to_gpu(pwe.coefficients[1, :, :].real.astype(np.float32))
        im_g_tm_d = gpuarray.to_gpu(pwe.coefficients[1, :, :].imag.astype(np.float32))
        
        re_em_x_d = gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
        im_em_x_d = gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
        re_em_y_d = gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
        im_em_y_d = gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
        re_em_z_d = gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
        im_em_z_d = gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
        
        
        if field_type == 'electric':
            kernel_source = cu_src.pwe_periodic_electric_field_evaluation_code%(xr.size, len(pwe.k_parallel), 
                                                                   len(pwe.azimuthal_angles), (1/pwe.k).real, 
                                                                   (1/pwe.k).imag)
            kernel_function = SourceModule(kernel_source).get_function("electric_field") 
            
            cuda_blocksize = 128
            cuda_gridsize = (xr.size + cuda_blocksize - 1) // cuda_blocksize
            
            kernel_function(re_kp_d, im_kp_d, re_kz_d, im_kz_d, alpha_d, xr_d, yr_d, zr_d, re_g_te_d, im_g_te_d,
                            re_g_tm_d, im_g_tm_d, re_em_x_d, im_em_x_d, re_em_y_d, im_em_y_d, re_em_z_d, im_em_z_d,
                            block=(cuda_blocksize,1,1), grid=(cuda_gridsize,1))
             
            emx[pwe.valid(x, y, z)] = re_em_x_d.get() + 1j * im_em_x_d.get()
            emy[pwe.valid(x, y, z)] = re_em_y_d.get() + 1j * im_em_y_d.get()
            emz[pwe.valid(x, y, z)] = re_em_z_d.get() + 1j * im_em_z_d.get()
            
        elif field_type == 'magnetic':
            kernel_source = cu_src.pwe_periodic_magnetic_field_evaluation_code%(xr.size, len(pwe.k_parallel), 
                                                                   len(pwe.azimuthal_angles), (pwe.k).real, 
                                                                   (pwe.k).imag)
            kernel_function = SourceModule(kernel_source).get_function("magnetic_field")
            
            cuda_blocksize = 128
            cuda_gridsize = (xr.size + cuda_blocksize - 1) // cuda_blocksize
            
            kernel_function(re_kp_d, im_kp_d, re_kz_d, im_kz_d, alpha_d, xr_d, yr_d, zr_d, re_g_te_d, im_g_te_d,
                            re_g_tm_d, im_g_tm_d, re_em_x_d, im_em_x_d, re_em_y_d, im_em_y_d, re_em_z_d, im_em_z_d,
                            block=(cuda_blocksize,1,1), grid=(cuda_gridsize,1))
             
            emx[pwe.valid(x, y, z)] = 1 / omega * (re_em_x_d.get() + 1j * im_em_x_d.get())
            emy[pwe.valid(x, y, z)] = 1 / omega * (re_em_y_d.get() + 1j * im_em_y_d.get())
            emz[pwe.valid(x, y, z)] = 1 / omega * (re_em_z_d.get() + 1j * im_em_z_d.get())
        else:
            assert False, 'Demanded field type must be electric or magnetic!'
            
    else: # run calculation of CPU
        kpgrid = pwe.k_parallel_grid()
        agrid = pwe.azimuthal_angle_grid()
        kx = kpgrid * np.cos(agrid)
        ky = kpgrid * np.sin(agrid)
        kz = pwe.k_z_grid()
    
        em_x_flat = np.zeros(xr.size, dtype=np.complex64)
        em_y_flat = np.zeros(xr.size, dtype=np.complex64)
        em_z_flat = np.zeros(xr.size, dtype=np.complex64)
        
        for i_chunk in range(np.ceil(xr.size / chunksize)):
            chunk_idcs = range(i_chunk * chunksize, min((i_chunk + 1) * chunksize, xr.size))
            xr_chunk = xr.flatten()[chunk_idcs]
            yr_chunk = yr.flatten()[chunk_idcs]
            zr_chunk = zr.flatten()[chunk_idcs]
    
            kr = np.zeros((len(xr_chunk), len(pwe.k_parallel), len(pwe.azimuthal_angles)), dtype=np.complex64)
            kr += np.tensordot(xr_chunk, kx, axes=0)
            kr += np.tensordot(yr_chunk, ky, axes=0)
            kr += np.tensordot(zr_chunk, kz, axes=0)
    
            eikr = np.exp(1j * kr)
            
            summand_x = np.zeros((len(xr_chunk), len(pwe.k_parallel), len(pwe.azimuthal_angles)),
                                   dtype=np.complex64)
            summand_y = np.zeros((len(yr_chunk), len(pwe.k_parallel), len(pwe.azimuthal_angles)),
                                   dtype=np.complex64)
            summand_z = np.zeros((len(zr_chunk), len(pwe.k_parallel), len(pwe.azimuthal_angles)),
                                   dtype=np.complex64)
            
            if field_type == 'electric':
                # pol=0
                summand_x += (-np.sin(agrid) * pwe.coefficients[0, :, :])[None, :, :] * eikr
                summand_y += (np.cos(agrid) * pwe.coefficients[0, :, :])[None, :, :] * eikr
                # pol=1
                summand_x += (np.cos(agrid) * kz / pwe.k * pwe.coefficients[1, :, :])[None, :, :] * eikr
                summand_y += (np.sin(agrid) * kz / pwe.k * pwe.coefficients[1, :, :])[None, :, :] * eikr
                summand_z += (-kpgrid / pwe.k * pwe.coefficients[1, :, :])[None, :, :] * eikr
                
                em_x_flat[chunk_idcs] = np.sum(np.sum(summand_x, axis=2), axis=1)
                em_y_flat[chunk_idcs] = np.sum(np.sum(summand_y, axis=2), axis=1)
                em_z_flat[chunk_idcs] = np.sum(np.sum(summand_z, axis=2), axis=1)        
            elif field_type == 'magnetic':
                 # pol=0
                summand_x += (- kz * np.cos(agrid) * pwe.coefficients[0, :, :])[None, :, :] * eikr
                summand_y += (- kz * np.sin(agrid) * pwe.coefficients[0, :, :])[None, :, :] * eikr
                summand_z += (kpgrid * pwe.coefficients[0, :, :])[None, :, :] * eikr
                # pol=1
                summand_x += (- np.sin(agrid) * pwe.k * pwe.coefficients[1, :, :])[None, :, :] * eikr
                summand_y += (np.cos(agrid) * pwe.k * pwe.coefficients[1, :, :])[None, :, :] * eikr 
                
                em_x_flat[chunk_idcs] = 1 / omega * np.sum(np.sum(summand_x, axis=2), axis=1)
                em_y_flat[chunk_idcs] = 1 / omega * np.sum(np.sum(summand_y, axis=2), axis=1)
                em_z_flat[chunk_idcs] = 1 / omega * np.sum(np.sum(summand_z, axis=2), axis=1)                
            else:
                assert False, 'Demanded field type must be electric or magnetic!'
                    
        emx[pwe.valid(x, y, z)] = em_x_flat.reshape(xr.shape)
        emy[pwe.valid(x, y, z)] = em_y_flat.reshape(xr.shape)
        emz[pwe.valid(x, y, z)] = em_z_flat.reshape(xr.shape)
    
    return emx, emy, emz


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
    
    intens = (k / (2 * omega) * np.cos(polar_angles)[None, :, None] \
              * abs(plane_wave_expansion.coefficients[:, :len(polar_angles), :]) ** 2)
        
    srt_idcs = np.argsort(polar_angles)  # reversing order in case of downgoing
    ff = FarField(polar_angles=polar_angles[srt_idcs], azimuthal_angles=azimuthal_angles)
    ff.signal = intens[:, srt_idcs, :]
    return ff


def scattered_periodic_ff_power_per_area(far_field):
    """ Compute the scattered farfield power per area of a periodic particle arangement.
    Args:
        far_field (smuthi.postprocessing.far_field.FarField):   Farfield object of a periodic, scattered SWE
    Returns:
        Scattered farfield power per area.
    """
    return abs(np.sum(far_field.signal))


# not explicitly connected to periodic particle arrangement and therefore could 
# be a method of the class smuthi.initial_field.PlaneWave
def initial_plane_wave_power_per_area(initial_field, layer_system):
    """Compute the power carried by the initial plane wave into vertical direction per area.
    Args:
        initial_field (smuthi.initial_field.PlaneWave):     initial plane wave object
        layer_system (smuthi.layer.LayerSystem):            stratified medium
    Returns:
        Initial fields power per area.
    """

    if 0 <= initial_field.polar_angle <= np.pi/2:
        n0 = layer_system.refractive_indices[0]
    elif np.pi / 2 <= initial_field.polar_angle < np.pi:
        n0 = layer_system.refractive_indices[-1]
    k0 = 2 * np.pi * n0 / initial_field.vacuum_wavelength
    omega = initial_field.angular_frequency()
    return  k0 / (2 * omega) * np.cos(initial_field.polar_angle) * initial_field.amplitude ** 2



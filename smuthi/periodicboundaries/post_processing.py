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
from tqdm import tqdm
import sys


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
                                                                    scattered field of a periodic particle arrangement 
        pwe_down (smuthi.fields.expansions.PlaneWaveExpansion):     downwards propagating plane wave expansion of the
                                                                    scattered field of a periodic particle arrangement
    """
    
    # prepare kparallel and alpha discretization of the PWE
    i_swe = layer_system.layer_number(particle.position[2]) 
    if initial_field.polar_angle < np.pi:
        pfe = initial_field.piecewise_field_expansion(layer_system).expansion_list[2 * i_swe]
    else:
        pfe = initial_field.piecewise_field_expansion(layer_system).expansion_list[2 * i_swe + 1]
    k = pfe.k
    k0t = np.array([pfe.k_parallel[0] * np.cos(pfe.azimuthal_angles),
                    pfe.k_parallel[0] * np.sin(pfe.azimuthal_angles)])
    
    A = np.linalg.norm(np.cross(a1, a2))    
    k_parallel, azimuthal_angles, kpar_alpha_comb, gmax = pbcoup.discrete_scattering_angles(k, k0t, pfe.azimuthal_angles[0], a1, a2)
    reference_point = np.array(particle.scattered_field.reference_point)
    
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
    reference_point = [0, 0, layer_system.reference_z(i_swe)]
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
    
    pwe_up.reference_point = reference_point
    pwe_up.coefficients = pwe_up.coefficients * ejkrSiS_up[None, :, :]
    pwe_down.reference_point = reference_point
    pwe_down.coefficients = pwe_down.coefficients * ejkrSiS_down[None, :, :]
    
    return pwe_up, pwe_down    


def periodic_scattered_field_piecewise_expansion(initial_field, particle_list, layer_system, a1, a2,
                                                 layer_numbers=None):
    """Compute a piecewise field expansion of the scattered field of a periodic particle arrangement.
    Args:
        initial_field (smuthi.initial_field.PlaneWave):     initial plane wave object
        particle_list (list):                               list of smuthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem):           stratified medium
        a1 (numpy.ndarray):                                 lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):                                 lattice vector 2 in Carthesian coordinates 
        layer_numbers (list):                               if specified, append only plane wave expansions for these layers
    Returns:
        periodic scattered field as smuthi.field_expansion.PiecewiseFieldExpansion object

    """
    if layer_numbers is None:
        layer_numbers = range(layer_system.number_of_layers())

    sfld = fldex.PiecewiseFieldExpansion()
    
    # direct field ---------------------------------------------------------------------------------------------
    for particle in particle_list:
        pwe_up, pwe_down = periodic_swe_to_pwe_conversion(initial_field, layer_system, particle, a1, a2)
        pwe_up.validity_conditions.append(particle.is_outside)
        pwe_down.validity_conditions.append(particle.is_outside)
        sfld.expansion_list.append(pwe_up)
        sfld.expansion_list.append(pwe_down)

    # layer mediated scattered field --------------------------------------------------------------------------- 
    for i in layer_numbers:        
        for ip in tqdm(range(len(particle_list)), desc='Scatt. field expansion (%i)'%i, file=sys.stdout,
                                        bar_format='{l_bar}{bar}| elapsed: {elapsed} ' 'remaining: {remaining}'):
            
            pwe_up_exc = sfld.expansion_list[2 * ip]
            pwe_down_exc =  sfld.expansion_list[2 * ip + 1]
            i_sca = layer_system.layer_number(pwe_up_exc.lower_z) 
            if ip == 0:
                pwe_up, pwe_down = layer_system.response((pwe_up_exc, pwe_down_exc), i_sca, i)
            else:
                add_up, add_down = layer_system.response((pwe_up_exc, pwe_down_exc), i_sca, i)
                
                pwe_up = pwe_up + add_up
                pwe_down = pwe_down + add_down
                
            pwe_up.validity_conditions.append(particle.is_outside)
            pwe_down.validity_conditions.append(particle.is_outside)
            
        # in bottom_layer, suppress upgoing waves, and in top layer, suppress downgoing waves
        if i > 0:
            sfld.expansion_list.append(pwe_up)
        if i < layer_system.number_of_layers()-1:
            sfld.expansion_list.append(pwe_down)

    return sfld


def total_periodic_field_piecewise_expansion(initial_field, particle_list, layer_system, a1, a2,
                                             layer_numbers=None):
    """Compute a piecewise field expansion of the total field of a periodic particle arrangement.
    Args:
        initial_field (smuthi.initial_field.PlaneWave):     initial plane wave object
        particle_list (list):                               list of smuthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem):           stratified medium
        a1 (numpy.ndarray):                                 lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):                                 lattice vector 2 in Carthesian coordinates 
        layer_numbers (list):                               if specified, append only plane wave expansions for these layers
    Returns:
        periodic total field as smuthi.fields.expansions.PiecewiseFieldExpansion object
    """
    if layer_numbers is None:
        layer_numbers = range(layer_system.number_of_layers())
        
    # piecewise field expansion of the initial field 
    pfe_init = initial_field.piecewise_field_expansion(layer_system)
    
    # piecewise field expansion of the periodic scattered fields
    pfe_sca = periodic_scattered_field_piecewise_expansion(initial_field=initial_field,
                            particle_list=particle_list, layer_system=layer_system,
                            a1=a1, a2=a2, layer_numbers=layer_numbers)
    
    # total field
    pfe_total = copy.deepcopy(pfe_sca)
    idx_a = np.argwhere(pfe_sca.expansion_list[0].azimuthal_angles == pfe_init.expansion_list[0].azimuthal_angles)[0][0]
    idx_kpar = np.argwhere(pfe_sca.expansion_list[0].k_parallel == pfe_init.expansion_list[0].k_parallel)[0][0]
    for i in layer_numbers:   
        exp_idx = 2 * len(particle_list) + 2 * i
        if i > 0:
            pfe_total.expansion_list[exp_idx - 1].coefficients[:, idx_kpar, idx_a] += pfe_init.expansion_list[2 * i].coefficients.reshape(2)
        if i < layer_system.number_of_layers()-1:
            pfe_total.expansion_list[exp_idx].coefficients[:, idx_kpar, idx_a] += pfe_init.expansion_list[2 * i + 1].coefficients.reshape(2)
    
    # add initial field contributions from outside
    pfe_total.expansion_list.insert(2 * len(particle_list), pfe_init.expansion_list[0])
    pfe_total.expansion_list.append(pfe_init.expansion_list[-1])
            
    return pfe_total


def transmitted_plane_wave_expansion(initial_field, particle_list, layer_system, a1, a2):
    """Compute the plane wave expansion of the total transmitted field.
    Args:
        initial_field (smuthi.initial_field.PlaneWave):     initial plane wave object
        particle_list (list):                               list of smuthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem):           stratified medium
        a1 (numpy.ndarray):                                 lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):                                 lattice vector 2 in Carthesian coordinates 
    Returns:
        smuthi.fields.expansions.PlaneWaveExpansion object of the total transmitted field
    """
    if 0 <= initial_field.polar_angle < np.pi / 2: 
        i_T = layer_system.number_of_layers() - 1
        pwe_init_T = initial_field.plane_wave_expansion(layer_system, i_T)[0]
    else:
        i_T = 0
        pwe_init_T = initial_field.plane_wave_expansion(layer_system, i_T)[1]
    
    pfe_sca = periodic_scattered_field_piecewise_expansion(initial_field, particle_list, layer_system, a1, a2, layer_numbers=[i_T])
        
    # layer response of all scattered fields
    pwe_T = copy.deepcopy(pfe_sca.expansion_list[-1])
    idx_a = np.argwhere(pwe_T.azimuthal_angles == pwe_init_T.azimuthal_angles)[0][0]
    idx_kpar = np.argwhere(pwe_T.k_parallel == pwe_init_T.k_parallel)[0][0]
    
    # add the initial field contribution
    pwe_T.coefficients[:, idx_kpar, idx_a] += pwe_init_T.coefficients.reshape(2)
    
    # if particles are located in the outer layer, add the direct scattered field contributions
    i_sca = layer_system.layer_number(particle_list[0].position[2]) # all particles have to be in the same layer
    if i_sca == i_T:
        for ip in range(len(particle_list)):
            if i_T == 0:
                pwe_T = pwe_T + pfe_sca.expansion_list[2 * ip + 1]
            else:
                pwe_T = pwe_T + pfe_sca.expansion_list[2 * ip]
 
    return pwe_T


def reflected_plane_wave_expansion(initial_field, particle_list, layer_system, a1, a2):
    """Compute the plane wave expansion of the total reflected field.
    Args:
        initial_field (smuthi.initial_field.PlaneWave):     initial plane wave object
        particle_list (list):                               list of smuthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem):           stratified medium
        a1 (numpy.ndarray):                                 lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):                                 lattice vector 2 in Carthesian coordinates 
    Returns:
        smuthi.fields.expansions.PlaneWaveExpansion object of the total reflected field
    """
    if 0 <= initial_field.polar_angle < np.pi / 2: 
        i_R = 0
        pwe_init_R = initial_field.plane_wave_expansion(layer_system, i_R)[1]
    else:
        i_R = layer_system.number_of_layers() - 1
        pwe_init_R = initial_field.plane_wave_expansion(layer_system, i_R)[0]
    
    pfe_sca = periodic_scattered_field_piecewise_expansion(initial_field, particle_list, layer_system, a1, a2, layer_numbers=[i_R])
    
    # layer response of all scattered fields
    pwe_R = copy.deepcopy(pfe_sca.expansion_list[-1])
    idx_a = np.argwhere(pwe_R.azimuthal_angles == pwe_init_R.azimuthal_angles)[0][0]
    idx_kpar = np.argwhere(pwe_R.k_parallel == pwe_init_R.k_parallel)[0][0]
    
    # add the initial field contribution
    pwe_R.coefficients[:, idx_kpar, idx_a] += pwe_init_R.coefficients.reshape(2)
    
    # if particles are located in the outer layer, add the direct scattered field contributions
    i_sca = layer_system.layer_number(particle_list[0].position[2]) # all particles have to be in the same layer
    if i_sca == i_R:
        for ip in range(len(particle_list)):
            if i_R == 0:
                pwe_R = pwe_R + pfe_sca.expansion_list[2 * ip + 1]
            else:
                pwe_R = pwe_R + pfe_sca.expansion_list[2 * ip]
 
    return pwe_R


# Field evaluation follows the procedure of Smuthi's "old (legacy)" electric field and magnetic field routines 
# basically identical to the field evaluation of any plane wave, but as a sum of discrete propagation angles (instead of an integral)
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
        Tuple of (E_x, E_y, E_z) or (H_x, H_y, H_z) numpy.ndarray objects with the Cartesian coordinates of
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
        
        for i_chunk in range(int(np.ceil(xr.size / chunksize))):
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


def conjugated_poynting_vector(E, H):
    """ 
    Args:
        E (tuple):  electric field tuple (E_x, E_y, E_z)
        H (tuple):  magnetic field tuple (H_x, H_y, H_z)
    Returns:
        Returns the complex conjugated Poynting vector S*.
    """
    ex_conj = np.conjugate(E[0])
    ey_conj = np.conjugate(E[1])
    ez_conj = np.conjugate(E[2])
    
    hx = H[0]
    hy = H[1]
    hz = H[2]
    
    sx = 1 / 2 * (ey_conj * hz - ez_conj * hy)
    sy = 1 / 2 * (ez_conj * hx - ex_conj * hz)
    sz = 1 / 2 * (ex_conj * hy - ey_conj * hx)
    
    return sx, sy, sz


# Note: The FarField.signal corresponds to the plane wave's power per interface area (Theobald 2021 dissertation, eq.(6.29))
# It differs from the smuthi.postprocessing.far_field.FarField object that is constructed from a PlaneWaveExpansion
# with a continuous angular distribution.
# Introduces the FarField.signal_type "normalized power".
def periodic_pwe_to_ff_conversion(plane_wave_expansion, initial_field, layer_system):
    """Compute the far field of a plane wave expansion object of discrete plane waves.
    Args:
        plane_wave_expansion (smuthi.fields.expansions.PlaneWaveExpansion):     
                                    plane wave expansion of the scattered field of a periodic particle
                                    arrangement that is to be convert into a far field object.
        initial_field (smuthi.initial_field.PlaneWave):     initial plane wave object
        layer_system (smuthi.layer.LayerSystem):            stratified medium
    Returns:
        A smuthi.postprocessing.far_field.FarField object.
    """
    omega = flds.angular_frequency(initial_field.vacuum_wavelength)
    k = plane_wave_expansion.k
    kp = plane_wave_expansion.k_parallel
    if plane_wave_expansion.kind == 'upgoing':
        polar_angles = np.arcsin(kp[kp < k] / k)
    elif plane_wave_expansion.kind == 'downgoing':
        polar_angles = np.pi - np.arcsin((kp / k)[kp / k < 1])
    else:
        raise ValueError('PWE type not specified')
    if any(polar_angles.imag):
        raise ValueError('complex angles are not allowed')
    polar_angles = polar_angles[np.logical_not(np.isnan(polar_angles))]
    azimuthal_angles = plane_wave_expansion.azimuthal_angles
    
    # power per interface area (Theobald 2021 dissertation, eq.(6.29))
    normalized_power = (k / (2 * omega) * np.cos(polar_angles)[None, :, None] \
              * abs(plane_wave_expansion.coefficients[:, :len(polar_angles), :]) ** 2)
        
    srt_idcs = np.argsort(polar_angles)  # reversing order in case of downgoing
    ff = FarField(polar_angles=polar_angles[srt_idcs], azimuthal_angles=azimuthal_angles,
                  signal_type='normalized power')
    ff.signal = normalized_power[:, srt_idcs, :]
    return ff


def scattered_periodic_ff_power_per_area(far_field):
    """ Compute the scattered far field power per area of a periodic particle arangement.
    Args:
        far_field (smuthi.postprocessing.far_field.FarField):   Farfield object of a periodic, scattered SWE
    Returns:
        Scattered farfield power per area.
    """
    if not far_field.signal_type == 'normalized power':
        raise TypeError('Far field does not consist of discrete plane waves.')
    return abs(np.sum(far_field.signal))


# not explicitly connected to periodic particle arrangement and therefore could be a method of the class smuthi.initial_field.PlaneWave
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







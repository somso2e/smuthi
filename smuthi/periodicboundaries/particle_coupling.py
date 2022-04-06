"""This module contains functions to evaluate the coupling between a 
single particle and a periodic particle arrangement."""

import smuthi.periodicboundaries.ewald_lattice_sums as pbels
import smuthi.periodicboundaries.ewald_helper as pbeh
import smuthi.periodicboundaries.coupling_helper as pbch
import smuthi.fields as flds
import smuthi.utility.math as sf
import smuthi.linearsystem.particlecoupling.layer_mediated_coupling as laycoup
from numba import njit, prange
import numpy as np

@njit()
def discrete_scattering_angles(k, k0t, azimuthal_angle_init, a1, a2):
    """Evaluates the discrete scattering angles of a plane wave exciting a periodic
    particle arrangement. 
    In particular, combinations of in-plane wavenumber and azimuthal angle. 
    Args:
        k (float):                      initial field's wavenumber
        k0t (numpy.ndarray):            complex in-plane wave vector in Carthesian coordinates 
        azimuthal_angle_init (float):   azimuthal angle of the initial field
        a1 (numpy.ndarray):             lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):             lattice vector 2 in Carthesian coordinates
    Returns:
        k_parallel (numpy.ndarray):         in-plane wavenumbers
        azimuthal_angles (numpy.ndarray):   azimuthal angles
        kpar_alpha_comb (list):             list of length len(k_parallel) that connects indices 
                                            of k_parallel with indices of azimuthal_angles 
        gmax (int):                         maximal order of unit cells to which the reciprocal
                                            lattive vector g is evaluated (k = k_inc - g) 
    """
    b1, b2 = pbeh.reciprocal_lattice_vec(a1, a2) 
    
     # find all kpar / alpha combinations and form a matrix
    kpar_list = np.array([np.linalg.norm(k0t)], np.float64)
    alpha_list = np.array([azimuthal_angle_init], np.float64)
    kpar_alpha_list = [np.array([0], np.int64)]
      
    flag, num_evan = True, 0
    gmax = 0
    # ensures that all propagating scattering orders + 3 evanescent orders are taken into account
    # 3 evanescent orders are somewhat arbitrary (could be set from outside)
    while flag == True or num_evan < 3:
        gmax += 1
        n1, n2 = pbeh.n1n2_indices(gmax)
        g = n1.reshape(n1.shape[0], 1) * b1.reshape(1, 2) + n2.reshape(n2.shape[0], 1) * b2.reshape(1, 2)
        ksca = k0t.reshape(2) - g    
        
        kpar = np.zeros(ksca.shape[0], np.float64)
        alpha = np.zeros(ksca.shape[0], np.float64)
        for idx in range(ksca.shape[0]):
            kpar[idx] = np.linalg.norm(ksca[idx, :])
            alpha[idx] = np.arctan2(ksca[idx, 1], ksca[idx, 0]) 
        for idx, ang in enumerate(alpha):
            if ang < 0:
                alpha[idx] = 2 * np.pi + ang
            
        if np.min(kpar) > k.real: # all scattered waves are evanescent
            flag = False
            num_evan += 1
        
        for idx in range(kpar.shape[0]):     
            isclose_a = np.zeros(alpha_list.shape[0], np.int32)
            for ii in range(alpha_list.shape[0]):
                isclose_a[ii] = np.round(alpha[idx], 5) == np.round(alpha_list[ii], 5)
            if not isclose_a.any():
                alpha_list = np.append(alpha_list, alpha[idx])
                alpha_idx = alpha_list.shape[0] - 1
            else:
                alpha_idx = np.where(isclose_a == 1)[0][0] 
                
            isclose_k = np.zeros(kpar_list.shape[0], np.int32)
            for ii in range(kpar_list.shape[0]):
                isclose_k[ii] = np.round(kpar[idx], 5) == np.round(kpar_list[ii], 5)
            if not isclose_k.any():
                kpar_list = np.append(kpar_list, kpar[idx])
                kpar_alpha_list.append(np.array([alpha_idx], np.int64))  
            else:
                kpar_idx = np.where(isclose_k == 1)[0][0]
                kpar_alpha_list[kpar_idx] = np.append(kpar_alpha_list[kpar_idx], alpha_idx)
       
    k_parallel = np.sort(kpar_list)
    azimuthal_angles = np.sort(alpha_list)
    
    for idx, kpar in enumerate(k_parallel):
        kpar_idx = np.where(kpar == kpar_list)[0][0]
        for idx_2, alpha_idx in enumerate(kpar_alpha_list[kpar_idx]):
            alpha = alpha_list[alpha_idx]
            alpha_idx_new = np.where(alpha == azimuthal_angles)[0][0]           
            if idx_2 == 0:
                arr = [alpha_idx_new]
            else:
                arr.append(alpha_idx_new)
        arr.sort()        
        if idx == 0:
            kpar_alpha_comb = [arr]
        else:
            kpar_alpha_comb.append(arr)

    return k_parallel, azimuthal_angles, kpar_alpha_comb, gmax


@njit(parallel=True)
def periodic_coupling_matrix(vacuum_wavelength, k0t, azimuthal_angle, layer_thicknesses, layer_refractive_indices,
                            i_sca, positions, radii, lmax_array, mmax_array, a1, a2, eta, a5b5_mat, mmax_global):
    """ Coupling matrix between one particle and a periodic particle arrangement.
        "prange()" enables the multithreading of the Numba code. 
    Args:
        vacuum_wavelength (float):          vacuum wavelength in length unit   
        k0t (numpy.ndarray):                complex in-plane wave vector in Carthesian coordinates  
        azimuthal_angle (numpy.ndarray):    azimuthal angle of excitation 
        layer_thicknesses (numpy.ndarray):  layer thicknesses of stratified medium
        layer_refractive_indices (numpy.ndarray): layer refractive indices of stratified medium
        i_sca (int):                        layer number that contains all scattering particles
        positions (numpy.ndarray):          position of each particle
        radii (numpy.ndarray):              radius of each particle
        lmax_array (numpy.ndarray):         lmax of each particle 
        mmax_array (numpy.ndarray):         mmax of each particle
        a1 (numpy.ndarray):                 lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):                 lattice vector 2 in Carthesian coordinates
        eta (float):                        Ewald sum separation parameter
        a5b5_mat (numpy.ndarray):           a5, b5 coefficient lookup table
        mmax_global (int):                  global maximal multipole order of all particles 
    Returns:
        sum of direct and layer mediated coupling matrices
    """    
    s_max = positions.shape[0]
    shape = pbch.index_block(s_max - 1, lmax_array, mmax_array)[1]    
    coup_mat = np.zeros((shape, shape), np.complex128)
    
    # alternative to status bar, to access the progress of the prange-numba-loop
    progress = np.zeros((s_max, s_max), np.int64)
    sections = np.zeros(10, np.float64) # in 10% steps
    print('Coupling matrix computation finished: 0.0 %.')
    
    for s1 in prange(s_max):
        idx1 = pbch.index_block(s1, lmax_array, mmax_array)
        for s2 in prange(s_max):
            idx2 = pbch.index_block(s2, lmax_array, mmax_array)
            coup_mat[idx1[0]:idx1[1], idx2[0]:idx2[1]] = (
                                    periodic_direct_coupling_block(vacuum_wavelength, layer_refractive_indices[i_sca],
                                            k0t, azimuthal_angle, positions[s1], lmax_array[s1], mmax_array[s1], radii[s1],
                                            positions[s2], lmax_array[s2], mmax_array[s2], radii[s2], a1, a2, eta, a5b5_mat, mmax_global)
                                    + periodic_layer_mediated_coupling_block(vacuum_wavelength, layer_thicknesses,
                                              layer_refractive_indices, i_sca, k0t, azimuthal_angle, positions[s1], lmax_array[s1],
                                              mmax_array[s1], radii[s1], positions[s2], lmax_array[s2], mmax_array[s2],
                                              radii[s2], a1, a2))
            progress[s1, s2] = 1 
            if (np.sum(progress) / s_max ** 2) >= np.sum(sections) + 0.1:
                sections[:int(np.floor(np.sum(progress) / s_max ** 2 * 1e1))] = 0.1
                print('Coupling matrix computation finished:', np.round(np.sum(sections) * 100), '%.')
                
    return coup_mat



@njit()
def periodic_direct_coupling_block(vacuum_wavelength, refractive_index, k0t, azimuthal_angle, 
                                   rs1, l_max1, m_max1, radius1, rs2, l_max2, m_max2, radius2,
                                   a1, a2, eta, a5b5_lookup, mmax_global):
    """  Direct particle coupling matrix :math:`W^R` between one (receiving) particle and
        a (emitting) periodic particle arrangement.
        Based on a intermediate plane wave expansion or an Ewald lattice sum. 
    Args:
        vacuum_wavelength (float):          vacuum wavelength in length unit   
        refractive_index (complex):         refractive index of the scattering layer 
        k0t (numpy.ndarray):                complex in-plane wave vector in Carthesian coordinates  
        azimuthal_angle (numpy.ndarray):    azimuthal angle of excitation 
        rs1 (numpy.ndarray):                position of the receiving particle in Carthesian coordinates 
        l_max1 (int):                       maximal multipole degree of the receiving particle
        m_max1 (int):                       maximal multipole order of the receiving particle
        radius1 (float):                    radius of the circumscribing sphere of the receiving particle
        rs2 (numpy.ndarray):                position of the emitting particle in Carthesian coordinates 
        l_max2 (int):                       maximal multipole degree of the emitting particle
        m_max2 (int):                       maximal multipole order of the emitting particle
        radius2 (float):                    radius of the circumscribing sphere of the emitting particle
        a1 (numpy.ndarray):                 lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):                 lattice vector 2 in Carthesian coordinates
        eta (float):                        Ewald sum separation parameter
        a5b5_lookup (numpy.ndarray):        lookup table for a5, b5 coefficients
        mmax_global (int):                  global maximal multipole order of all particles 
    
    Returns:
        w (numpy.ndarray):                  direct coupling block 
    """
    blocksize1 = flds.blocksize(l_max1, m_max1)
    blocksize2 = flds.blocksize(l_max2, m_max2)
    w = np.zeros((blocksize1, blocksize2), np.complex128)
    
    k = complex(2 * np.pi * refractive_index / vacuum_wavelength)
    
    # intermediate plane wave expansion
    if np.abs(rs1[2] - rs2[2]) > radius1 + radius2: # no vertical overlap of particles' circumscribing spheres
        k_parallel, azimuthal_angles, kpar_alpha_comb, gmax = discrete_scattering_angles(k, k0t, azimuthal_angle, a1, a2) 
        
        rs2s1 = rs1 - rs2
        rhos2s1 = np.linalg.norm(rs2s1[0:2])
        phis2s1 = np.arctan2(rs2s1[1], rs2s1[0])
        
        kz = np.sqrt(k ** 2 - k_parallel ** 2 + 0j)
        kz = (kz.imag >= 0) * kz + (kz.imag < 0) * (-kz)  # Branch cut to prohibit negative imaginary part
        
        # transformation coefficients; either upwards or downwards
        B = [np.zeros((2, blocksize1, k_parallel.shape[0]), np.complex128),
             np.zeros((2, blocksize2, k_parallel.shape[0]), np.complex128)]
        # list index: particle, np indices: pol, n, kpar_idx
        
        m_vec = [np.zeros(blocksize1, np.int32), np.zeros(blocksize2, np.int32)]
        
        # precompute spherical functions
        ct = kz / k
        st = k_parallel / k
        if rs2s1[2] > 0:
            _, pilm, taulm = sf.legendre_normalized_numbed(ct, st, l_max1)
        else:
            _, pilm, taulm = sf.legendre_normalized_numbed(-ct, st, l_max1)
            
        for tau in range(2):
            for m in range(-m_max1, m_max1 + 1):
                for l in range(max(1, abs(m)), l_max1 + 1):
                    n = flds.multi_to_single_index(tau, l, m, l_max1, m_max1)
                    m_vec[0][n] = m
                    for pol in range(2):
                        B[0][pol, n, :] = pbch.transformation_coefficients_vwf(tau, l, m, pol, kp=k_parallel, kz=kz,
                                                                          pilm=pilm, taulm=taulm, dagger=True)                  
        if rs2s1[2] > 0:
            _, pilm, taulm = sf.legendre_normalized_numbed(ct, st, l_max2)
        else:
            _, pilm, taulm = sf.legendre_normalized_numbed(-ct, st, l_max2)
                         
        for tau in range(2):
            for m in range(-m_max2, m_max2 + 1):
                for l in range(max(1, abs(m)), l_max2 + 1):
                    n = flds.multi_to_single_index(tau, l, m, l_max2, m_max2)
                    m_vec[1][n] = m
                    for pol in range(2):
                        B[1][pol, n, :] = pbch.transformation_coefficients_vwf(tau, l, m, pol, kp=k_parallel, kz=kz,
                                                                          pilm=pilm, taulm=taulm, dagger=False)
                        
        # k_parallel (and idx), azimuthal_angles arrays that reduce the summation to non-zero elements 
        alpha_vec = np.zeros((2 * gmax + 1) ** 2, np.float64)
        kpar_vec = np.zeros((2 * gmax + 1) ** 2, np.complex128)
        k_idx = np.zeros((2 * gmax + 1) ** 2, np.int32)
        s = 0
        for idx in range(k_parallel.shape[0]):
            num = len(kpar_alpha_comb[idx])
            alpha_vec[s:s+num] = azimuthal_angles[np.array(kpar_alpha_comb[idx], np.int32)]
            kpar_vec[s:s+num] = k_parallel[idx]
            k_idx[s:s+num] = idx
            s += num
            
        kz_vec = np.sqrt(k ** 2 - kpar_vec ** 2 + 0j)
        kz_vec = (kz_vec.imag >= 0) * kz_vec + (kz_vec.imag < 0) * (-kz_vec)  # Branch cut  to prohibit negative imaginary part
        jacobi_vector = 1 / (kz_vec * k)
        if rs2s1[2] < 0:
            kz_vec = -kz_vec
            
        prefac = 8 * np.pi / np.linalg.norm(np.cross(a1, a2))
        m2_minus_m1 = m_vec[1] - m_vec[0].reshape(m_vec[0].shape[0], 1)
        eikrs2s1= np.exp(1j * (kpar_vec * rhos2s1 * np.cos(alpha_vec - phis2s1) + kz_vec * rs2s1[2]))
        
        for n1 in range(blocksize1):
            for n2 in range(blocksize2):
                BB = np.zeros(k_parallel.shape[0], np.complex128) 
                eimpma = np.exp(1j * m2_minus_m1[n1, n2] * alpha_vec)
                for pol in range(2):
                    BB += B[0][pol, n1, :] * B[1][pol, n2, :]
                w[n1, n2] = np.sum(eimpma * eikrs2s1 * jacobi_vector * BB[k_idx])               
        w = prefac * w 
    
    # Ewald lattice sum 
    else:
        c = rs1 - rs2 
        for l1 in range(1, l_max1 + 1):
            for l2 in range(1, l_max2 + 1):
                for m1 in range(-min(m_max1, l1), min(m_max1, l1) + 1):
                    for m2 in range(-min(m_max2, l2), min(m_max2, l2) + 1):
                        M = m2 - m1
                        prefac = pbeh.normalization_spherical_harmonics(M)
                        A, B = complex(0), complex(0) 
                        for L in range(max(abs(l1 - l2), abs(M)), l1 + l2 + 1): # if L < abs(M) then P=0 
                            Dlm = pbels.D_LM(L, M, k, k0t, a1, a2, eta, c)
                            A += a5b5_lookup[0, l1 - 1, l2 - 1, m1 + mmax_global, m2 + mmax_global, L] * Dlm
                            B += a5b5_lookup[1, l1 - 1, l2 - 1, m1 + mmax_global, m2 + mmax_global, L] * Dlm
                        for tau1 in range(2):
                            n1 = flds.multi_to_single_index(tau1, l1, m1, l_max1, m_max1)
                            for tau2 in range(2):
                                n2 = flds.multi_to_single_index(tau2, l2, m2, l_max2, m_max2)
                                if tau1 == tau2:
                                    w[n1, n2] = prefac * A
                                else:
                                    w[n1, n2] = prefac * B 
    return w


@njit()
def periodic_layer_mediated_coupling_block(vacuum_wavelength, layer_thicknesses, layer_refractive_indices,
                                           i_sca, k0t, azimuthal_angle, rs1, l_max1, m_max1, radius1,
                                           rs2, l_max2, m_max2, radius2, a1, a2):
    """ Layer-system mediated particle coupling matrix :math:`W^R` between one (receiving) particle
        and a (emitting) periodic particle arrangement.
    Args:
        vacuum_wavelength (float):          vacuum wavelength in length unit
        layer_thicknesses (numpy.ndarray):  layer thicknesses of stratified medium
        layer_refractive_indices (numpy.ndarray): layer refractive indices of stratified medium
        i_sca (int):                        layer number that contains all scattering particles
        k0t (numpy.ndarray):                complex in-plane wave vector in Carthesian coordinates  
        azimuthal_angle (numpy.ndarray):    azimuthal angle of excitation 
        rs1 (numpy.ndarray):                position of the receiving particle in Carthesian coordinates 
        l_max1 (int):                       maximal multipole degree of the receiving particle
        m_max1 (int):                       maximal multipole order of the receiving particle
        radius1 (float):                    radius of the circumscribing sphere of the receiving particle
        rs2 (numpy.ndarray):                position of the emitting particle in Carthesian coordinates 
        l_max2 (int):                       maximal multipole degree of the emitting particle
        m_max2 (int):                       maximal multipole order of the emitting particle
        radius2 (float):                    radius of the circumscribing sphere of the emitting particle
        a1 (numpy.ndarray):                 lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):                 lattice vector 2 in Carthesian coordinates
        a5b5_mat (numpy.ndarray):           a5, b5 coefficient lookup table
        mmax_global (int):                  global maximal multipole order of all particles 
    
    Returns:
        wr (numpy.ndarray):                 layer mediated coupling block 
    """
    blocksize1 = flds.blocksize(l_max1, m_max1)
    blocksize2 = flds.blocksize(l_max2, m_max2)
    wr = np.zeros((blocksize1, blocksize2), np.complex128)
    
    omega = 2 * np.pi / vacuum_wavelength
    k = complex(omega * layer_refractive_indices[i_sca])
    k_parallel, azimuthal_angles, kpar_alpha_comb, gmax = discrete_scattering_angles(k, k0t, azimuthal_angle, a1, a2)
    kz = np.sqrt(k ** 2 - k_parallel ** 2 + 0j)
    kz = (kz.imag >= 0) * kz + (kz.imag < 0) * (-kz)  # Branch cut to prohibit negative imaginary part
    
    rs2s1 = rs1 - rs2
    rhos2s1 = np.linalg.norm(rs2s1[0:2])
    phis2s1 = np.arctan2(rs2s1[1], rs2s1[0])
    ziss1 = rs1[2] - np.sum(layer_thicknesses[:i_sca])
    ziss2 = rs2[2] - np.sum(layer_thicknesses[:i_sca])
    
    # phase factors
    ejkz = np.zeros((2, 2, k_parallel.shape[0]), np.complex128)  # indices are: particle, plus/minus, kpar_idx
    ejkz[0, 0, :] = np.exp(1j * kz * ziss1)
    ejkz[0, 1, :] = np.exp(- 1j * kz * ziss1)
    ejkz[1, 0, :] = np.exp(1j * kz * ziss2)
    ejkz[1, 1, :] = np.exp(- 1j * kz * ziss2)
    
    # remove inf from ejkz
    # (otherwise evanescent waves cause overflow in case of very thick layers)
    idces =  np.where(np.isinf(ejkz[0, 1, :]))[0]
    for ii in idces:
        ejkz[0, 1, ii] = 1.79769313e+308
    idces =  np.where(np.isinf(ejkz[1, 1, :]))[0]
    for ii in idces:
        ejkz[1, 1, ii] = 1.79769313e+308
    
    # layer response
    L = np.zeros((2, 2, 2, k_parallel.shape[0]), np.complex128)  # polarization, pl/mn1, pl/mn2, kpar_idx
    for pol in range(2):
        L[pol, :, :, :] = pbch.layersystem_response_matrix(pol, layer_thicknesses, layer_refractive_indices,
                                                      k_parallel, omega, i_sca, i_sca)
        
    # transformation coefficients; either upwards or downwards
    B = [np.zeros((2, 2, blocksize1, k_parallel.shape[0]), np.complex128),
         np.zeros((2, 2, blocksize2, k_parallel.shape[0]), np.complex128)]
    # list index: particle, np indices: pol, plus/minus, n, kpar_idx
    
    m_vec = [np.zeros(blocksize1, np.int32), np.zeros(blocksize2, np.int32)]
    
    # precompute spherical functions
    ct = kz / k
    st = k_parallel / k    
    _, pilm_list_pl, taulm_list_pl = sf.legendre_normalized_numbed(ct, st, l_max1)
    _, pilm_list_mn, taulm_list_mn = sf.legendre_normalized_numbed(-ct, st, l_max1)
    pilm = (pilm_list_pl, pilm_list_mn)
    taulm = (taulm_list_pl, taulm_list_mn)
    
    for tau in range(2):
        for m in range(-m_max1, m_max1 + 1):
            for l in range(max(1, abs(m)), l_max1 + 1):
                n = flds.multi_to_single_index(tau, l, m, l_max1, m_max1)
                m_vec[0][n] = m
                for iplmn in range(2):
                    for pol in range(2):
                        B[0][pol, iplmn, n, :] = pbch.transformation_coefficients_vwf(tau, l, m, pol, kp=k_parallel, kz=kz,
                                                                                 pilm=pilm[iplmn],
                                                                                 taulm=taulm[iplmn], dagger=True)

    _, pilm_list_pl, taulm_list_pl = sf.legendre_normalized_numbed(ct, st, l_max2)
    _, pilm_list_mn, taulm_list_mn = sf.legendre_normalized_numbed(-ct, st, l_max2)
    pilm = (pilm_list_pl, pilm_list_mn)
    taulm = (taulm_list_pl, taulm_list_mn)
                       
    for tau in range(2):
        for m in range(-m_max2, m_max2 + 1):
            for l in range(max(1, abs(m)), l_max2 + 1):
                n = flds.multi_to_single_index(tau, l, m, l_max2, m_max2)
                m_vec[1][n] = m
                for iplmn in range(2):
                    for pol in range(2):
                        B[1][pol, iplmn, n, :] = pbch.transformation_coefficients_vwf(tau, l, m, pol, kp=k_parallel, kz=kz,
                                                                                 pilm=pilm[iplmn],
                                                                                 taulm=taulm[iplmn], dagger=False)

    # k_parallel (and idx), azimuthal_angles arrays that reduce the summation to non-zero elements 
    alpha_vec = np.zeros((2 * gmax + 1) ** 2, np.float64)
    kpar_vec = np.zeros((2 * gmax + 1) ** 2, np.complex128)
    k_idx = np.zeros((2 * gmax + 1) ** 2, np.int32)
    s = 0
    for idx in range(k_parallel.shape[0]):
        num = len(kpar_alpha_comb[idx])
        alpha_vec[s:s+num] = azimuthal_angles[np.array(kpar_alpha_comb[idx], np.int32)]
        kpar_vec[s:s+num] = k_parallel[idx]
        k_idx[s:s+num] = idx
        s += num
        
    prefac = 8 * np.pi / np.linalg.norm(np.cross(a1, a2))
    jacobi_vector = 1 / (kz * k)
    m2_minus_m1 = m_vec[1] - m_vec[0].reshape(m_vec[0].shape[0], 1)
    eikparrs2s1 = np.exp(1j * kpar_vec * rhos2s1 * np.cos(alpha_vec - phis2s1))
    
    for n1 in range(blocksize1):
        BeL = np.zeros((2, 2, k_parallel.shape[0]), np.complex128)  # indices are: pol, plmn2, n1, kpar_idx
        for iplmn1 in range(2):
            for pol in range(2):
                BeL[pol, :, :] += (L[pol, iplmn1, :, :]
                                        * B[0][pol, iplmn1, n1, :]
                                        * ejkz[0, iplmn1, :])
        for n2 in range(blocksize2):
            eimpma = np.exp(1j * m2_minus_m1[n1, n2] * alpha_vec)
            BeLBe = np.zeros((k_parallel.shape[0]), np.complex128)
            laycoup.eval_BeLBe(BeLBe, BeL, B[1], ejkz, n2)           
            wr[n1, n2] = np.sum(eimpma * eikparrs2s1 * jacobi_vector[k_idx] * BeLBe[k_idx])  
    return prefac * wr    
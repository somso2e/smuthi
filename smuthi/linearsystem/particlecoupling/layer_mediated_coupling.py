"""This module contains functions to compute the layer mediated particle 
coupling coefficients."""

import numpy as np
import scipy.special
from numba import complex128, int64, jit
import smuthi.fields
import smuthi.fields.transformations as trf
import smuthi.layers as lay
import smuthi.utility.math as sf


@jit(complex128(complex128[:], complex128[:]),
     nopython=True, cache=True, nogil=True)
def numba_trapz(y, x):
    out = 0.0 + 0.0j
    #TODO implement some (optional) advanced summation?
    #e.g. https://github.com/nschloe/accupy/blob/master/accupy/sums.py
    #or better Sum2  from https://doi.org/10.1137/030601818 (Algorithm 4.4)
    #Note, that this may need to have exact summation for x and y, and exact product.
    for i in range( len(y) - 2 ):
        out += (x[i+1]-x[i]) * (y[i+1] + y[i])/2.0
    return out


@jit((complex128[:], complex128[:,:,:],
      complex128[:,:,:,:],complex128[:,:,:], int64),
     nopython=True,cache=True
     ,nogil=True
     # , parallel=True
    )
def eval_BeLBe(BeLBe, BeL, B1, ejkz, n2):
    for k in range(len(BeLBe)):
        for iplmn2 in range(2):
            for pol in range(2):
                BeLBe[k] += BeL[pol, iplmn2, k] * B1[pol, iplmn2, n2, k] * ejkz[1, 1 - iplmn2, k]
        

def g_function(vacuum_wavelength, receiving_particle, emitting_particle, layer_system, k_parallel):
    """
    This function returns the function g_n(k_\rho) as defined in equation (16) of the paper "A quick way to approximate
    a Sommerfeld-Weyl_type Sommerfeld integral" by Chew (1988) for the Sommerfeld integral for the layer-mediated
    particle coupling (compare Amos Egel's dissertation, equation (3.46)).
    The purpose of this function is to allow for stationary phase approximation in the case of  large lateral distances.

    The use of this function is constrained to the important special case of particles above a (possibly layered)
    substrate.

    Args:
        vacuum_wavelength (float):                          Vacuum wavelength :math:`\lambda` (length unit)
        receiving_particle (smuthi.particles.Particle):     Particle that receives the scattered field (must be in top
                                                            layer)
        emitting_particle (smuthi.particles.Particle):      Particle that emits the scattered field (must be in top
                                                            layer)
        layer_system (smuthi.layers.LayerSystem):           Stratified medium in which the coupling takes place
        k_parallel (float):                                 In-plane wavenumber

    Returns:
        a tuple with the following items:
        - g: a matrix of g_function as entries (dimension is blocksize x blocksize)
        - rho: radial distance between particles
        - delta_z: z-argument of e^(ik_z z) term
        - bessel_order: matrix of positive numbers to be used as the order of Hankel functions
    """

    # angular frequency
    omega = smuthi.fields.angular_frequency(vacuum_wavelength)

    # index specs
    blocksize1 = smuthi.fields.blocksize(receiving_particle.l_max, receiving_particle.m_max)
    blocksize2 = smuthi.fields.blocksize(emitting_particle.l_max, emitting_particle.m_max)

    # cylindrical coordinates of relative position vectors
    rs1 = np.array(receiving_particle.position)
    rs2 = np.array(emitting_particle.position)
    rs2s1 = rs1 - rs2
    rhos2s1 = np.linalg.norm(rs2s1[0:2])
    phis2s1 = np.arctan2(rs2s1[1], rs2s1[0])
    is1 = layer_system.layer_number(rs1[2])
    ziss1 = rs1[2] - layer_system.reference_z(is1)
    is2 = layer_system.layer_number(rs2[2])
    ziss2 = rs2[2] - layer_system.reference_z(is2)

    # assert that both particles are in top layer (this function applies only to this special case)
    assert (is1 == (layer_system.number_of_layers() - 1) and is2 == (layer_system.number_of_layers() - 1))

    # wave numbers
    k = omega * layer_system.refractive_indices[is1]
    kz = smuthi.fields.k_z(k_parallel=k_parallel, k=k)

    L = []
    for pol in range(2):
        # layer response
        # (we pick the [0, 1]-component corresponding to downgoing excitation and upgoing response
        # - should be the only non-zero)
        L.append(lay.layersystem_response_matrix(pol, layer_system.thicknesses,
                                                 layer_system.refractive_indices,
                                                 k_parallel, omega, is2, is1)[0, 1])

    # transformation coefficients, indices: pol, n
    Bdag = np.zeros((2, blocksize1), dtype=complex)
    B = np.zeros((2, blocksize2), dtype=complex)

    # multipole order
    m_vec = [np.zeros(blocksize1, dtype=int), np.zeros(blocksize2, dtype=int)]

    for tau in range(2):
        for m in range(-receiving_particle.m_max, receiving_particle.m_max + 1):
            for l in range(max(1, abs(m)), receiving_particle.l_max + 1):
                n = smuthi.fields.multi_to_single_index(tau, l, m, receiving_particle.l_max, receiving_particle.m_max)
                m_vec[0][n] = m
                for pol in range(2):
                    Bdag[pol, n] = trf.transformation_coefficients_vwf(tau, l, m, pol, kp=k_parallel, kz=kz, dagger=True)

    for tau in range(2):
        for m in range(-emitting_particle.m_max, emitting_particle.m_max + 1):
            for l in range(max(1, abs(m)), emitting_particle.l_max + 1):
                n = smuthi.fields.multi_to_single_index(tau, l, m, emitting_particle.l_max, emitting_particle.m_max)
                m_vec[1][n] = m
                for pol in range(2):
                    B[pol, n] = trf.transformation_coefficients_vwf(tau, l, m, pol, kp=k_parallel, kz=-kz, dagger=False)

    # factor
    m2_minus_m1 = m_vec[1] - m_vec[0][np.newaxis].T
    wr_const = 4 * (1j) ** abs(m2_minus_m1) * np.exp(1j * m2_minus_m1 * phis2s1) / (kz * k)

    # result
    g = np.zeros((blocksize1, blocksize2), dtype=complex)
    for pol in range(2):
        g += L[pol] * Bdag[pol, :][np.newaxis].T * B[pol, :]
    g *= wr_const

    # bessel function order
    bessel_order = abs(m2_minus_m1)
    delta_z = ziss1 + ziss2

    return g, rhos2s1, delta_z, bessel_order


def layer_mediated_coupling_block_stat_phase_approx(vacuum_wavelength, receiving_particle, emitting_particle,
                                                    layer_system):
    """
    Compute the layer mediated coupling coefficients by means of the stationary phase approximation, as presented in
    "A quick way to approximate a Sommerfeld-Weyl_type Sommerfeld integral" by W.C. Chew (1988).

    The stationary phase approximation is expected to yield good results for particles with a large lateral distance.

    ********************************************************************************************************************
    Note: This function assumes that both particles (emitter and receiver) are located in the top layer of the layered
    medium. For other cases, this function does not apply.
    ********************************************************************************************************************

    Args:
        vacuum_wavelength (float):                          Vacuum wavelength :math:`\lambda` (length unit)
        receiving_particle (smuthi.particles.Particle):     Particle that receives the scattered field
        emitting_particle (smuthi.particles.Particle):      Particle that emits the scattered field
        layer_system (smuthi.layers.LayerSystem):           Stratified medium in which the coupling takes place

    Returns:
        Stationary phase approximation for the layer mediated coupling matrix block as numpy array.

    """

    # index specs
    lmax1 = receiving_particle.l_max
    mmax1 = receiving_particle.m_max
    lmax2 = emitting_particle.l_max
    mmax2 = emitting_particle.m_max
    blocksize1 = smuthi.fields.blocksize(lmax1, mmax1)
    blocksize2 = smuthi.fields.blocksize(lmax2, mmax2)

    # cylindrical coordinates of relative position vectors
    rs1 = np.array(receiving_particle.position)
    rs2 = np.array(emitting_particle.position)
    rs2s1 = rs1 - rs2
    rhos2s1 = np.linalg.norm(rs2s1[0:2])
    phis2s1 = np.arctan2(rs2s1[1], rs2s1[0])
    is1 = layer_system.layer_number(rs1[2])
    ziss1 = rs1[2] - layer_system.reference_z(is1)
    is2 = layer_system.layer_number(rs2[2])
    ziss2 = rs2[2] - layer_system.reference_z(is2)

    # assert that both particles are in top layer (this function applies only to this special case)
    assert (is1 == (layer_system.number_of_layers() - 1) and is2 == (layer_system.number_of_layers() - 1))

    # wavenumber
    omega = smuthi.fields.angular_frequency(vacuum_wavelength)
    k = omega * layer_system.refractive_indices[is1]

    # stationary phase point
    r = np.sqrt(rhos2s1**2 + (ziss1 + ziss2)**2)
    sin_theta = rhos2s1 / r
    kpar0 = sin_theta * k
    kz0 = k * np.sqrt(1 - sin_theta**2)

    # g function (see [Chew1988])
    g, _, _, bess_ord = g_function(vacuum_wavelength=vacuum_wavelength,
                                   receiving_particle=receiving_particle,
                                   emitting_particle=emitting_particle,
                                   layer_system=layer_system,
                                   k_parallel=kpar0)

    # stationary phase approximation
    hankel = np.zeros((blocksize1, blocksize2), dtype=complex)

    for i1 in range(blocksize1):
        for i2 in range(blocksize2):
            hankel[i1, i2] = scipy.special.hankel1(bess_ord[i1, i2], kpar0 * rhos2s1)

    wr0 = -1j * g * hankel / scipy.special.hankel1(0, kpar0 * rhos2s1) * kz0 * np.exp(1j * k * r) / r

    return wr0


def layer_mediated_coupling_block(vacuum_wavelength, receiving_particle, emitting_particle, layer_system,
                                  k_parallel='default', show_integrand=False):
    """Layer-system mediated particle coupling matrix :math:`W^R` for two particles. This routine is explicit, but slow.

    Args:
        vacuum_wavelength (float):                          Vacuum wavelength :math:`\lambda` (length unit)
        receiving_particle (smuthi.particles.Particle):     Particle that receives the scattered field
        emitting_particle (smuthi.particles.Particle):      Particle that emits the scattered field
        layer_system (smuthi.layers.LayerSystem):           Stratified medium in which the coupling takes place
        k_parallel (numpy ndarray):                         In-plane wavenumbers for Sommerfeld integral
                                                            If 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
        show_integrand (bool):                              If True, the norm of the integrand is plotted.

    Returns:
        Layer mediated coupling matrix block as numpy array.
    """
    if type(k_parallel) == str and k_parallel == 'default':
        k_parallel = smuthi.fields.default_Sommerfeld_k_parallel_array
       
    omega = smuthi.fields.angular_frequency(vacuum_wavelength)

    # index specs
    lmax1 = receiving_particle.l_max
    mmax1 = receiving_particle.m_max
    lmax2 = emitting_particle.l_max
    mmax2 = emitting_particle.m_max
    blocksize1 = smuthi.fields.blocksize(lmax1, mmax1)
    blocksize2 = smuthi.fields.blocksize(lmax2, mmax2)

    # cylindrical coordinates of relative position vectors
    rs1 = np.array(receiving_particle.position)
    rs2 = np.array(emitting_particle.position)
    rs2s1 = rs1 - rs2
    rhos2s1 = np.linalg.norm(rs2s1[0:2])
    phis2s1 = np.arctan2(rs2s1[1], rs2s1[0])
    is1 = layer_system.layer_number(rs1[2])
    ziss1 = rs1[2] - layer_system.reference_z(is1)
    is2 = layer_system.layer_number(rs2[2])
    ziss2 = rs2[2] - layer_system.reference_z(is2)

    # wave numbers
    kis1 = omega * layer_system.refractive_indices[is1]
    kis2 = omega * layer_system.refractive_indices[is2]
    kzis1 = smuthi.fields.k_z(k_parallel=k_parallel, k=kis1)
    kzis2 = smuthi.fields.k_z(k_parallel=k_parallel, k=kis2)

    # phase factors
    ejkz = np.zeros((2, 2, len(k_parallel)), dtype=complex)  # indices are: particle, plus/minus, kpar_idx
    ejkz[0, 0, :] = np.exp(1j * kzis1 * ziss1)
    ejkz[0, 1, :] = np.exp(- 1j * kzis1 * ziss1)
    ejkz[1, 0, :] = np.exp(1j * kzis2 * ziss2)
    ejkz[1, 1, :] = np.exp(- 1j * kzis2 * ziss2)

    # layer response
    L = np.zeros((2, 2, 2, len(k_parallel)), dtype=complex)  # polarization, pl/mn1, pl/mn2, kpar_idx
    for pol in range(2):
        L[pol, :, :, :] = lay.layersystem_response_matrix(pol, layer_system.thicknesses,
                                                          layer_system.refractive_indices, k_parallel, omega, is2, is1)

    # transformation coefficients
    B = [np.zeros((2, 2, blocksize1, len(k_parallel)), dtype=complex),
         np.zeros((2, 2, blocksize2, len(k_parallel)), dtype=complex)]
    # list index: particle, np indices: pol, plus/minus, n, kpar_idx

    m_vec = [np.zeros(blocksize1, dtype=int), np.zeros(blocksize2, dtype=int)]
   
    # precompute spherical functions
    ct = kzis1 / kis1
    st = k_parallel / kis1
    _, pilm_list_pl, taulm_list_pl = sf.legendre_normalized(ct, st, lmax1)
    _, pilm_list_mn, taulm_list_mn = sf.legendre_normalized(-ct, st, lmax1)
    pilm = (pilm_list_pl, pilm_list_mn)
    taulm = (taulm_list_pl, taulm_list_mn)
   
    for tau in range(2):
        for m in range(-mmax1, mmax1 + 1):
            for l in range(max(1, abs(m)), lmax1 + 1):
                n = smuthi.fields.multi_to_single_index(tau, l, m, lmax1, mmax1)
                m_vec[0][n] = m
                for iplmn in range(2):
                    for pol in range(2):
                        B[0][pol, iplmn, n, :] = trf.transformation_coefficients_vwf(tau, l, m, pol, pilm_list=pilm[iplmn],
                                                                                     taulm_list=taulm[iplmn], dagger=True)
   
    ct = kzis2 / kis2
    st = k_parallel / kis2
    _, pilm_list_pl, taulm_list_pl = sf.legendre_normalized(ct, st, lmax2)
    _, pilm_list_mn, taulm_list_mn = sf.legendre_normalized(-ct, st, lmax2)
    pilm = (pilm_list_pl, pilm_list_mn)
    taulm = (taulm_list_pl, taulm_list_mn)
                       
    for tau in range(2):
        for m in range(-mmax2, mmax2 + 1):
            for l in range(max(1, abs(m)), lmax2 + 1):
                n = smuthi.fields.multi_to_single_index(tau, l, m, lmax2, mmax2)
                m_vec[1][n] = m
                for iplmn in range(2):
                    for pol in range(2):
                        B[1][pol, iplmn, n, :] = trf.transformation_coefficients_vwf(tau, l, m, pol, pilm_list=pilm[iplmn],
                                                                                     taulm_list=taulm[iplmn], dagger=False)

    # bessel function and jacobi factor
    bessel_list = []
    for dm in range(lmax1 + lmax2 + 1):
        bessel_list.append(scipy.special.jv(dm, k_parallel * rhos2s1))
    jacobi_vector = k_parallel / (kzis2 * kis2)
    m2_minus_m1 = m_vec[1] - m_vec[0][np.newaxis].T
    wr_const = 4 * (1j) ** abs(m2_minus_m1) * np.exp(1j * m2_minus_m1 * phis2s1) 

    integral = np.zeros((blocksize1, blocksize2), dtype=complex) 
    for n1 in range(blocksize1):
        BeL = np.zeros((2, 2, len(k_parallel)), dtype=complex)  # indices are: pol, plmn2, n1, kpar_idx
        for iplmn1 in range(2):
            for pol in range(2):
                BeL[pol, :, :] += (L[pol, iplmn1, :, :]
                                        * B[0][pol, iplmn1, n1, :]
                                        * ejkz[0, iplmn1, :])
        for n2 in range(blocksize2):
            bessel_full = bessel_list[abs(m_vec[0][n1] - m_vec[1][n2])]
            BeLBe = np.zeros((len(k_parallel)), dtype=complex)
            eval_BeLBe(BeLBe, BeL, B[1], ejkz, n2)
            integrand = bessel_full * jacobi_vector * BeLBe
            integral[n1,n2] = numba_trapz(integrand, k_parallel)
    wr = wr_const * integral

    return wr


def layer_mediated_coupling_matrix(vacuum_wavelength, particle_list, layer_system, k_parallel='default'):
    """Layer system mediated particle coupling matrix W^R for a particle collection in a layered medium.

    Args:
        vacuum_wavelength (float):                                  Wavelength in length unit
        particle_list (list of smuthi.particles.Particle obejcts:   Scattering particles
        layer_system (smuthi.layers.LayerSystem):                   The stratified medium
        k_parallel (numpy.ndarray or str):                          In-plane wavenumber for Sommerfeld integrals.
                                                                    If 'default', smuthi.fields.default_Sommerfeld_k_parallel_array
   
    Returns:
        Ensemble coupling matrix as numpy array.
    """
   
    # indices
    blocksizes = [smuthi.fields.blocksize(particle.l_max, particle.m_max) for particle in particle_list]

    # initialize result
    wr = np.zeros((sum(blocksizes), sum(blocksizes)), dtype=complex)

    for s1, particle1 in enumerate(particle_list):
        idx1 = np.array(range(sum(blocksizes[:s1]), sum(blocksizes[:s1]) + blocksizes[s1]))
        for s2, particle2 in enumerate(particle_list):
            idx2 = range(sum(blocksizes[:s2]), sum(blocksizes[:s2]) + blocksizes[s2])
            wr[idx1[:, None], idx2] = layer_mediated_coupling_block(vacuum_wavelength, particle1, particle2,
                                                                    layer_system, k_parallel)

    return wr

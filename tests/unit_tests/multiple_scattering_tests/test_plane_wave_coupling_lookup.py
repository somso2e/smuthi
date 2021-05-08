''' test pwe coupling lookup table '''

import smuthi.fields as flds
import smuthi.linearsystem.particlecoupling.prepare_lookup as look
import smuthi.linearsystem.particlecoupling.direct_coupling as pacou
import smuthi.utility.cuda as cu
import smuthi.fields
import smuthi.particles as part
import smuthi.layers as lay
import numpy as np

wl = 550
n_lay =  1.5
rho_max = 400
lmax = 3
k_is = 2 * np.pi * n_lay / wl
resolution = 10
 
neffmax = 3
neffimag = 0.01
neff_waypoints = [0, 0.8, 0.8 - 1j * neffimag, 2.1 - 1j * neffimag, 2.1, neffmax]
neff_resolution = 1e-3 
k_parallel = flds.create_k_parallel_array(wl, neff_waypoints, neff_resolution)



''' compare CPU and GPU lookup '''
cu.enable_gpu(enable=False)
w_CPU, rho_array = look.radial_coupling_lookup_table_pwe_correction(vacuum_wavelength=wl,
                                                                rho_max=rho_max,
                                                                l_max=lmax,
                                                                k_is=k_is,
                                                                k_parallel=k_parallel,
                                                                resolution=resolution)


cu.enable_gpu(enable=True)
w_CUDA, rho_array_CUDA = look.radial_coupling_lookup_table_pwe_correction(vacuum_wavelength=wl,
                                                                rho_max=rho_max,
                                                                l_max=lmax,
                                                                k_is=k_is,
                                                                k_parallel=k_parallel,
                                                                resolution=resolution)




''' compare direct evaluation with lookup '''
lay_sys = lay.LayerSystem([0, 0], [n_lay, n_lay])

cylinder1 = part.FiniteCylinder(position=[-100, 0, 400], polar_angle=0, azimuthal_angle=0, refractive_index=2 + 0j,
                 cylinder_radius=100, cylinder_height=200, l_max=lmax, m_max=lmax)
cylinder2 = part.FiniteCylinder(position=[0, 111.80339887498948, 400], polar_angle=0, azimuthal_angle=0, refractive_index=2 + 0j,
                 cylinder_radius=100, cylinder_height=200, l_max=lmax, m_max=lmax)
particle_list = [cylinder1, cylinder2]


emitter = 0
receiver = 1

p1 = np.array(particle_list[receiver].position)  
p2 = np.array(particle_list[emitter].position)  
p1p2 = p2 - p1
azimuth = np.arctan2(p1p2[1], p1p2[0])
elevation = np.arctan2(p1p2[2], (p1p2[0] ** 2 + p1p2[1] ** 2) ** 0.5)

beta = (-np.pi / 2) + elevation
alpha = -azimuth

w_pvwf = pacou.direct_coupling_block_pvwf_mediated(vacuum_wavelength=wl,
                                                   receiving_particle=particle_list[receiver],
                                                   emitting_particle=particle_list[emitter],
                                                   layer_system=lay_sys, 
                                                   k_parallel=k_parallel,
                                                   alpha=alpha,
                                                   beta=beta)





def test_pwe_coupling_radial_lookup():
    relerr = np.linalg.norm(w_CUDA - w_CPU) / np.linalg.norm(w_CPU)
    print('relative difference between CPU and CUDA implementation: ', relerr)
    assert relerr < 5e-4
    
def test_lookup_vs_explicit_evaluation():
    distance = np.linalg.norm(p1p2)
    idx = np.argwhere(rho_array == distance)[0][0]
    
    blocksize = smuthi.fields.blocksize(lmax, lmax)
    # copied from smuthi.linearsystem.linear_system.CouplingMatrixRadialLookupCPU()
    x_array = np.array([particle.position[0] for particle in particle_list])
    y_array = np.array([particle.position[1] for particle in particle_list])
    particle_phi_array = np.arctan2(y_array[:, None] - y_array[None, :], x_array[:, None] - x_array[None, :])
    
    m_list = [None for i in range(blocksize)]
    for i, particle in enumerate(particle_list):
        for m in range(-particle.m_max, particle.m_max + 1):
            for l in range(max(1, abs(m)), particle.l_max + 1):
                for tau in range(2):
                    n_lookup = flds.multi_to_single_index(tau=tau, l=l, m=m, l_max=lmax, m_max=lmax)
                    m_list[n_lookup] = m
    
    phi = particle_phi_array[receiver, emitter]
    w_pvwf_lookup = np.zeros((blocksize, blocksize), dtype=complex)
    for n1 in range(blocksize):
        m1 = m_list[n1]
        for n2 in range(blocksize):
            m2 = m_list[n2]
            M = w_CPU[idx, n1, n2]
            w_pvwf_lookup[n1, n2] = M * np.exp(1j * (m2 - m1) * phi)
            
    relerr = np.linalg.norm(w_pvwf_lookup - w_pvwf) / np.linalg.norm(w_pvwf)
    print('relative difference between direct evaluation and lookup: ', relerr)
    assert relerr < 5e-4
    
    
if __name__ == '__main__':
    test_pwe_coupling_radial_lookup()
    test_lookup_vs_explicit_evaluation()
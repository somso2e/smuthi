""" Test simulation run and far field intensities of periodic particle arrangements """

import smuthi.simulation as sim
import smuthi.initial_field as init
import smuthi.layers as lay
import smuthi.particles as part
import smuthi.fields as flds
import smuthi.periodicboundaries as pb
import smuthi.periodicboundaries.post_processing as pbpost
import numpy as np
from tqdm import tqdm
import sys
import os
import matplotlib.pyplot as plt


# Comsol reference data
filename = 'test_15spheres_laysys_polarsweep_comsol_values.txt'  
with open(filename, 'r') as n_file:
    lines = n_file.readlines()
theta_FEM = np.zeros(len(lines[5:]), dtype=float)
T_FEM = np.zeros(len(lines[5:]), dtype=float)
for idx, line in enumerate(lines[5:]):
    split_line = line.split()
    theta_FEM[idx] = float(split_line[0])
    T_FEM[idx] = np.round(float(split_line[1]), 4)



# set up parameter sweep
wl = 500
polar_angles = np.arange(0, 90, 1) / 180 * np.pi
T = np.zeros(len(polar_angles))
R = np.zeros(len(polar_angles))
A = np.zeros(len(polar_angles))
alpha = 3 / 2 * np.pi

# layer system
layer_system = lay.LayerSystem([0, 3000, 0], [1, 1.3, 1.5])

# particle list
lmax = 3
mmax = 3
rfidx = 2
particle_list = [part.Sphere(position=[-200, -200, 1000], refractive_index=rfidx, radius=150, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[200, 200, 2000], refractive_index=rfidx, radius=100, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[0, 30, 1500], refractive_index=rfidx, radius=120, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[-300, -200, 2700], refractive_index=rfidx, radius=110, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[110, 150, 2720], refractive_index=rfidx, radius=140, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[200, 0, 390], refractive_index=rfidx, radius=120, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[-100, 380, 400], refractive_index=rfidx, radius=100, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[200, -100, 1900], refractive_index=rfidx, radius=100, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[-320, 300, 2050], refractive_index=rfidx, radius=110, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[-360, -360, 220], refractive_index=rfidx, radius=100, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[100, 300, 1160], refractive_index=rfidx, radius=100, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[-150, 300, 2360], refractive_index=rfidx, radius=120, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[-270, -300, 2200], refractive_index=rfidx, radius=130, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[-170, -180, 1470], refractive_index=rfidx, radius=100, l_max=lmax, m_max=mmax),
                  part.Sphere(position=[340, -300, 1500], refractive_index=rfidx, radius=100, l_max=lmax, m_max=mmax)]


neffmax = 3
neffimag = 0.01
waypoints = [0, 0.8, 0.8 - 1j * neffimag, 2.1 - 1j * neffimag, 2.1, neffmax]
neff_discr = 1e-3  
flds.default_Sommerfeld_k_parallel_array = flds.reasonable_Sommerfeld_kpar_contour(vacuum_wavelength=wl,
                                                                                   neff_waypoints=waypoints,
                                                                                   neff_resolution=neff_discr)
flds.angular_arrays(angular_resolution=np.pi/360)


def test_polar_angle_sweep():

    for idx, theta in enumerate(tqdm(polar_angles, desc='Wavelength iterations  ', file=sys.stdout,
                                bar_format='{l_bar}{bar}| elapsed: {elapsed} ' 'remaining: {remaining}')):
        # suppress consol output
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        old_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')
    
        # initial field
        initial_field = init.PlaneWave(vacuum_wavelength=wl,
                                       polar_angle=theta,
                                       azimuthal_angle=alpha,
                                       polarization=0,
                                       amplitude=1,
                                       reference_point=[0, 0, 0])
        
        # define unit cell
        a1 = np.array([1001, 0, 0], dtype=float)
        a2 = np.array([0, 1001, 0], dtype=float) 
        pb.default_Ewald_sum_separation = pb.set_ewald_sum_separation(a1, a2)
        
        # run simulation
        simulation = sim.Simulation(initial_field=initial_field,
                                    layer_system=layer_system,
                                    particle_list=particle_list,
                                    periodicity=(a1, a2),
                                    ewald_sum_separation_parameter='default',
                                    number_of_threads=-2)
        simulation.run()
        
        # total field plane wave expansion in the outer layers
        pwe_total_T = pbpost.transmitted_plane_wave_expansion(initial_field,
                                                              particle_list,
                                                              layer_system,
                                                              a1, a2)
        pwe_total_R = pbpost.reflected_plane_wave_expansion(initial_field,
                                                            particle_list,
                                                            layer_system,
                                                            a1, a2)
        
        
        # farfield objects       
        ff_T = pbpost.periodic_pwe_to_ff_conversion(pwe_total_T,
                                                    simulation.initial_field,
                                                    simulation.layer_system)
        ff_R = pbpost.periodic_pwe_to_ff_conversion(pwe_total_R,
                                                    simulation.initial_field,
                                                    simulation.layer_system)
            
        # power flow per area
        initial_power = pbpost.initial_plane_wave_power_per_area(initial_field, layer_system)
        transmitted_power = pbpost.scattered_periodic_ff_power_per_area(ff_T)
        reflected_power = pbpost.scattered_periodic_ff_power_per_area(ff_R)
        # T, R, A
        T[idx] = transmitted_power / initial_power
        R[idx] = reflected_power / initial_power
        A[idx] = 1 - T[idx] - R[idx]
        
        
        # restore consol output settings
        sys.stdout = old_stdout  
        sys.stderr = old_stderr
    
        fig, ax = plt.subplots()
        ax.plot(polar_angles[:idx] / np.pi * 180, T[:idx], 
                color='C0', label='Transmittance')
        ax.plot(polar_angles[:idx] / np.pi * 180, R[:idx],
                color='C1', label='Reflectance')
        ax.plot(polar_angles[:idx] / np.pi * 180, A[:idx],
                color='C2', label='Absorbance')      
        ax.plot(theta_FEM / np.pi * 180, T_FEM,
                color='C4', linestyle='--', label='Transmittance, Comsol')
        ax.set_xlim(-2, 92)
        ax.set_ylim(-0.05, 1.05)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        ax.set_ylabel('intensity')
        ax.set_xlabel('polar angle (°)')
        ax.legend()  
        plt.show()

if __name__ == '__main__':
    test_polar_angle_sweep()
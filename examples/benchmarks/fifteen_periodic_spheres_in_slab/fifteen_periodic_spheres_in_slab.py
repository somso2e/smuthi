#*****************************************************************************#
# This benchmark compares the transmittance of a plane wave impinging onto a  #
# periodic particle arrangement under various incident angles as computed     #
# with Smuthi to that computed with COMSOL (FEM).                             #
# The tetragonal Bravais lattice hosts fifteen spheres of various sizes per   #
# unit cell that are embedded in a dielectic thin film between a glass        #
# substrate and air.                                                          #
# The script runs with Smuthi version 1.2                                     #
#*****************************************************************************#

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


# set up parameter sweep
wl = 500
polar_angles = np.arange(0, 90, 1) / 180 * np.pi
transmittance = np.zeros(len(polar_angles))
reflectance = np.zeros(len(polar_angles))
absorbance = np.zeros(len(polar_angles))
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



for idx, beta in enumerate(tqdm(polar_angles, desc='Polar angle iterations  ', file=sys.stdout,
                            bar_format='{l_bar}{bar}| elapsed: {elapsed} ' 'remaining: {remaining}')):
    # suppress consol output
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    old_stderr = sys.stderr
    sys.stderr = open(os.devnull, 'w')

    # initial field
    initial_field = init.PlaneWave(vacuum_wavelength=wl,
                                   polar_angle=beta,
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
                                number_of_threads_periodic=-2) # all but 2 threads
    simulation.run()
    
    # plane wave expansion of total transmitted field
    pwe_total_T = pbpost.transmitted_plane_wave_expansion(initial_field,
                                                          particle_list,
                                                          layer_system,
                                                          a1, a2)
    # plane wave expansion of total reflected field
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
    transmittance[idx] = transmitted_power / initial_power
    reflectance[idx] = reflected_power / initial_power
    absorbance[idx] = 1 - transmittance[idx] - reflectance[idx]
    
    
    # restore consol output settings
    sys.stdout = old_stdout  
    sys.stderr = old_stderr


# load COMSOL results 
filename = 'comsol_data_fifteen_periodic_spheres_in_slab.txt'
comsol_data = np.loadtxt(filename, comments='%')
beta_FEM = comsol_data[:, 0]
transmittance_FEM = comsol_data[:, 1]

    
fig, ax = plt.subplots()
ax.plot(polar_angles[:idx] / np.pi * 180, transmittance[:idx],
        linewidth=2, color='C0', label='SMUTHI')
ax.plot(beta_FEM / np.pi * 180, transmittance_FEM,
        linestyle='--', linewidth=2, color='C1', label='COMSOL')

ax.set_xlim(-2, 92)
ax.set_ylim(-0.05, 1.05)
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax.set_ylabel('transmittance')
ax.set_xlabel('polar angle (°)')
ax.legend(loc='lower left')
plt.show()


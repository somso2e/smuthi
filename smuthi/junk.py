#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:55:59 2022

@author: parkerwray
"""


import matplotlib.pyplot as plt
import numpy as np
from smuthi.particles import Sphere
from smuthi.simulation import Simulation
from smuthi.initial_field import PlaneWave
from smuthi.layers import LayerSystem
from smuthi.postprocessing.graphical_output import plot_particles, show_scattered_far_field, show_total_far_field
from smuthi.postprocessing.far_field import scattering_cross_section

radius = 390
period = 1000
max_range = 5
spheres= []
vals = np.arange(-max_range, max_range+1, 1)
for x in vals :
    for y in vals :
        spheres.append(Sphere(position=[x*period, y*period, radius], refractive_index=3.5 + 0j, radius=radius, l_max=4, m_max=4))

plt.figure()
plot_particles(-max_range*period*25, max_range*period*25, -max_range*period, max_range*period, radius, radius, spheres,
                   draw_circumscribing_sphere = True, fill_particle=True)
plt.xlim([-max_range*period-radius,max_range*period+radius])
plt.ylim([-max_range*period-radius,max_range*period+radius])
plt.show()

spheres = [Sphere(position=[x*period, y*period, radius], refractive_index=1 + 0j, radius=radius, l_max=4, m_max=4)]
layers = []
layers = LayerSystem(thicknesses=[0, 0],refractive_indices=[1.45, 1])


plane_wave = []
plane_wave = PlaneWave(vacuum_wavelength = 1340, 
                            polar_angle = 0,
                            azimuthal_angle = 0,
                            polarization = 0, 
                            amplitude=1,
                            reference_point=None)


simulation = []
simulation = Simulation(layer_system=layers,
                        particle_list=spheres,
                        initial_field=plane_wave,
                        solver_type='LU',
                        store_coupling_matrix=True,
                        log_to_terminal=True,
                        check_circumscribing_spheres = True,
                        identical_particles = True)

simulation.run()

dscs = scattering_cross_section(plane_wave, spheres, layers,)

#%%
show_scattered_far_field(simulation)
show_total_far_field(simulation)
test = dscs.azimuthal_integral()
test = np.sum(test, axis = 0)


#%%
plt.figure()
plt.plot(test)


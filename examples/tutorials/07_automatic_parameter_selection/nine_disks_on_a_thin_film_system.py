#*****************************************************************************#
# This is an example script for Smuthi v1.0.0                                 #
# It evaluates the scattered far field from nine metallic cylindric disks on  #
# a substrate covered by a thin dielectric film                               #
# The purpose of this script is to illustrate the process of automatic        #
# parameter selection with the help of a test balloon simulation              #
#*****************************************************************************#

import numpy as np
import smuthi.simulation
import smuthi.initial_field
import smuthi.layers
import smuthi.particles
import smuthi.postprocessing.far_field
import smuthi.postprocessing.graphical_output
import smuthi.utility.automatic_parameter_selection

# optical system
vacuum_wavelength = 550
cylinder_radius = 150
cylinder_height = 100
cylinder_refractive_index = 6.6+1.01j
cylinder_distance = 400
thin_film_thickness = 150
thin_film_refractive_index = 2.1 + 0.03j

# layer system
three_layers = smuthi.layers.LayerSystem(thicknesses=[0, thin_film_thickness, 0],
                                         refractive_indices=[1.52, thin_film_refractive_index, 1])

# initial field
plane_wave = smuthi.initial_field.PlaneWave(vacuum_wavelength=vacuum_wavelength,
                                            polar_angle=np.pi*0.8,
                                            azimuthal_angle=0,
                                            polarization=0)

# particles
# particle list for test balloon simulation
one_disk = [smuthi.particles.FiniteCylinder(position=[0, 0, thin_film_thickness + 0.5*cylinder_height],
                                            refractive_index=cylinder_refractive_index,
                                            cylinder_radius=cylinder_radius,
                                            cylinder_height=cylinder_height,
                                            l_max=3)]   # this setting will be overwritten
                                                        # by the automatic parameter selection

# particle list for actual simulation
nine_disks = []
for i in range(3):
    for j in range(3):
        nine_disks.append(smuthi.particles.FiniteCylinder(position=[(i-1)*cylinder_distance, 
                                                                    (j-1)*cylinder_distance, 
                                                                    thin_film_thickness + 0.5*cylinder_height],
                                                          refractive_index=cylinder_refractive_index,
                                                          cylinder_radius=cylinder_radius,
                                                          cylinder_height=cylinder_height,
                                                          l_max=3))   # this setting will be overwritten later

# test balloon simulation
test_balloon_simulation = smuthi.simulation.Simulation(layer_system=three_layers,
                                                       particle_list=one_disk,
                                                       initial_field=plane_wave)

# run automatic parameter selection
smuthi.utility.automatic_parameter_selection.select_numerical_parameters(test_balloon_simulation,
                                                                         tolerance=1e-2)

# show the resulting parameters
print("Automatic parameter selection finished.")
print("Results for test balloon simulation:")
print("l_max=", one_disk[0].l_max)
print("m_max=", one_disk[0].m_max)
print("neff_max=", test_balloon_simulation.neff_max.real)
print("neff_resolution=", test_balloon_simulation.neff_resolution)
input("Continue with actual simulation? Press any key ...")

# transfer parameters to actual simulation
for disk in nine_disks:
    disk.l_max = one_disk[-1].l_max
    disk.m_max = one_disk[-1].m_max

# actual simulation
simulation = smuthi.simulation.Simulation(layer_system=three_layers,
                                          particle_list=nine_disks,
                                          initial_field=plane_wave,
                                          neff_max=test_balloon_simulation.neff_max,
                                          neff_resolution=test_balloon_simulation.neff_resolution)
simulation.run()

# show far field
smuthi.postprocessing.graphical_output.show_scattering_cross_section(simulation,
                                                                     show_plots=True,
                                                                     show_opts=[{'label':'scattering_cross_section',
                                                                                 'vmin':1e2,    # play with these parameters
                                                                                 'vmax':1e5}],  # to get an appealing result
                                                                     log_scale=True)

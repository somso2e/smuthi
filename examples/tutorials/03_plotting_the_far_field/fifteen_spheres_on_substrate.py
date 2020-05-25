#*****************************************************************************#
# This is a simple example script for Smuthi v1.0.0                           #
# It evaluates the scattered far field from fifteen spheres on a substrate    #
# excited by a plane wave under oblique incidence.                            #
#*****************************************************************************#

import numpy as np
import smuthi.simulation
import smuthi.initial_field
import smuthi.layers
import smuthi.particles
import smuthi.postprocessing.graphical_output as go


# In this file, all lengths are given in nanometers

# Initialize the layer system object containing the substrate (glass) half 
# space and the ambient (air) half space. The coordinate system is such that 
# the interface between the first two layers defines the plane z=0.
# Note that semi infinite layers have thickness 0!
two_layers = smuthi.layers.LayerSystem(thicknesses=[0, 0],
                                       refractive_indices=[1.52, 1])

# particles
def vogel_spiral(number_of_spheres):
    spheres_list = []
    for i in range(1, number_of_spheres):
        r = 200 * np.sqrt(i)
        theta = i * 137.508 * np.pi/180
        spheres_list.append(smuthi.particles.Sphere(position=[r*np.cos(theta),
                                                              r*np.sin(theta),
                                                              100],
                                                    refractive_index=1.52,
                                                    radius=100,
                                                    l_max=3))
    return spheres_list

# Initial field
plane_wave = smuthi.initial_field.PlaneWave(vacuum_wavelength=550,
                                            polar_angle=0.8*np.pi,  # oblique incidence from top
                                            azimuthal_angle=0,
                                            polarization=0)         # 0=TE 1=TM

# Initialize and run simulation
simulation = smuthi.simulation.Simulation(layer_system=two_layers,
                                          particle_list=vogel_spiral(15),
                                          initial_field=plane_wave)
simulation.run()

# plot the scattered far field
go.show_scattering_cross_section(simulation,
                                 log_scale=True,
                                 show_opts=[{'vmin':1e2,   # play with these parameters
                                             'vmax':1e5}]) # to get an appealing result

# note: if the initial field is not a plane wave, we would not use "show_scattering_cross_section", but
# "show_scattered_far_field". E.g.:
# go.show_scattered_far_field(simulation,
#                             log_scale=True,
#                             show_opts=[{'vmin':1e2,   # play with these parameters
#                                         'vmax':1e5}]) # to get an appealing result

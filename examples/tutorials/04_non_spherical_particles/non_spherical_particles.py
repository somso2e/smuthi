#*****************************************************************************#
# This is a simple example script for Smuthi v1.1.0                           #
# It evaluates the total scattering cross section fro three non-spherical     #
# particles located on a glass substrate, excited by a plane wave under       #
# normal incidence.                                                           #
#*****************************************************************************#

import numpy as np
import smuthi.simulation
import smuthi.initial_field
import smuthi.layers
import smuthi.particles
import smuthi.postprocessing.far_field as ff
import smuthi.linearsystem.tmatrix.nfmds.stlmanager

# In this file, all lengths are given in nanometers

# Initialize the layer system object containing the substrate (glass) half 
# space and the ambient (air) half space.
two_layers = smuthi.layers.LayerSystem(thicknesses=[0, 0],
                                       refractive_indices=[1.52, 1])

# Scattering particles
spheroid = smuthi.particles.Spheroid(position=[-500, 0, 100],
                                     refractive_index=2.1+0.01j,
                                     semi_axis_c=100,
                                     semi_axis_a=25,
                                     l_max=5)

cylinder = smuthi.particles.FiniteCylinder(position=[0, 0, 100],
                                           refractive_index=1.52,
                                           cylinder_height=200,
                                           cylinder_radius=100,
                                           l_max=5)


# list of all scattering particles
three_particles = [spheroid, cylinder]

# Initial field
plane_wave = smuthi.initial_field.PlaneWave(vacuum_wavelength=550,
                                            polar_angle=np.pi,    # from top
                                            azimuthal_angle=0,
                                            polarization=0)       # 0=TE 1=TM

# Initialize and run simulation
simulation = smuthi.simulation.Simulation(layer_system=two_layers,
                                          particle_list=three_particles,
                                          initial_field=plane_wave)
simulation.run()

# evaluate the scattering cross section
scs = ff.total_scattering_cross_section(initial_field=plane_wave,
                                        particle_list=three_particles,
                                        layer_system=two_layers)

print("\n****************************************************")
print("Scattering cross section: %e µm^2"%(scs/1e6))
print("****************************************************")

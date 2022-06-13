#*****************************************************************************#
# This is an example script for the use of plane wave coupling in Smuthi      #
# It simulates light scattering by two cylinders on a substrate.              #
# The cylinders are so close that their circumscribing spheres overlap.       #
# For a description of the underlying theory, see
#*****************************************************************************#

import numpy as np
import smuthi.simulation
import smuthi.initial_field
import smuthi.layers
import smuthi.particles
import smuthi.postprocessing.far_field as ff

# In this file, all lengths are given in nanometers

# Initialize the layer system
two_layers = smuthi.layers.LayerSystem(thicknesses=[0, 0],
                                       refractive_indices=[1.52, 1])

# Scattering particles
cylinder1 = smuthi.particles.FiniteCylinder(position=[0, 0, 40],
                                            refractive_index=1+6j,
                                            cylinder_radius=30,
                                            cylinder_height=80,
                                            l_max=7)

# the second cylinder is placed next to the first one, with a gap of only 10nm
cylinder2 = smuthi.particles.FiniteCylinder(position=[70, 0, 40],
                                            refractive_index=1+6j,
                                            cylinder_radius=30,
                                            cylinder_height=80,
                                            l_max=7)

# Initial field
plane_wave = smuthi.initial_field.PlaneWave(vacuum_wavelength=550,
                                            polar_angle=np.pi,
                                            azimuthal_angle=0,
                                            polarization=0)

# Initialize and run simulation
simulation = smuthi.simulation.Simulation(layer_system=two_layers,
                                          particle_list=[cylinder1, cylinder2],
                                          initial_field=plane_wave,
                                          use_pvwf_coupling=True,
                                          pvwf_coupling_neff_max=3,
                                          pvwf_coupling_neff_resolution=1e-3)

simulation.run()

# evaluate the scattering cross section
scs = ff.total_scattering_cross_section(simulation=simulation)

print("\n****************************************************")
print("Scattering cross section: %e µm^2"%(scs/1e6))
print("****************************************************")

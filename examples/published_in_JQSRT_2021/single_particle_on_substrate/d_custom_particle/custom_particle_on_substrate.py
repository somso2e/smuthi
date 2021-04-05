#*****************************************************************************#
# This is a script for Smuthi v1.1.0                                          #
# It computes the total scattering cross section of a particle on a substrate #
#*****************************************************************************#


import numpy as np
import matplotlib.pyplot as plt
from smuthi.simulation import Simulation
from smuthi.initial_field import PlaneWave
from smuthi.layers import LayerSystem
from smuthi.particles import Sphere, Spheroid, FiniteCylinder, CustomParticle
from smuthi.postprocessing.far_field import total_scattering_cross_section
from smuthi.postprocessing.far_field import extinction_cross_section
from smuthi.utility.automatic_parameter_selection import select_numerical_parameters
from smuthi.utility.optical_constants import  read_refractive_index_from_yaml

# In this file, all lengths are given in micrometers

# Layer system (glass substrate under vacuum)
layers = LayerSystem(thicknesses=[0, 0],
                     refractive_indices=[1.52, 1])

# Particle
scale = 1
particle = CustomParticle(position=[0, 0, 0.326],
                          geometry_filename="freeform_particle_vol0.268_minz-0.326.stl",
                          refractive_index=1.52,
                          scale=scale,
                          l_max=1,
                          m_max=1)

# Initial field
plane_wave = PlaneWave(vacuum_wavelength=0.4,
                       polar_angle=np.pi,
                       azimuthal_angle=0,
                       polarization=0)

# Simulation
simulation = Simulation(layer_system=layers,
                        particle_list=[particle],
                        initial_field=plane_wave)

# Find suitable numerical settings
select_numerical_parameters(simulation,
                            tolerance=5e-3,
                            tolerance_factor=0.5,
                            max_iter=30,
                            neff_max_increment=0.2,
                            neff_max_offset=0.5,
                            show_plot=True)

# Launch --------------------------------------------------------
wavelengths = np.arange(0.370, 0.850, 0.01)
scs_list = []
ecs_list = []
for wl in wavelengths:
    print("running wl = ", wl)
    simulation.k_parallel = "default"
    simulation.initial_field.vacuum_wavelength = wl
    simulation.run()
    scs = total_scattering_cross_section(simulation)
    ecs = extinction_cross_section(simulation)
    scs_list.append(scs.real)
    ecs_list.append(ecs.real)

np.savetxt("results.dat", np.transpose([wavelengths, ecs_list, scs_list]), header="wl(mum) ecs(mum2) scs(mum2)")

# Plot
benchmark_data = np.loadtxt("fdtd.dat")
plt.figure(figsize=(4,3))
plt.plot(wavelengths, scs_list)
plt.plot(benchmark_data[:, 0], benchmark_data[:, 1], ".")
plt.xlabel("wavelength (micron)")
plt.ylabel("cross section (square micron)")
plt.legend(["SCS (smuthi)", "SCS (FDTD)"])
plt.tight_layout()
plt.savefig("SCS.png")
plt.show()

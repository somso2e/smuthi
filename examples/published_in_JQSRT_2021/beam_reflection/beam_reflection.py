#*****************************************************************************#
# This is a script for Smuthi v1.1.0                                          #
# It simulates a Gaussian beam impinging on a micro paint layer on an iron    #
# substrate.                                                                  # 
# Warining: depending on the particle number, this script can take a long     #
#           to run                                                            #
#*****************************************************************************#

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from smuthi.simulation import Simulation
from smuthi.initial_field import GaussianBeam
from smuthi.layers import LayerSystem
from smuthi.particles import Sphere
from smuthi.postprocessing.far_field import total_far_field
from smuthi.postprocessing.graphical_output import show_near_field, plot_particles
from smuthi.utility.automatic_parameter_selection import select_numerical_parameters
from smuthi.utility import cuda
import smuthi.linearsystem.linear_system as linsys

cuda.enable_gpu()

# all lengths are given in micrometers
paint_thickness = 1
mean_radius = 0.1
std_radius = 0.01
volume_fraction = 0.2
resin_refractive_index = 1.50
tio2_refractive_index = 2.67
iron_refractive_index = 2.95 + 2.93j
wavelength = 0.55
num_particles = 10000
output_dir = "results%ipart%ivol" % (num_particles, 100 * volume_fraction)

if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

vol = num_particles * 4 / 3 * np.pi * mean_radius**3
rho_max = np.sqrt(vol / paint_thickness / np.pi / volume_fraction)
def random_sequential_addition():
    with open(output_dir + "/particles.dat", 'w') as file:
        particles = []
        volume = 0
        while len(particles) < num_particles:
            rad = np.random.normal(mean_radius, std_radius)
            while True:
                print("adding %i of %i"%(len(particles), num_particles))
                x = np.random.uniform(-rho_max, rho_max)
                y = np.random.uniform(-rho_max, rho_max)
                z = np.random.uniform(rad, paint_thickness - rad)
                if (x ** 2 + y ** 2 > rho_max ** 2):
                    print("outa space")
                    continue
                any_overlap = False
                for part in particles:
                    distance2 = (x - part.position[0])**2 + (y - part.position[1])**2 + (z - part.position[2])**2
                    if distance2 < (rad + part.radius) ** 2:
                        any_overlap = True
                        print("overlap:", distance2, rad, part.radius, rho_max, z)
                        break
                if any_overlap:
                    continue
                sph = Sphere(position=[x,y,z],
                             refractive_index=tio2_refractive_index,
                             radius=rad,
                             l_max=3,
                             m_max=3)
                particles.append(sph)
                volume = volume + 4 / 3 * np.pi * sph.radius ** 3
                file.write("%f,%f,%f,%f\n" % (x, y, z, rad))
                break
        print("Volume fraction:", volume / (paint_thickness * np.pi * rho_max ** 2))
    return particles

layers = LayerSystem(thicknesses=[0, paint_thickness, 0],
                     refractive_indices=[iron_refractive_index, resin_refractive_index, 1])

beam = GaussianBeam(vacuum_wavelength=wavelength,
                    polar_angle=np.pi*3/4,
                    azimuthal_angle=0,
                    polarization=0,
                    beam_waist=2)

simulation = Simulation(layer_system=layers,
                        particle_list=random_sequential_addition(),
                        initial_field=beam,
                        length_unit='micron',
                        solver_type='GCROTMK',
                        solver_tolerance=5e-3,
                        store_coupling_matrix=False,
                        coupling_matrix_lookup_resolution=0.01,
                        coupling_matrix_interpolator_kind="linear",
                        log_to_file=True)

xmax = rho_max + 0.4
zmin = -0.1
zmax = 2

if False:
    plot_particles(xmin=-xmax, xmax=xmax, ymin=0, ymax=0, zmin=-0.1, zmax=1.5,
                   particle_list=simulation.particle_list,
                   draw_circumscribing_sphere=False, fill_particle=True)
    plt.gca().set_aspect('equal')
    plt.xlim([-xmax, xmax])
    plt.ylim([zmin, zmax])
    plt.xlabel("x [micron]")
    plt.ylabel("z [micron]")
    plt.show()

simulation.save("before_run.smuthi")
simulation.run()
simulation.save("after_run.smuthi")

polar_angles = np.arange(0, np.pi/2, np.pi / 360)  # 0.5 degree resolution
total_ff, init_ff, scat_ff = total_far_field(initial_field=simulation.initial_field,
                                             particle_list=simulation.particle_list,
                                             layer_system=simulation.layer_system,
                                             polar_angles=polar_angles)
total_reflected_power = total_ff.integral()[0] + total_ff.integral()[1]
initial_reflected_power = init_ff.integral()[0] + init_ff.integral()[1]
initial_power = beam.initial_intensity(layer_system=layers).integral()[0] + beam.initial_intensity(layer_system=layers).integral()[1]

total_reflection = total_reflected_power / initial_power
initial_reflection = initial_reflected_power / initial_power

print("total reflection:", total_reflection)
print("initial reflection:", initial_reflection)

# warning, near field evaluation takes long time
if False:
    show_near_field(simulation=simulation,
                    show_plots=True,
                    quantities_to_plot=["norm(E)", "E_x", "E_y"],
                    show_opts=[{},
                               {}],
                    save_plots=True,
                    save_opts=[{'format': 'pdf', 'dpi': 200},
                               {'format': 'gif'}],
                    outputdir=output_dir,
                    xmin=-xmax,
                    xmax=xmax,
                    ymin=0,
                    ymax=0,
                    zmin=zmin,
                    zmax=zmax,
                    resolution_step=0.05,
                    show_internal_field=True,
                    save_data=True,
                    data_format='ascii')

if False:
    show_far_field(far_field=total_ff[0])

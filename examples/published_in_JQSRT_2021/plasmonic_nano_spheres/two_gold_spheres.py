#*****************************************************************************#
# This is a script for Smuthi v1.1.0                                          #
# It computes the field enhancement in the gap between two gold nano spheres  #
# on a substrate                                                              # 
#*****************************************************************************#

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from smuthi.simulation import Simulation
from smuthi.initial_field import PlaneWave
from smuthi.layers import LayerSystem
from smuthi.particles import Sphere
from smuthi.postprocessing.graphical_output import show_near_field
from smuthi.postprocessing.scattered_field import scattered_field_piecewise_expansion
from smuthi.postprocessing.far_field import extinction_cross_section
from smuthi.utility.automatic_parameter_selection import select_numerical_parameters
from smuthi.utility.optical_constants import read_refractive_index_from_yaml
from smuthi.utility import cuda


# In this file, all lengths are given in micrometers

# Decide what to compute:
compute_gap_field = True    # calculate the field enhancement in the middle of the gap
compute_near_field = False  # create a field distribution image
compute_lmax = False        # run a convergence analysis with regard to l_max

if compute_near_field:
    cuda.enable_gpu()

radius = 0.03

def set_simulation(wavelength, distance, lmax=12, mmax=12):

    n_au = read_refractive_index_from_yaml("johnson.yml", wavelength)[1]
    n_si = read_refractive_index_from_yaml("aspnes.yml", wavelength)[1]

    print("Wavelength: ", wavelength)
    print("Gold refractive index: ", n_au)
    print("Silicon refractive index: ", n_si)

    particle1 = Sphere(position=[-radius - distance / 2, 0, radius],
                       refractive_index=n_au,
                       radius=radius,
                       l_max=lmax,
                       m_max=mmax)

    particle2 = Sphere(position=[radius + distance / 2, 0, radius],
                       refractive_index=n_au,
                       radius=radius,
                       l_max=lmax,
                       m_max=mmax)

    particle_list = [particle1, particle2]

    layers = LayerSystem(thicknesses=[0, 0],
                         refractive_indices=[n_si, 1.33])

    plane_wave = PlaneWave(vacuum_wavelength=wavelength,
                           polar_angle=np.pi,  # normal incidence from top
                           azimuthal_angle=0,
                           polarization=1)

    return Simulation(layer_system=layers,
                      particle_list=particle_list,
                      initial_field=plane_wave,
                      neff_max=7,
                      length_unit='micron',
                      solver_type='lu')


def gap_field(simulation):
    pfe = scattered_field_piecewise_expansion(vacuum_wavelength=simulation.initial_field.vacuum_wavelength,
                                              particle_list=simulation.particle_list,
                                              layer_system=simulation.layer_system,
                                              layer_numbers=[1])
    pfe = pfe + simulation.initial_field.piecewise_field_expansion(simulation.layer_system)
    ex, ey, ez = pfe.electric_field(0, 0, radius)
    return np.sqrt(abs(ex)**2 + abs(ey)**2 + abs(ez)**2)


def run(wavelengths, distances, autoparam=True):
    simulation = set_simulation(wavelength=wavelengths[0], distance=distances[0])
    if autoparam:
        select_numerical_parameters(simulation,
                                    detector=gap_field,
                                    select_neff_max=True,
                                    tolerance=0.005,
                                    tolerance_factor=0.25,
                                    max_iter=100,
                                    neff_max_offset=1,
                                    show_plot=True)

    l_max = simulation.particle_list[0].l_max
    m_max = simulation.particle_list[0].m_max
    neff_max = simulation.neff_max

    gap_list = []
    for dist in distances:
        print("running d = ", dist)
        simulation = set_simulation(wavelengths[0], dist)
        for part in simulation.particle_list:
            part.l_max = l_max
            part.m_max = m_max
        simulation.neff_max = neff_max
        simulation.run()
        gap_list.append(gap_field(simulation))

    return gap_list


if compute_gap_field:
    distances = [0.02, 0.01, 0.005]
    wavelengths = np.arange(0.4, 0.801, 0.01)
    for distance in distances:
        gap_fields = []
        wl_list = []
        ece_list = []
        for wavelength in wavelengths:
            print("running d = ", distance, "wl = ", wavelength)
            simulation = set_simulation(wavelength=wavelength, distance=distance)
            simulation.run()
            gap_fields.append(gap_field(simulation))
            wl_list.append(wavelength)
            ece_list.append(extinction_cross_section(simulation)["top"] + extinction_cross_section(simulation)["bottom"])
            np.savetxt("results/gap_fields_dist%i.dat" % (distance * 1000), np.transpose([wl_list, gap_fields, ece_list]))

if compute_near_field:
    simulation = set_simulation(wavelength=0.625, distance=0.005)
    simulation.run()
    vmax = 50
    show_near_field(simulation=simulation,
                    show_plots=False,
                    quantities_to_plot=["norm(E)", "E_x"],
                    show_opts=[{'vmin': -vmax, 'vmax': vmax},
                               {'vmin': 0, 'vmax': vmax},
                               {'vmin': -vmax, 'vmax': vmax}],
                    save_plots=True,
                    save_opts=[{'format': 'pdf', 'dpi': 200},
                               {'format': 'png', 'dpi': 200},
                               {'format': 'gif'}],
                    outputdir="results",
                    xmin=-0.1,
                    xmax=0.1,
                    ymin=0,
                    ymax=0,
                    zmin=-0.02,
                    zmax=0.12,
                    resolution_step=0.001,
                    show_internal_field=True,
                    save_data=True,
                    data_format='ascii')


if compute_lmax:
    distances = [0.100, 0.050, 0.020, 0.010, 0.005, 0.002, 0.001]
    wavelength = 0.635
    for dist in distances:
        simulation = set_simulation(wavelength=wavelength, distance=dist, lmax=1, mmax=1)
        simulation.run()
        E = gap_field(simulation)
        ecs = extinction_cross_section(simulation)['top'] + extinction_cross_section(simulation)['bottom']
        lmax_list = [1]
        e_list = [E]
        de_list = [1]
        ecs_list = [ecs]
        decs_list = [1]
        for lmax in range(2, 21):
            simulation = set_simulation(wavelength=wavelength, distance=dist, lmax=lmax, mmax=lmax)
            simulation.run()
            E = gap_field(simulation)
            ecs = extinction_cross_section(simulation)['top'] + extinction_cross_section(simulation)['bottom']
            lmax_list.append(lmax)
            de_list.append(abs(E - e_list[-1]) / abs(E))
            e_list.append(E)
            decs_list.append(abs(ecs - ecs_list[-1]) / abs(ecs))
            ecs_list.append(ecs)
            np.savetxt("results/lmax_dist%i.dat" % (dist * 1000), np.transpose([lmax_list, e_list, de_list, ecs_list, decs_list]))

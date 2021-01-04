import matplotlib.pyplot as plt
import unittest
import numpy as np
import smuthi
import smuthi.simulation
import smuthi.initial_field
import smuthi.layers
import smuthi.particles
import smuthi.postprocessing.far_field as ff
import smuthi.utility.optical_constants as opt
import io
import yaml

        
# for single wavelength
def get_extinctions(wavelength, sphere_refractive_index, 
                                layers_thicknesses, layers_refractive_indeces, polarization, radius):
    spacer = 2 # nm

    layers = smuthi.layers.LayerSystem(thicknesses=layers_thicknesses,
                                       refractive_indices=layers_refractive_indeces)

    # Scattering particle
    sphere = smuthi.particles.Sphere(position=[0, 0, - radius - spacer],
                                 refractive_index=sphere_refractive_index,
                                 radius=radius,
                                 l_max=3)

    # list of all scattering particles (only one in this case)
    spheres_list = [sphere]

    #Initial field
    plane_wave = smuthi.initial_field.PlaneWave(vacuum_wavelength=wavelength,
                                    polar_angle=(np.pi/2-np.pi/7.2), # 25 grad to the surface
                                    azimuthal_angle=0,
                                    polarization=polarization)       # 0=TE 1=TM

    # Initialize and run simulation
    simulation = smuthi.simulation.Simulation(layer_system=layers,
                                            particle_list=spheres_list,
                                            initial_field=plane_wave,
                                            solver_type='LU',
                                            solver_tolerance=1e-5,
                                            store_coupling_matrix=True,
                                            coupling_matrix_lookup_resolution=None,
                                            coupling_matrix_interpolator_kind='cubic')
    simulation.run()

    # Since this moment we count extinction
    # Count total extinction of dipole (another multipoles are cut by only_l=1)
    general_ecs = ff.extinction_cross_section(initial_field=plane_wave,
                                        particle_list=spheres_list,
                                        layer_system=layers, only_l=1, only_pol=polarization)

    extinctions = [[], [], []]
    # extinctions[0] -- magnetic projections,
    # extinctions[1] -- electric projections
    # extinctions[2] -- total extinction
    for tau in [0, 1]:
        for m in [-1, 0, 1]:
            counted_extinction = ff.extinction_cross_section(initial_field=plane_wave,
                                        particle_list=spheres_list,
                                        layer_system=layers, only_l=1, only_pol=polarization, 
                                        only_tau=tau, only_m=m)
            extinctions[tau].append(counted_extinction)
    extinctions[2].append(general_ecs)

    return extinctions


def add_extinctions_to_output(wavelengths, extinctions, subplot, radius, tau,
                                color_for_single_line, color_for_twin_lines,
                                x_description, y_description, z_description):
    x_component = np.zeros(len(extinctions))
    z_component = np.zeros(len(extinctions))
    y_component = np.zeros(len(extinctions))
    total = np.zeros(len(extinctions))

    # Here we convert extinction from spherical to cartesian projection.
    for i in range(len(extinctions)):
        # extinctions[i][j][k], where [i] responses for wavelength,
        # [j] responses for magnetic/electric component (0 -- magnetic, 1 -- electric, 2 -- total),
        # [k] responses for projection ([-1, 0, 1] in our case)
        x_component[i] = (extinctions[i][tau][0] + \
                            extinctions[i][tau][2]).real / \
                            np.pi / radius**2

        y_component[i] = (extinctions[i][1 - tau][0] + \
                            extinctions[i][1 - tau][2]).real / \
                            np.pi / radius**2

        z_component[i] = (extinctions[i][tau][1]).real / \
                            np.pi / radius**2

        total[i] = (extinctions[i][2][0]).real / np.pi / radius**2

    subplot.plot(wavelengths, x_component, linestyle='dashed', color=color_for_twin_lines, \
                            label=x_description)
    subplot.plot(wavelengths, y_component, linestyle='dashed', color=color_for_single_line, \
                            label=y_description)
    subplot.plot(wavelengths, z_component, color=color_for_twin_lines, \
                            label=z_description)
    subplot.plot(wavelengths, total, color='black', label='total')
    subplot.legend()


wavelength_min = 500
wavelength_max = 800
total_points = 151
wavelengths = np.linspace(wavelength_min, wavelength_max, total_points)

index_Si = opt.read_refractive_index_from_yaml('Si-Green-2008.yml', wavelengths, units="nm") 
index_Au = opt.read_refractive_index_from_yaml('Au-Johnson.yml', wavelengths, units="nm")
index_glass = 1.5
index_media = 1.0
radius = 170.0 / 2.0
Au_thickness = 40.0

extinctions_te_glass = []
extinctions_tm_glass = []
extinctions_te_gold = []
extinctions_tm_gold = []

for i in range(len(wavelengths)):
    # Without gold, TE polarization
    decomposed_te_glass_extinction = get_extinctions(wavelengths[i], index_Si[i][1], 
                            [0, 0], [index_media, index_glass], 0, radius)

    # Without gold, TM polarization
    decomposed_tm_glass_extinction = get_extinctions(wavelengths[i], index_Si[i][1], 
                            [0, 0], [index_media, index_glass], 1, radius)

    # With gold, TE polarization
    decomposed_te_gold_extinction = get_extinctions(wavelengths[i], index_Si[i][1], 
                            [0, Au_thickness, 0], [index_media, index_Au[i][1], index_glass], 0, radius)

    # With gold, TM polarization
    decomposed_tm_gold_extinction = get_extinctions(wavelengths[i], index_Si[i][1], 
                            [0, Au_thickness, 0], [index_media, index_Au[i][1], index_glass], 1, radius)

    extinctions_te_glass.append(decomposed_te_glass_extinction)
    extinctions_tm_glass.append(decomposed_tm_glass_extinction)
    extinctions_te_gold.append(decomposed_te_gold_extinction)
    extinctions_tm_gold.append(decomposed_tm_gold_extinction)

# Since this moment we just output the counted extinction decompositions.
f, axarr = plt.subplots(2, 2, sharex=True)
plt.annotate(text="TE", xy=(0.25, 0.9), xycoords='figure fraction')
plt.annotate(text="TM", xy=(0.75, 0.9), xycoords='figure fraction')
plt.annotate(text="glass substrate", xy=(0.05, 0.6), xycoords='figure fraction', rotation=90)
plt.annotate(text="gold substrate", xy=(0.05, 0.2), xycoords='figure fraction', rotation=90)

# Now we use tau=0 for TE-polarization and tau=1 for TM, 
# because it defines which component (magnetic or electric) we will use in plotting x and z axis,
# and which component -- for y axis.
add_extinctions_to_output(wavelengths, extinctions_te_glass, axarr[0, 0], radius, 0, 'red', 'blue', \
                                    '$m_x$', '$p_y$', '$m_z$')
add_extinctions_to_output(wavelengths, extinctions_te_gold, axarr[1, 0], radius, 0, 'red', 'blue', \
                                    '$m_x$', '$p_y$', '$m_z$')
add_extinctions_to_output(wavelengths, extinctions_tm_glass, axarr[0, 1], radius, 1, 'blue', 'red', \
                                    '$p_x$', '$m_y$', '$p_z$')
add_extinctions_to_output(wavelengths, extinctions_tm_gold, axarr[1, 1], radius, 1, 'blue', 'red', \
                                    '$p_x$', '$m_y$', '$p_z$')

axarr[0, 0].set_ylim(ymin=0, ymax=9)
axarr[0, 1].set_ylim(ymin=0, ymax=9)
axarr[1, 0].set_ylim(ymin=-1, ymax=18)
axarr[1, 1].set_ylim(ymin=-1, ymax=18)
plt.show()

# You can compare the output results with the same calculation in paper 
# "Laser Photonics Rev.10,No. 5, 799–806 (2016)/DOI10.1002/lpor.201600055"
# (graphics in page 803)
# -*- coding: utf-8 -*-
"""This script runs a simulation for two spheres in a slab waveguide, illuminated by a plane wave."""

import sys
import numpy as np
import smuthi.particles as part
import smuthi.layers as lay
import smuthi.initial_field as init
import smuthi.simulation as simul
import smuthi.postprocessing.far_field as farf
import smuthi.fields as flds


# Parameter input ----------------------------
vacuum_wavelength = 550
plane_wave_polar_angle = np.pi
plane_wave_azimuthal_angle = 0
plane_wave_polarization = 0
plane_wave_amplitude = 1
lmax = 3
neff_waypoints = [0, 0.5, 0.8-0.01j, 2-0.01j, 2.5, 20]
neff_discr = 1e-3

# --------------------------------------------

# initialize particle object
part1 = part.Sphere(position=[100,100,150], refractive_index=2.4+0.0j, radius=120, l_max=lmax)
part2 = part.Sphere(position=[-100,-100,250], refractive_index=1.9+0.1j, radius=120, l_max=lmax)

# initialize layer system object
lay_sys1 = lay.LayerSystem([0, 400, 0], [1.5, 1.7, 1])
lay_sys2 = lay.LayerSystem([0, 200, 200, 0], [1.5, 1.7, 1.7, 1])

# initialize initial field object
plane_wave = init.PlaneWave(vacuum_wavelength=vacuum_wavelength, polar_angle=plane_wave_polar_angle,
                            azimuthal_angle=plane_wave_azimuthal_angle, polarization=plane_wave_polarization,
                            amplitude=plane_wave_amplitude, reference_point=[0, 0, 400])

# initialize simulation object
simulation1 = simul.Simulation(layer_system=lay_sys1, particle_list=[part1,part2], initial_field=plane_wave,
                               neff_waypoints=neff_waypoints, neff_resolution=neff_discr,
                               log_to_terminal=(not sys.argv[0].endswith('nose2')))  # suppress output if called by nose
simulation1.run()

simulation2 = simul.Simulation(layer_system=lay_sys2, particle_list=[part1,part2], initial_field=plane_wave,
                               neff_waypoints=neff_waypoints, neff_resolution=neff_discr,
                               log_to_terminal=(not sys.argv[0].endswith('nose2')))  # suppress output if called by nose
simulation2.run()

ff = farf.scattered_far_field(vacuum_wavelength=vacuum_wavelength, particle_list=simulation1.particle_list,
                              layer_system=simulation1.layer_system)


def test_equivalent_layer_systems():
    relerr = (np.linalg.norm(simulation1.linear_system.coupling_matrix.linear_operator.A
                             - simulation2.linear_system.coupling_matrix.linear_operator.A)
              / np.linalg.norm(simulation1.linear_system.coupling_matrix.linear_operator.A))
    print(relerr)
    assert relerr < 1e-3
    print('Tests Passed!')

def test_against_prototype():
    b00 = -0.8609045 + 0.4019615j
    assert abs((simulation1.particle_list[0].scattered_field.coefficients[0] - b00) / b00) < 1e-4

    b01 = 0.0223519 + 0.0675757j
    assert abs((simulation1.particle_list[0].scattered_field.coefficients[1] - b01) / b01) < 1e-4

    b027 = -0.0258275 + 0.1145384j
    assert abs((simulation1.particle_list[0].scattered_field.coefficients[27] - b027) / b027) < 1e-4

    b10 = 0.2065701 + 0.1197903j
    assert abs((simulation1.particle_list[1].scattered_field.coefficients[0] - b10) / b10) < 1e-4

    top_power_flux = 4.3895865e+02
    print(abs((sum(ff.top().integral()) - top_power_flux) / top_power_flux))
    assert abs((sum(ff.top().integral()) - top_power_flux) / top_power_flux) < 1e-3

    bottom_power_flux = 2.9024410e+04
    print(abs((sum(ff.bottom().integral()) - bottom_power_flux) / bottom_power_flux))
    assert abs((sum(ff.bottom().integral()) - bottom_power_flux) / bottom_power_flux) < 1e-3
    print('Tests Passed!')

if __name__ == '__main__':
    test_equivalent_layer_systems()
    test_against_prototype()

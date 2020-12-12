import smuthi.initial_field as init
import smuthi.particles as part
import smuthi.simulation as simul
import smuthi.layers as lay
import smuthi.postprocessing.scattered_field as sf
import numpy as np


def test_electric_field_spa():
    ld = 550

    # initialize particles
    sphere1 = part.Sphere(position=[0, 0, 110], refractive_index=2.4 + 0.0j, radius=110, l_max=1, m_max=1)
    sphere2 = part.Sphere(position=[-2000, -2000, 120], refractive_index=2.4 + 0.0j, radius=120, l_max=5, m_max=4)
    sphere3 = part.Sphere(position=[-2000, 2000, 90], refractive_index=2.5 + 0.0j, radius=90, l_max=2, m_max=2)
    part_list = [sphere1, sphere2, sphere3]

    # water over aluminum
    lay_sys = lay.LayerSystem([0, 0], [1+6j, 1])

    # initial field
    plane_wave = init.PlaneWave(vacuum_wavelength=ld, polar_angle=np.pi, azimuthal_angle=0, polarization=0)

    # run simulation
    simulation = simul.Simulation(layer_system=lay_sys, particle_list=part_list, initial_field=plane_wave,
                                  neff_resolution=1e-4)
    simulation.run()

    r = 5e4
    polar = np.pi * 1/5
    azim = np.pi * 0.9
    x = r * np.sin(polar) * np.cos(azim)
    y = r * np.sin(polar) * np.sin(azim)
    z = r * np.cos(polar)

    scat_fld_exp = sf.scattered_field_piecewise_expansion(ld, part_list, lay_sys, layer_numbers=[1])
    ex, ey, ez = scat_fld_exp.electric_field(x, y, z)
    ex_spa, ey_spa, ez_spa = sf.evaluate_scattered_field_stat_phase_approx(x, y, z, ld, part_list, lay_sys)

    absolute_error = np.linalg.norm([ex - ex_spa, ey - ey_spa, ez - ez_spa])
    relative_error = absolute_error / np.linalg.norm([ex, ey, ez])
    print("relative error:", relative_error)
    assert relative_error < 2e-3


if __name__ == '__main__':
    test_electric_field_spa()

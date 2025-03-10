import smuthi.initial_field as init
import smuthi.particles as part
import smuthi.simulation as simul
import smuthi.layers as lay
import smuthi.postprocessing.scattered_field as sf
import smuthi.postprocessing.far_field as ff
import smuthi.fields as flds
import numpy as np


def test_electric_field_ffaprox():
    """
    With this routine we test equation (3.10) of Bohren & Huffman, i.e.,
    the far field approximation E ~ exp(ikr) / (-ikr) A
    for large kr
    """
    ld = 550
    k = flds.angular_frequency(ld)

    # initialize particles
    sphere1 = part.Sphere(position=[0, 0, 110], refractive_index=2.4 + 0.0j, radius=110, l_max=1, m_max=1)
    sphere2 = part.Sphere(position=[-2000, -2000, 120], refractive_index=2.4 + 0.0j, radius=120, l_max=5, m_max=4)
    sphere3 = part.Sphere(position=[-2000, 2000, 90], refractive_index=2.5 + 0.0j, radius=90, l_max=2, m_max=2)
    part_list = [sphere1, sphere2, sphere3]

    # water over aluminum
    lay_sys = lay.LayerSystem([0, 0], [1, 1])

    # initial field
    plane_wave = init.PlaneWave(vacuum_wavelength=ld, polar_angle=np.pi, azimuthal_angle=0, polarization=0)

    # run simulation
    simulation = simul.Simulation(layer_system=lay_sys, particle_list=part_list, initial_field=plane_wave,
                                  neff_resolution=5e-4)
    simulation.run()
    
    scat_ff = ff.scattered_far_field(ld, part_list, lay_sys)
    
    ith = 12
    iph = 18
    polar = scat_ff.polar_angles[ith]
    azim = scat_ff.azimuthal_angles[iph]
    r = 1e8
    x = r * np.sin(polar) * np.cos(azim)
    y = r * np.sin(polar) * np.sin(azim)
    z = r * np.cos(polar)

    scat_fld_exp = sf.scattered_field_piecewise_expansion(ld, part_list, lay_sys, layer_numbers=[1])
    ex, ey, ez = scat_fld_exp.electric_field(x, y, z)
    
    Ax, Ay, Az = scat_ff.electric_field_amplitude()
    prefac = np.exp(1j * k * r) / (-1j * k * r)
    
    ex2 = prefac * Ax[ith, iph]
    ey2 = prefac * Ay[ith, iph]
    ez2 = prefac * Az[ith, iph]

    absolute_error = np.linalg.norm([ex - ex2, ey - ey2, ez - ez2])
    relative_error = absolute_error / np.linalg.norm([ex, ey, ez])
    print("relative error:", relative_error)
    assert relative_error < 5e-3
    print("OK.")


if __name__ == '__main__':
    test_electric_field_ffaprox()

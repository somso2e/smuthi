import smuthi.initial_field as init
import smuthi.particles as part
import smuthi.simulation as simul
import smuthi.layers as lay
import smuthi.postprocessing.scattered_field as sf
import smuthi.postprocessing.far_field as ff
import smuthi.fields as flds
import numpy as np


def test_pwe_translation():
    """
    With this routine we test if setting a new rerefenence point to a plane wave expansion
    yields the expected behaviour.
    """
    ld = 550
    k = flds.angular_frequency(ld)

    # initialize particles
    sphere1 = part.Sphere(position=[0, 0, 110], refractive_index=2.4 + 0.0j, radius=110, l_max=1, m_max=1)
    sphere2 = part.Sphere(position=[-2000, -2000, 120], refractive_index=2.4 + 0.0j, radius=120, l_max=3, m_max=2)
    sphere3 = part.Sphere(position=[-2000, 2000, 90], refractive_index=2.5 + 0.0j, radius=90, l_max=2, m_max=2)
    part_list = [sphere1, sphere2, sphere3]

    # water over aluminum
    lay_sys = lay.LayerSystem([0, 0], [1+6j, 1])

    # initial field
    plane_wave = init.PlaneWave(vacuum_wavelength=ld, polar_angle=np.pi, azimuthal_angle=0, polarization=0)

    # run simulation
    simulation = simul.Simulation(layer_system=lay_sys, particle_list=part_list, initial_field=plane_wave,
                                  neff_resolution=5e-4)
    simulation.run()
    
    scat_pwe, _ = sf.scattered_field_pwe(ld, part_list, lay_sys, 1)
       
    x = np.array(200.0)
    y = np.array(-140.0)
    z = np.array(460.0)

    ex, ey, ez = scat_pwe.electric_field(x, y, z)
    
    scat_pwe.set_reference_point([300, 260, -100])
    ex2, ey2, ez2 = scat_pwe.electric_field(x, y, z)
    
    absolute_error = np.linalg.norm([ex - ex2, ey - ey2, ez - ez2])
    relative_error = absolute_error / np.linalg.norm([ex, ey, ez])
    print("relative error:", relative_error)
    assert relative_error < 5e-3
    print("OK.")


if __name__ == '__main__':
    test_pwe_translation()

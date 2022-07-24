import smuthi.initial_field as init
import smuthi.particles as part
import smuthi.simulation as simul
import smuthi.layers as lay
import smuthi.postprocessing.scattered_field as sf
import smuthi.postprocessing.far_field as ff
import smuthi.fields as flds
import numpy as np
import matplotlib.pyplot as plt


def test_ff_power_consistency():
    """
    With this routine we test equation (3.14) of Bohren & Huffman, i.e.,
    the far field intensity
    I = k / (2 omega mu) |A|^2 / k^2
      = |A|^2 / (n k^2)
    """

    ld = 550
    n0 = 1.52
    n1 = 1.33
    omega = flds.angular_frequency(ld)
    k = omega * n1

    # initialize particles
    sphere1 = part.Sphere(position=[0, 0, 110], refractive_index=2.4 + 0.0j,
                          radius=110, l_max=1, m_max=1)
    sphere2 = part.Sphere(position=[-2000, -2000, 120],
                          refractive_index=2.4 + 0.0j, radius=120, l_max=5,
                          m_max=4)
    sphere3 = part.Sphere(position=[-2000, 2000, 90],
                          refractive_index=2.5 + 0.0j, radius=90, l_max=2,
                          m_max=2)
    part_list = [sphere1, sphere2, sphere3]

    # water over aluminum
    lay_sys = lay.LayerSystem([0, 0], [n0, n1])

    # initial field
    plane_wave = init.PlaneWave(vacuum_wavelength=ld, polar_angle=np.pi * 4/5,
                                azimuthal_angle=0.1, polarization=0)

    # run simulation
    simulation = simul.Simulation(layer_system=lay_sys, particle_list=part_list,
                                  initial_field=plane_wave,
                                  neff_resolution=5e-4)
    simulation.run()

    scat_ff = ff.scattered_far_field(ld, part_list, lay_sys).top()

    iph = 18  # some arbitrary azimuthal angle value

    thetas_degree = []
    intens_te = []
    intens_tm = []
    intens_te_from_amplitude = []
    intens_tm_from_amplitude = []

    for ith in range(len(scat_ff.polar_angles)):

        polar = scat_ff.polar_angles[ith]
        azim = scat_ff.azimuthal_angles[iph]

        thetas_degree.append(polar * 180 / np.pi)

        intens_te.append(scat_ff.signal[0, ith, iph])
        intens_tm.append(scat_ff.signal[1, ith, iph])

        # alternatively, we compute the intensity according to (3.14) of Bohren&Huffman
        intens_te_from_amplitude.append(k / (2 * omega) * abs(scat_ff.amplitude[0, ith, iph])**2 / k**2)
        intens_tm_from_amplitude.append(k / (2 * omega) * abs(scat_ff.amplitude[1, ith, iph])**2 / k**2)

        relerr_te = abs(intens_te[-1] - intens_te_from_amplitude[-1]) / intens_te[-1]
        relerr_tm = abs(intens_tm[-1] - intens_tm_from_amplitude[-1]) / intens_tm[-1]

        print("===============================")
        print("polar angle: ", thetas_degree[-1])
        print("relative error TE: ", relerr_te)
        print("relative error TM: ", relerr_tm)

        assert relerr_te < 1e-8
        assert relerr_tm < 1e-8

        print ("OK.")

if __name__ == '__main__':
    test_ff_power_consistency()

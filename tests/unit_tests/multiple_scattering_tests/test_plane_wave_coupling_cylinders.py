# -*- coding: utf-8 -*-
import numpy as np
import smuthi.particles as part
import smuthi.linearsystem.particlecoupling.direct_coupling as pacou
import smuthi.layers as lay
import smuthi.fields as flds
import smuthi.initial_field as init
import smuthi.simulation

wl = 550
lay_sys = lay.LayerSystem([0, 800, 0], [1, 1, 1])
plane_wave = init.PlaneWave(vacuum_wavelength=wl, polar_angle=0, azimuthal_angle=0, polarization=0)

cylinder1 = part.FiniteCylinder(position=[0, 0, 200],  refractive_index=2.4 + 0.0j,
                                cylinder_height=100, cylinder_radius=50, l_max=5)
  
cylinder2 = part.FiniteCylinder(position=[300, 100, 200],  refractive_index=2.4 + 0.0j,
                                cylinder_height=100, cylinder_radius=50, l_max=5)


def test_W_block_pvwf_coupling():
    simul = smuthi.simulation.Simulation(layer_system=lay_sys,
                                         initial_field=plane_wave,
                                         particle_list=[cylinder1, cylinder2],
                                         use_pvwf_coupling=True,
                                         pvwf_coupling_neff_max=8,
                                         pvwf_coupling_neff_resolution=1e-3)

    simul.run()
    blocksize1 = flds.blocksize(cylinder1.l_max, cylinder1.m_max)
    blocksize2 = flds.blocksize(cylinder2.l_max, cylinder2.m_max)
    w = np.zeros([blocksize1 + blocksize2, blocksize1 + blocksize2])
    for ii in range(blocksize1 + blocksize2):
        vec = np.zeros([blocksize1+blocksize2])
        vec[ii] = 1
        w[:, ii] = simul.linear_system.coupling_matrix.linear_operator.matvec(vec)

    simul = smuthi.simulation.Simulation(layer_system=lay_sys,
                                         initial_field=plane_wave,
                                         particle_list=[cylinder1, cylinder2])

    simul.run()
    w2 = np.zeros([blocksize1 + blocksize2, blocksize1 + blocksize2])
    for ii in range(blocksize1 + blocksize2):
        vec = np.zeros([blocksize1+blocksize2])
        vec[ii] = 1
        w2[:, ii] = simul.linear_system.coupling_matrix.linear_operator.matvec(vec)

    diff = np.linalg.norm(w - w2)
    print("Relative difference: ", diff / np.linalg.norm(w2))
    assert diff / np.linalg.norm(w2) < 1e-4


if __name__ == '__main__':
    test_W_block_pvwf_coupling()

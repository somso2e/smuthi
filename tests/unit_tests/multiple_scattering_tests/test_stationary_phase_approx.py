# -*- coding: utf-8 -*-
"""Test the stationary phase approximation for the layer mediated particle coupling"""

import numpy as np
from smuthi.linearsystem.particlecoupling.layer_mediated_coupling \
    import layer_mediated_coupling_block, layer_mediated_coupling_block_stat_phase_approx
import smuthi.layers as lay
import smuthi.particles as part
import smuthi.fields as flds

wl = 550
rho = 200000
phi = np.pi / 3

# dummy particles
emitter = part.Sphere(position=[0, 0, 200], l_max=4, m_max=3)
receiver = part.Sphere(position=[rho * np.cos(phi), rho * np.sin(phi), 300], l_max=2, m_max=2)

# dielectric substrate in vacuum
lsys = lay.LayerSystem(thicknesses=[0, 0], refractive_indices=[1.5, 1])

# high accuracy Sommerfeld contour for reference
kpar = flds.reasonable_Sommerfeld_kpar_contour(vacuum_wavelength=wl, layer_refractive_indices=lsys.refractive_indices,
                                               neff_imag=1e-4, neff_resolution=1e-4, neff_max=4)


def test_stationary_phase_wr_against_numeric_integral():
    wr0 = layer_mediated_coupling_block_stat_phase_approx(vacuum_wavelength=wl,
                                                          receiving_particle=receiver,
                                                          emitting_particle=emitter,
                                                          layer_system=lsys)

    wr_check = layer_mediated_coupling_block(vacuum_wavelength=wl,
                                             receiving_particle=receiver,
                                             emitting_particle=emitter,
                                             layer_system=lsys,
                                             k_parallel=kpar)

    rel_error = np.linalg.norm(wr0 - wr_check) / np.linalg.norm(wr_check)

    print("Stationary phase approximation: Relative deviation = ", rel_error)
    assert rel_error < 1e-2


if __name__ == '__main__':
    test_stationary_phase_wr_against_numeric_integral()

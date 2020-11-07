import numpy as np
from matplotlib import pyplot as plt
import smuthi.linearsystem.particlecoupling.layer_mediated_coupling as pc
import smuthi.particles as part
import smuthi.layers as lay
import smuthi.fields as flds

# vacuum wavelength
wl = 550

# dummy particles
spher1 = part.Sphere(position=[0, 0, 100], l_max=1, m_max=1)
spher2 = part.Sphere(position=[1000, 0, 100], l_max=1, m_max=1)

# layer system (ambient = air, substrate = glass)
laysys = lay.LayerSystem(thicknesses=[0, 0], refractive_indices=[1.5, 1])

# lateral distances
rho_array = 10**np.linspace(0, 5.7, 100)

def compute_coupling(neff_resolution, neff_imag):
    kpar = flds.reasonable_Sommerfeld_kpar_contour(550, layer_refractive_indices=laysys.refractive_indices,
                                                   neff_max=2.5, neff_resolution=neff_resolution, neff_imag=neff_imag)
    si_list = []
    for rho in rho_array:
        spher2.position = [rho, 0, 100]
        wr = pc.layer_mediated_coupling_block(vacuum_wavelength=wl,
                                              receiving_particle=spher2,
                                              emitting_particle=spher1,
                                              layer_system=laysys,
                                              k_parallel=kpar)
        si_list.append(np.linalg.norm(wr))
    return si_list

# show results
plt.figure

# accurate reference result (high resolution, moderate deflection into imaginary)
plt.loglog(rho_array, compute_coupling(neff_resolution=2e-4, neff_imag=1e-3), label="neff_resol=2e-4, neff_imag=1e-3")
plt.loglog(rho_array, compute_coupling(neff_resolution=2e-3, neff_imag=1e-3), label="neff_resol=2e-3, neff_imag=1e-3 (aliasing)")
plt.loglog(rho_array, compute_coupling(neff_resolution=2e-4, neff_imag=1e-2), label="neff_resol=2e-4, neff_imag=1e-2 (complex Bessel problem)")
plt.legend()
plt.xlabel("lateral distance (nm)")
plt.ylabel("substrate-mediated particle coupling (a.u.)")
plt.show()

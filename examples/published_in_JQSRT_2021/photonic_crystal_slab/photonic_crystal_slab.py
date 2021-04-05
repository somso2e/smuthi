#*****************************************************************************#
# This is a script for Smuthi v1.1.0                                          #
# It computes the resonant Purcell factor enhancement of a dipole placed in   #
# the center of an L3 photonic crystal slab cavity                            # 
#*****************************************************************************#

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from smuthi.simulation import Simulation
from smuthi.initial_field import DipoleSource
from smuthi.layers import LayerSystem
from smuthi.particles import FiniteCylinder
from smuthi.postprocessing.graphical_output import show_near_field, plot_particles
from smuthi.utility.automatic_parameter_selection import select_numerical_parameters
from smuthi.utility import cuda
import smuthi.linearsystem.linear_system as linsys

cuda.enable_gpu()

# all lengths are given in micrometers
lattice_constant = 0.42
slab_thickness = 0.6 * lattice_constant
hole_radius = 0.29 * lattice_constant
silicon_ri = 3.473
thickness = 7      # the cavity walls are that many holes thick
shift = 0.1        # see https://doi.org/10.1038/nature02063
hexagonal = False  # if true, cut corners off the pattern
run_purcell = False
show_field_distribution = True

if hexagonal:
    outputdir = "results/thickness%i_hex"%thickness
else:
    outputdir = "results/thickness%i_nohex"%thickness

if not os.path.isdir(outputdir):
    os.makedirs(outputdir)


def set_simulation(wavelength=1.6, lmax=5, mmax=2, thickness=7, neff_max=5.67,
                   test_balloon=False, hexagonal=True, shift=0,
                   tolerance=1e-3):
    holes_list = []
    if test_balloon:
      x = 0.5 * lattice_constant
      y = lattice_constant * np.sqrt(3) / 2
      hole = FiniteCylinder(position=[x, y, slab_thickness / 2],
                            refractive_index=1.0,
                            cylinder_height=slab_thickness,
                            cylinder_radius=hole_radius,
                            l_max=lmax,
                            m_max=mmax)
      holes_list.append(hole)
      x = 1.5 * lattice_constant
      hole = FiniteCylinder(position=[x, y, slab_thickness / 2],
                            refractive_index=1.0,
                            cylinder_height=slab_thickness,
                            cylinder_radius=hole_radius,
                            l_max=lmax,
                            m_max=mmax)
      holes_list.append(hole)
    elif hexagonal:
      for iy in range(-thickness, thickness+1):
        ixmax = 1 + thickness - abs(iy)
        for ix in range(-thickness-1, ixmax+1):            
          if (not (iy == 0 and ix in [-2, -1, 0, 1, 2])):
            x = (ix + abs(iy) / 2) * lattice_constant
            y = iy * lattice_constant * np.sqrt(3)/2
            hole = FiniteCylinder(position=[x, y, slab_thickness/2],
                                  refractive_index=1.0,
                                  cylinder_height=slab_thickness,
                                  cylinder_radius=hole_radius,
                                  l_max=lmax,
                                  m_max=mmax)
            holes_list.append(hole)                                     
    else:
      for iy in range(-thickness, thickness+1):        
        if iy%2:
          xrange = range(-thickness-1, thickness+1)
        else:
          xrange = range(-thickness, thickness+1)
        for ix in xrange:            
          if (not (iy == 0 and ix in [-2, -1, 0, 1, 2])):
            x = (ix + iy%2/2) * lattice_constant
            y = iy * lattice_constant * np.sqrt(3)/2
            hole = FiniteCylinder(position=[x, y, slab_thickness/2],
                                  refractive_index=1.0,
                                  cylinder_height=slab_thickness,
                                  cylinder_radius=hole_radius,
                                  l_max=lmax,
                                  m_max=mmax)
            holes_list.append(hole)

    if not test_balloon:
      x = - (2 + shift) * lattice_constant
      hole = FiniteCylinder(position=[x, 0, slab_thickness/2],
                            refractive_index=1.0,
                            cylinder_height=slab_thickness,
                            cylinder_radius=hole_radius,
                            l_max=lmax,
                            m_max=mmax)
      holes_list.append(hole)
      x = (2 + shift) * lattice_constant
      hole = FiniteCylinder(position=[x, 0, slab_thickness/2],
                            refractive_index=1.0,
                            cylinder_height=slab_thickness,
                            cylinder_radius=hole_radius,
                            l_max=lmax,
                            m_max=mmax)
      holes_list.append(hole)


    layers = LayerSystem(thicknesses=[0, slab_thickness, 0],
                         refractive_indices=[1, silicon_ri, 1])

    dipole = DipoleSource(vacuum_wavelength=wavelength,
                          dipole_moment=[0, 1j, 0],
                          position=[0, 0, slab_thickness / 2])

    return Simulation(layer_system=layers,
                      particle_list=holes_list,
                      initial_field=dipole,
                      length_unit='micron',
                      neff_max=neff_max,
                      neff_resolution=1e-3,
                      neff_imag=1e-3,
                      solver_type='GCROTMK',
                      solver_tolerance=tolerance,
                      store_coupling_matrix=False,
                      coupling_matrix_lookup_resolution=0.005,
                      log_to_file=True,
                      identical_particles=True)


def visualize_geometry(thickness=7, hexagonal=True, shift=0):
    demo_simulation = set_simulation(thickness=thickness, hexagonal=False,
                                     shift=0.15)
    print("Number of particles:", len(demo_simulation.particle_list))
    plt.figure()
    xmax=4
    plot_particles(xmin=-xmax, xmax=xmax, ymin=-xmax, ymax=xmax, zmin=0, zmax=0,
                   particle_list=demo_simulation.particle_list,
                   draw_circumscribing_sphere=False, fill_particle=True)
    plt.gca().set_aspect('equal')
    plt.xlim([-xmax,xmax])
    plt.ylim([-xmax,xmax])
    plt.xlabel("x [micron]")
    plt.ylabel("y [micron]")
    plt.savefig("thickness%i.png"%thickness)
    plt.show()


def purcell_factor(simulation):
    dip = simulation.initial_field
    diss_pow = dip.dissipated_power(particle_list=simulation.particle_list,
                                    layer_system=simulation.layer_system)
    diss_pow_hom = dip.dissipated_power_homogeneous_background(simulation.layer_system)
    return diss_pow / diss_pow_hom


def get_autoparams():
    simulation = set_simulation(test_balloon=True)
    select_numerical_parameters(simulation,
                                detector=purcell_factor,
                                tolerance=0.001,
                                tolerance_factor=0.1,
                                max_iter=30,
                                neff_max_increment=0.5,
                                neff_max_offset=1,
                                show_plot=False)
    return (simulation.particle_list[0].l_max, 
            simulation.particle_list[0].m_max,
            simulation.neff_max)


def Lorentzian(x, amp, cen, wid):
    return amp * wid**2 / ((x - cen)**2 + wid**2)


shifts = [0, 0.05, 0.1, 0.15]
# expected outcome (from previous runs):
expected_center_wls = [1.60190480, 1.60687863, 1.61092022, 1.61401999]
expected_half_widths = [1.76294018e-04, 1.15822133e-04, 6.94904799e-05, 4.28346599e-05]

if run_purcell:
    for i, shift in enumerate(shifts):
        wls = [expected_center_wls[i]]
        for iwl in range(1, 6):
            dwl = (iwl / 5) * 3 * expected_half_widths[i]
            wls.append(expected_center_wls[i] - dwl)
            wls.append(expected_center_wls[i] + dwl)
        print("wls:", wls)
        wl_list = []
        purc_list = []
        tol_list = []
        for wl in wls:
            for tol in [1e-3, 2e-3, 3e-3, 4e-3, 5e-3]:
                print("running shift = ", shift, ", wl = ", wl, ", tol = ", tol)
                simulation = set_simulation(wavelength=wl,
                                            thickness=thickness,
                                            hexagonal=hexagonal,
                                            shift=shift,
                                            tolerance=tol)
                simulation.run()
                pf = purcell_factor(simulation)
                if 0 < pf < 1e5:
                    break
            wl_list.append(wl)
            purc_list.append(pf)
            tol_list.append(simulation.solver_tolerance)
            np.savetxt(outputdir + "/purcell_shift%i.txt"%(100*shift),
                       np.transpose([wl_list, purc_list, tol_list]))


if show_field_distribution:
    simulation = set_simulation(wavelength=expected_center_wls[0], hexagonal=hexagonal)
    simulation = set_simulation(wavelength=expected_center_wls[0], 
                                thickness=thickness, 
                                hexagonal=hexagonal, 
                                shift=0)
    simulation.run()
    vmax = 3e3
    show_near_field(simulation=simulation,
                    show_plots=False,
                    quantities_to_plot=["norm(E)", "E_y"],
                    show_opts=[{'vmin': -vmax, 'vmax': vmax},
                               {'vmin': 0, 'vmax': vmax},
                               {'vmin': -vmax, 'vmax': vmax}],
                    save_plots=True,
                    save_opts=[{'format': 'pdf', 'dpi': 200},
                               {'format': 'eps', 'dpi': 200},
                               {'format': 'gif'}],
                    outputdir=outputdir,
                    xmin=-3.5,
                    xmax=3.5,
                    ymin=-3,
                    ymax=3,
                    zmin=slab_thickness * 0.6,
                    zmax=slab_thickness * 0.6,
                    resolution_step=0.02,
                    show_internal_field=False,
                    save_data=True,
                    data_format='ascii')
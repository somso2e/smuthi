#*****************************************************************************#
# In this Smuthi example script the power flow of a plane wave through a      #
# 200nm thick glass sheet is visualized that is perforated by a periodic      #
# arrangement of nanoholes.                                                   #
# The transmittance of the plane wave strongly depends on the ambient medium. #
# As the refractive index of the ambient medium is increased, the power flow  #
# localizes into the sheet's nanoholes.                                       #
# The script runs with Smuthi version 1.2                                     #
#*****************************************************************************#

import smuthi.simulation as sim
import smuthi.initial_field as init
import smuthi.layers as lay
import smuthi.particles as part
import smuthi.fields as flds
import smuthi.periodicboundaries as pb
import smuthi.periodicboundaries.post_processing as pbpost
import smuthi.utility.cuda
import numpy as np
from tqdm import tqdm
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import tempfile
import shutil
import imageio

# try to enable GPU calculations
smuthi.utility.cuda.enable_gpu()

# in this file, all lengths are given in nanometers
wl = 500

# set up parameter sweep
refractive_indices = np.arange(1, 2.5 + 0.1, 0.1)

neffmax = 3
neffimag = 0.01
waypoints = [0, 0.8, 0.8 - 1j * neffimag, 2.6 - 1j * neffimag, 2.6, neffmax]
neff_discr = 1e-3  
flds.default_Sommerfeld_k_parallel_array = flds.reasonable_Sommerfeld_kpar_contour(vacuum_wavelength=wl,
                                                                                   neff_waypoints=waypoints,
                                                                                   neff_resolution=neff_discr)
flds.angular_arrays(angular_resolution=np.pi/360)


# define initial field
initial_field = init.PlaneWave(vacuum_wavelength=wl,
                               polar_angle=0,
                               azimuthal_angle=3/2*np.pi,
                               polarization=0,
                               amplitude=1,
                               reference_point=[0, 0, 0])

# define unit cell
a1 = np.array([320, 0, 0], dtype=float)
a2 = np.array([160, 320, 0], dtype=float) 
pb.default_Ewald_sum_separation = pb.set_ewald_sum_separation(a1, a2)


images = []
tempdir = tempfile.mkdtemp()
for idx, rfidx in enumerate(tqdm(refractive_indices, desc='Refractive index iterations  ', file=sys.stdout,
                            bar_format='{l_bar}{bar}| elapsed: {elapsed} ' 'remaining: {remaining}')):
    # suppress consol output
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    old_stderr = sys.stderr
    sys.stderr = open(os.devnull, 'w')
    
    # define layer system
    rfidx_glass = 1.5
    layer_system = lay.LayerSystem([0, 200, 0], [rfidx, rfidx_glass, rfidx])

    # define particle list
    particle_list = [part.FiniteCylinder(position=[0, 0, 100], refractive_index=rfidx,
                                         cylinder_radius=100, cylinder_height=200,
                                         l_max=5, m_max=5)]
    
    # run simulation
    simulation = sim.Simulation(initial_field=initial_field,
                                layer_system=layer_system,
                                particle_list=particle_list,
                                periodicity=(a1, a2),
                                ewald_sum_separation_parameter='default',
                                number_of_threads_periodic=-2) # all but 2 threads
    simulation.run()
    
    # plane wave expansion of total transmitted field
    pwe_total_T = pbpost.transmitted_plane_wave_expansion(initial_field,
                                                          particle_list,
                                                          layer_system,
                                                          a1, a2)
    
    
    # specify field points for which the electric and magnetic field is evaluated
    # in this case 5x5 unit cells in the x-y-plane are defined
    discr = 3
    x = np.arange(-2.5 * a1[0], 2.5 * a1[0] + discr, discr)
    y = np.arange(-2.5 * a2[1], 2.5 * a2[1] + discr, discr)
    xgrid, ygrid = np.meshgrid(x, y)
    zgrid = xgrid - xgrid + np.sum(layer_system.thicknesses) + 0
    

    # evaluate the electric and magnetic field
    ex, ey, ez = pbpost.electromagnetic_nearfield(pwe_total_T, xgrid, ygrid, zgrid,
                    field_type='electric')
    hx, hy, hz = pbpost.electromagnetic_nearfield(pwe_total_T, xgrid, ygrid, zgrid,
                    field_type='magnetic', vacuum_wavelength=simulation.initial_field.vacuum_wavelength)
    # z-component of the Poynting vector
    Sz = pbpost.conjugated_poynting_vector((ex, ey, ez), (hx, hy, hz))[2]
    # power per area of the intial plane wave (to normalize the Poynting vector)
    P_init = pbpost.initial_plane_wave_power_per_area(simulation.initial_field,
                                                      simulation.layer_system)
    
    
    # Creat plots
    tempfig, ax = plt.subplots()
    plt.title('z-component of Poynting vector')
    plt.imshow(abs(Sz)/P_init, cmap='inferno', origin='lower', vmin=0, vmax=2,
              extent=[x[0], x[-1], y[0], y[-1]])

    # highlight central unit cell
    ax.plot([-240, 80], [-160, -160], color='w')
    ax.plot([-240, -80], [-160, 160], color='w')
    ax.plot([80, 240], [-160, 160], color='w')
    ax.plot([-80, 240], [160, 160], color='w')
    ax.text(x=-100, y=-220, s='$a_1$', color='w')
    ax.text(x=-270, y=-10, s='$a_2$', color='w')
    plt.text(x=-750, y=450, s='$n_{\mathrm{amb}}=%.1f$' % rfidx, color='w')

    ax.xaxis.set_major_locator(MultipleLocator(320))
    ax.set_xticklabels(['-3', '-2', '-1', '0', '1', '2'])
    ax.xaxis.set_minor_locator(MultipleLocator(160))
    plt.xlabel('x / P')
    ax.yaxis.set_major_locator(MultipleLocator(320))
    ax.set_yticklabels(['-3', '-2', '-1', '0', '1', '2'])
    ax.yaxis.set_minor_locator(MultipleLocator(160))
    plt.ylabel('y / P')
    plt.colorbar()
    
    # save gif
    tempfig_filename = tempdir + '/temp_rfidx_' + str(rfidx) + '.png'
    plt.savefig(tempfig_filename, dpi=200, bbox_inches='tight')
    plt.close(tempfig)
    images.append(imageio.imread(tempfig_filename))
export_filename = 'poynting_vector_z'
imageio.mimsave(export_filename + '.gif', images, duration=0.3)
shutil.rmtree(tempdir)    





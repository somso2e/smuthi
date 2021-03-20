"""Functions to generate plots and animations."""

# -*- coding: utf-8 -*-
import numpy as np
import smuthi.postprocessing.scattered_field as sf
import smuthi.postprocessing.internal_field as intf
import smuthi.postprocessing.far_field as ff
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Rectangle
from matplotlib.colors import Normalize, LogNorm, SymLogNorm
from itertools import cycle
import tempfile
import shutil
import imageio
import os
import warnings
import sys
from tqdm import tqdm


def plot_layer_interfaces(dim1min, dim1max, layer_system):
    """Add lines to plot to display layer system interfaces

    Args:
        dim1min (float):    From what x-value plot line
        dim1max (float):    To what x-value plot line
        layer_system (smuthi.layers.LayerSystem):   Stratified medium
    """
    for il in range(1, layer_system.number_of_layers()):
        plt.plot([dim1min, dim1max], [layer_system.reference_z(il), layer_system.reference_z(il)], 'g')


def plot_particles(xmin, xmax, ymin, ymax, zmin, zmax, particle_list,
                   draw_circumscribing_sphere, fill_particle=True):
    """Add circles, ellipses and rectangles to plot to display spheres, spheroids and cylinders.

    Args:
        xmin (float):   Minimal x-value of plot
        xmax (float):   Maximal x-value of plot
        ymin (float):   Minimal y-value of plot
        ymax (float):   Maximal y-value of plot
        zmin (float):   Minimal z-value of plot
        zmax (float):   Maximal z-value of plot
        particle_list (list): List of smuthi.particles.Particle objects
        draw_circumscribing_sphere (bool): If true (default), draw a circle indicating the circumscribing sphere of
                                           particles.
        fill_particle (bool):              If true, draw opaque particles.
    """

    ax = plt.gca()

    if fill_particle:
        particle_face_color = 'w'
    else:
        particle_face_color = 'none'

    if xmin == xmax:
        plane_coord = 0
        draw_coord = [1, 2]
    elif ymin == ymax:
        plane_coord = 1
        draw_coord = [0, 2]
    elif zmin == zmax:
        plane_coord = 2
        draw_coord = [0, 1]
    else:
        raise ValueError('Field points must define a plane')

    for particle in particle_list:
        pos = particle.position
        dis = abs((xmin, ymin, zmin)[plane_coord] - pos[plane_coord])

        if dis > particle.circumscribing_sphere_radius():
            continue

        circradius = np.sqrt(particle.circumscribing_sphere_radius()**2 - dis**2)
        if type(particle).__name__ == 'Sphere':
            ax.add_patch(Circle((pos[draw_coord[0]], pos[draw_coord[1]]), circradius,
                         facecolor=particle_face_color, edgecolor='k'))
        else:
            if not particle.euler_angles == [0, 0, 0]:
                warnings.warn("Drawing rotated particles currently not supported - drawing black disc with size"
                              + " based on the circumscribing sphere instead")
                ax.add_patch(Circle((pos[draw_coord[0]], pos[draw_coord[1]]), circradius, facecolor='k', edgecolor='k'))
                ax.text(pos[draw_coord[0]], pos[draw_coord[1]], 'rotated ' + type(particle).__name__,
                        verticalalignment='center', horizontalalignment='center', color='blue', fontsize=5)
            else:
                if draw_circumscribing_sphere:
                    ax.add_patch(Circle((pos[draw_coord[0]], pos[draw_coord[1]]), circradius,
                                 linestyle='dashed', facecolor='none', edgecolor='k'))

                if type(particle).__name__ == 'Spheroid':
                    a = particle.semi_axis_a
                    c = particle.semi_axis_c
                    if plane_coord in [0, 1]:
                        w = a * np.sqrt(1 - (dis/a)**2) if dis < a else 0
                        h = c * np.sqrt(1 - (dis/a)**2) if dis < a else 0
                    else:
                        w = a * np.sqrt(1 - (dis/c)**2) if dis < c else 0
                        h = w
                    ax.add_patch(Ellipse(xy=(pos[draw_coord[0]], pos[draw_coord[1]]), width=2*w, height=2*h,
                                         facecolor='w', edgecolor='k'))

                elif type(particle).__name__ == 'FiniteCylinder':
                    r = particle.cylinder_radius
                    h = particle.cylinder_height
                    if plane_coord in [0, 1]:
                        w = np.sqrt(r**2 - dis**2) if dis < r else 0
                        h *= np.sign(w)
                        ax.add_patch(Rectangle((pos[draw_coord[0]]-w, pos[draw_coord[1]]-h/2),
                                                2*w, h, facecolor='w', edgecolor='k'))
                    elif np.logical_and(plane_coord == 2, dis <= h/2):
                        ax.add_patch(Circle((pos[draw_coord[0]], pos[draw_coord[1]]), r,
                                            facecolor='w', edgecolor='k'))


def compute_near_field(simulation=None, X=None, Y=None, Z=None, type='scatt', chunksize=4096,
                       k_parallel='default', azimuthal_angles='default', angular_resolution=None):
    """Compute a certain component of the electric near field"""
    X, Y, Z = np.atleast_1d(X).T, np.atleast_1d(Y).T, np.atleast_1d(Z).T
    e_x, e_y, e_z = (np.zeros_like(X, dtype=np.complex128) for i in range(3))

    if 'sca' in type:
        min_laynum = simulation.layer_system.layer_number(Z.min())
        max_laynum = simulation.layer_system.layer_number(Z.max())
        layer_numbers = [i for i in range(min_laynum, max_laynum + 1)]
        scat_fld_exp = sf.scattered_field_piecewise_expansion(simulation.initial_field.vacuum_wavelength,
                                                              simulation.particle_list,
                                                              simulation.layer_system,
                                                              k_parallel,
                                                              azimuthal_angles,
                                                              angular_resolution,
                                                              layer_numbers)
    if 'int' in type:
        int_fld_exp = intf.internal_field_piecewise_expansion(simulation.initial_field.vacuum_wavelength, simulation.particle_list,
                                                              simulation.layer_system)

    descr = 'Compute '+type+' near-field'

    nchunks = np.ceil(X.size / chunksize).astype(int)
    e_x, e_y, e_z = ([None] * nchunks for i in range(3))

    Xc, Yc, Zc = np.array_split(X, nchunks), np.array_split(Y, nchunks), np.array_split(Z, nchunks)

    for c in tqdm(range(len(Xc)), desc=descr.ljust(26), file=sys.stdout,
                                  bar_format='{l_bar}{bar}| elapsed: {elapsed} ' 'remaining: {remaining}'):
        xarr, yarr, zarr = Xc[c], Yc[c], Zc[c]
        if 'sca' in type:
            e_x[c], e_y[c], e_z[c] = scat_fld_exp.electric_field(xarr, yarr, zarr)
        elif 'ini' in type:
            e_x[c], e_y[c], e_z[c] = simulation.initial_field.electric_field(xarr, yarr, zarr, simulation.layer_system)
        elif 'int' in type:
            e_x[c], e_y[c], e_z[c] = int_fld_exp.electric_field(xarr, yarr, zarr)

    return np.concatenate(e_x), np.concatenate(e_y), np.concatenate(e_z)


def show_near_field(simulation=None, quantities_to_plot=None,
                    show_plots=True, show_opts=None, save_plots=False, save_opts=None,
                    save_data=False, data_format='hdf5', outputdir='.',
                    xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0,
                    resolution_step=25, k_parallel='default', azimuthal_angles='default', angular_resolution=None,
                    draw_circumscribing_sphere=True, show_internal_field=False):
    """Plot the electric near field along a plane. To plot along the xy-plane, specify zmin=zmax and so on.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object
        quantities_to_plot:     List of strings that specify what to plot. Select from 'E_x', 'E_y', 'E_z', 'norm(E)'
                                The list may contain one or more of the following strings:

                                'E_x'           real part of x-component of complex total electric field
                                'E_y'           real part of y-component of complex total electric field
                                'E_z'           real part of z-component of complex total electric field
                                'norm(E)'       norm of complex total electric field

                                'E_scat_x'      real part of x-component of complex scattered electric field
                                'E_scat_y'      real part of y-component of complex scattered electric field
                                'E_scat_z'      real part of z-component of complex scattered electric field
                                'norm(E_scat)'  norm of complex scattered electric field

                                'E_init_x'      real part of x-component of complex initial electric field
                                'E_init_y'      real part of y-component of complex initial electric field
                                'E_init_z'      real part of z-component of complex initial electric field
                                'norm(E_init)'  norm of complex initial electric field
        show_plots (logical):   If True, plots are shown
        show_opts (dict list):  List of dictionaries containing options to be passed to imshow for plotting.
                                For each entry in quantities_to_plot, all show_opts dictionaries will be applied.
                                If save_plots=True, a 1:1 correspondence between show_opts and save_opts dictionaries
                                is assumed. For simplicity, one can also provide a single show_opts entry that will
                                be applied to all save_opts.
                                The following keys are made available (see matplotlib.pyplot.imshow documentation):
                                'cmap'          defaults to 'inferno' for norm quantities and 'RdYlBu' otherwise
                                'norm'          (None). If a norm is provided, its vmin and vmax take precedence
                                'aspect'        ('equal')
                                'interpolation' (None), also available: bilinear, bicubic, spline16, quadric, ...
                                'alpha'         (None)
                                'vmin'          (None), will be set to 0 for norm quantities and -vmax otherwise
                                'vmax'          initialized with the max of the quantity to plot
                                'origin'        ('lower')
                                'extent'        calculated automatically based on plotting coordinate limits
                                An optional extra key called 'label' of type string is shown in the plot title
                                and appended to the associated file if save_plots is True
                                Finally, an optional 'figsize' key is available to set the width and height of
                                the figure window (see matplotlib.pyplot.figure documentation)
        save_plots (logical):   If True, plots are exported to file.
        save_opts (dict list):  List of dictionaries containing options to be passed to savefig.
                                For each entry in quantities_to_plot, all save_opts dictionaries will be applied.
                                A 1:1 correspondence between save_opts and show_opts dictionaries is assumed. For
                                simplicity, one can also provide a single save_opts entry that will be applied to
                                all show_opts.
                                The following keys are made available (see matplotlib.pyplot.savefig documentation):
                                'dpi'           (None)
                                'orientation'   (None)
                                'format'        ('png'), also available: eps, jpeg, jpg, pdf, ps, svg, tif, tiff ...
                                'transparent'   (False)
                                'bbox_inches'   ('tight')
                                'pad_inches'    (0.1)
                                Passing 'gif' as one of the format values will result in an animation if the
                                quantity to plot is of non-norm type
        save_data (logical):    If True, raw data are exported to file
        data_format (str):      Output data format string, 'hdf5' and 'ascii' formats are available
        outputdir (str):        Path to directory where to save the export files
        xmin (float):           Plot from that x (length unit)
        xmax (float):           Plot up to that x (length unit)
        ymin (float):           Plot from that y (length unit)
        ymax (float):           Plot up to that y (length unit)
        zmin (float):           Plot from that z (length unit)
        zmax (float):           Plot up to that z (length unit)
        resolution_step (float):     Compute the field with that spatial resolution (length unit,
                                     distance between computed points), can be a tuple for [resx, resy, resz]
        k_parallel (numpy.ndarray or str):      in-plane wavenumbers for the plane wave expansion
                                                if 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
        azimuthal_angles (numpy.ndarray or str):azimuthal angles for the plane wave expansion
                                                if 'default', use smuthi.fields.default_azimuthal_angles
        angular_resolution (float):             If provided, angular arrays are generated with this angular
                                                resolution over the default angular range
        draw_circumscribing_sphere (bool):      If true (default), draw a circle indicating the circumscribing
                                                sphere of particles.
        show_internal_field (bool):             If true, compute also the field inside the particles (only for spheres)
    """

    if (not os.path.exists(outputdir)) and (save_plots or save_data):
        os.makedirs(outputdir)

    if quantities_to_plot is None:
        quantities_to_plot = ['norm(E)']
    if show_opts is None:
        show_opts = [{'interpolation':'none'}]
    if save_opts is None:
        save_opts = [{'format':'png'}]

    res = np.resize(resolution_step,3)
    x = np.arange(xmin, xmax + res[0]/2, res[0])
    y = np.arange(ymin, ymax + res[1]/2, res[1])
    z = np.arange(zmin, zmax + res[2]/2, res[2])
    xarr, yarr, zarr = np.squeeze(np.meshgrid(x, y, z, indexing='ij'))

    if xmin == xmax:
        dim1vec, dim2vec = y, z
        step1, step2 = res[1]/2, res[2]/2 # these are for proper registration of particles and layers to the pixel grid
        titlestr = '$x$ = {} '.format(xmin) + simulation.length_unit
        dim1name = '$y$ (' + simulation.length_unit + ')'
        dim2name = '$z$ (' + simulation.length_unit + ')'
    elif ymin == ymax:
        dim1vec, dim2vec = x, z
        step1, step2 = res[0]/2, res[2]/2
        titlestr = '$y$ = {} '.format(ymin) + simulation.length_unit
        dim1name = '$x$ (' + simulation.length_unit + ')'
        dim2name = '$z$ (' + simulation.length_unit + ')'
    elif zmin == zmax:
        dim1vec, dim2vec = x, y
        step1, step2 = res[0]/2, res[1]/2
        titlestr = '$z$ = {} '.format(zmin) + simulation.length_unit
        dim1name = '$x$ (' + simulation.length_unit + ')'
        dim2name = '$y$ (' + simulation.length_unit + ')'

    sys.stdout.write("Evaluate fields ...\n")
    sys.stdout.flush()
    e_x_scat, e_y_scat, e_z_scat = compute_near_field(simulation=simulation, X=xarr, Y=yarr, Z=zarr, type='scatt.')
    e_x_init, e_y_init, e_z_init = compute_near_field(simulation=simulation, X=xarr, Y=yarr, Z=zarr, type='initl.')

    if show_internal_field:
        e_x_int, e_y_int, e_z_int = compute_near_field(simulation=simulation, X=xarr, Y=yarr, Z=zarr, type='intrn.')

    sys.stdout.write("Generate plots ...\n")
    sys.stdout.flush()

    if not show_plots:
        import matplotlib
        default_backend = matplotlib.get_backend()
        matplotlib.use('Agg') # a non-gui backend to avoid opening figure windows

    for quantity in quantities_to_plot:

        for show_opt, save_opt in zip(cycle(show_opts), save_opts) if len(show_opts) < len(save_opts) else zip(show_opts, cycle(save_opts)):

            try:
                filename = 'E'
                label_str = show_opt.get('label','')

                if 'scat' in quantity:
                    e_x, e_y, e_z = e_x_scat, e_y_scat, e_z_scat
                    field_type_string = 'scat'
                elif 'init' in quantity:
                    e_x, e_y, e_z = e_x_init, e_y_init, e_z_init
                    field_type_string = 'init'
                elif 'intern' in quantity:
                    if not show_internal_field:
                        raise Exception("show_internal_field flag needs to be set to true!")
                    e_x, e_y, e_z = e_x_int, e_y_int, e_z_int
                    field_type_string = 'inter'
                else:
                    e_x, e_y, e_z = e_x_scat + e_x_init, e_y_scat + e_y_init, e_z_scat + e_z_init
                    if show_internal_field:
                        e_x, e_y, e_z = e_x + e_x_int, e_y + e_y_int, e_z + e_z_int
                    field_type_string = 'tot'

                fig = plt.figure(figsize=show_opt.get('figsize',[6.4, 4.8]))
                if not zmin == zmax:
                    plot_layer_interfaces(dim1vec[0], dim1vec[-1], simulation.layer_system)
                plot_particles(xmin, xmax, ymin, ymax, zmin, zmax, simulation.particle_list,
                               draw_circumscribing_sphere, not show_internal_field)
                if 'norm' in quantity:
                    e = np.sqrt(abs(e_x)**2 + abs(e_y)**2 + abs(e_z)**2)
                    vmax = np.abs(e).max()
                    color_norm = show_opt.get('norm', Normalize(vmin=show_opt.get('vmin',0), vmax=show_opt.get('vmax',vmax)))
                    plt.imshow(e,
                               alpha=show_opt.get('alpha'), norm=color_norm,
                               cmap=show_opt.get('cmap','inferno'), origin=show_opt.get('origin','lower'),
                               interpolation=show_opt.get('interpolation','none'),
                               extent=show_opt.get('extent', [dim1vec.min()-step1, dim1vec.max()+step1,
                                                              dim2vec.min()-step2, dim2vec.max()+step2]))
                    plt_title = '$|' + filename + '^{' + field_type_string + '}|$ at ' + titlestr \
                                + ('' if label_str == '' else ' (' + label_str + ')')
                    filename = 'norm_' + filename
                    plt.title(plt_title)
                else:
                    if '_x' in quantity:
                        e = e_x
                        filename = filename + '_x'
                    elif '_y' in quantity:
                        e = e_y
                        filename = filename + '_y'
                    elif '_z' in quantity:
                        e = e_z
                        filename = filename + '_z'
                    else:
                        print('Quantity:', quantity)
                        raise ValueError('field component not specified')

                    vmax = np.abs(e).max()
                    color_norm = show_opt.get('norm', Normalize(vmin=show_opt.get('vmin',-vmax), vmax=show_opt.get('vmax',vmax)))
                    plt.imshow(e.real,
                               alpha=show_opt.get('alpha'), norm=color_norm,
                               cmap=show_opt.get('cmap','RdYlBu'), origin=show_opt.get('origin','lower'),
                               interpolation=show_opt.get('interpolation','none'), aspect=show_opt.get('aspect','equal'),
                               extent=show_opt.get('extent', [dim1vec.min()-step1, dim1vec.max()+step1,
                                                              dim2vec.min()-step2, dim2vec.max()+step2]))
                    plt_title = '$ ' + filename + '^{' + field_type_string + '}$' + ' at ' + titlestr \
                                + ('' if label_str == '' else ' (' + label_str + ')')
                    plt.title(plt_title)

                plt.colorbar()
                plt.xlabel(dim1name)
                plt.ylabel(dim2name)

                label_str = '' if label_str == '' else '_' + label_str # prepend underscore
                filename = filename + '_' + field_type_string + label_str
                export_filename = outputdir + '/' + filename
                if save_plots:
                    if save_opt.get('format') == 'gif':
                        if 'norm' in quantity:
                            plt.close(fig)
                            continue # it seems like savefig does not accept 'gif' as a format, so we have to continue explicitly here
                        tempdir = tempfile.mkdtemp()
                        images = []
                        for i_t, t in enumerate(np.linspace(0, 1, 20, endpoint=False)):
                            tempfig = plt.figure(figsize=show_opt.get('figsize',[6.4, 4.8]))
                            if not zmin == zmax:
                                plot_layer_interfaces(dim1vec[0], dim1vec[-1], simulation.layer_system)
                            plot_particles(xmin, xmax, ymin, ymax, zmin, zmax, simulation.particle_list,
                                           draw_circumscribing_sphere, not show_internal_field)
                            e_t = e * np.exp(-1j * t * 2 * np.pi)
                            plt.imshow(e_t.real,
                                       alpha=show_opt.get('alpha'), norm=color_norm,
                                       cmap=show_opt.get('cmap','RdYlBu'), origin=show_opt.get('origin','lower'),
                                       interpolation=show_opt.get('interpolation','none'), aspect=show_opt.get('aspect','equal'),
                                       extent=show_opt.get('extent', [dim1vec.min()-step1, dim1vec.max()+step1,
                                                                      dim2vec.min()-step2, dim2vec.max()+step2]))
                            plt.title(plt_title)
                            plt.colorbar()
                            plt.xlabel(dim1name)
                            plt.ylabel(dim2name)

                            tempfig_filename = tempdir + '/temp_' + str(i_t) + '.png'
                            plt.savefig(tempfig_filename, dpi=save_opt.get('dpi'),
                                        orientation=save_opt.get('orientation','portrait'),
                                        transparent=save_opt.get('transparent',False),
                                        bbox_inches=save_opt.get('bbox_inches','tight'),
                                        pad_inches=save_opt.get('pad_inches',0.1))
                            plt.close(tempfig)
                            images.append(imageio.imread(tempfig_filename))
                        imageio.mimsave(export_filename + '.gif', images, duration=0.1)
                        shutil.rmtree(tempdir)
                    else:
                        file_ext = '.' + save_opt.get('format','png') # default to png if not specified
                        plt.savefig(export_filename + file_ext, dpi=save_opt.get('dpi'),
                                    orientation=save_opt.get('orientation','portrait'),
                                    transparent=save_opt.get('transparent',False),
                                    bbox_inches=save_opt.get('bbox_inches','tight'),
                                    pad_inches=save_opt.get('pad_inches',0.1))

                if not show_plots:
                    plt.close(fig)

            except Exception as e: # TODO: parse to understand if this is due to incompatible options or trying to plot non-2D data
                print("Skipping " + quantity + " for show_opt = ", show_opt, " and save_opt = ", save_opt)

    if not show_plots:
        matplotlib.use(default_backend) # now we can restore the original backend
    else:
        plt.show(block=False)

    if save_data:
        Ei = np.stack((e_x_init, e_y_init, e_z_init))
        Es = np.stack((e_x_scat, e_y_scat, e_z_scat))
        if show_internal_field:
            En = np.stack((e_x_int, e_y_int, e_z_int))

        if data_format.lower() == 'hdf5':
            import h5py
            # TODO: move metadata writing to auxiliary function
            # TODO: add metadata relative to particles and layer system ?
            if type(simulation.initial_field).__name__ == "DipoleSource":
                metadata = {'initial_field': type(simulation.initial_field).__name__,
                            'dipole_moment': simulation.initial_field.dipole_moment,
                            'dipole_position': simulation.initial_field.position,}
            elif type(simulation.initial_field).__name__ == "DipoleCollection":
                metadata = {'initial_field': type(simulation.initial_field).__name__,
                            'dipole_moments': [dip.dipole_moment for dip in simulation.initial_field.dipole_list],
                            'dipole_positions': [dip.position for dip in simulation.initial_field.dipole_list],}
            elif type(simulation.initial_field).__name__ == "PlaneWave":
                metadata = {'initial_field': type(simulation.initial_field).__name__,
                            'polar_angle': simulation.initial_field.polar_angle,
                            'azimuthal_angle': simulation.initial_field.azimuthal_angle,
                            'polarization': simulation.initial_field.polarization,}
            elif type(simulation.initial_field).__name__ == "GaussianBeam":
                metadata = {'initial_field': type(simulation.initial_field).__name__,
                            'polar_angle': simulation.initial_field.polar_angle,
                            'azimuthal_angle': simulation.initial_field.azimuthal_angle,
                            'polarization': simulation.initial_field.polarization,
                            'reference_point': simulation.initial_field.reference_point,
                            'beam_waist': simulation.initial_field.beam_waist,}

            with h5py.File(outputdir + '/data.hdf5', 'a') as f:
                # add attributes that are dependent on the initial field
                f.attrs.update(metadata)
                # additional attributes that are common to all simulations
                f.attrs['vacuum_wavelength'] = simulation.initial_field.vacuum_wavelength
                # write actual data
                g = f.require_group('near_field')
                g.create_dataset('x', data=x)
                g.create_dataset('y', data=y)
                g.create_dataset('z', data=z)
                g.create_dataset('init_electric_field', data=Ei, dtype='c16', compression="gzip")
                g.create_dataset('scat_electric_field', data=Es, dtype='c16', compression="gzip")
                if show_internal_field:
                    g.create_dataset('inte_electric_field', data=En, dtype='c16', compression="gzip")
        elif data_format.lower() == 'ascii':
            np.savetxt(outputdir + '/1st_dim_axis.dat', x, fmt='%g')
            np.savetxt(outputdir + '/2nd_dim_axis.dat', y, fmt='%g')
            np.savetxt(outputdir + '/3rd_dim_axis.dat', z, fmt='%g')
            fmt = list(np.repeat('%.15g %+.15gj', np.size(dim1vec)))
            np.savetxt(outputdir + '/e_init_x.dat', Ei[0,], fmt=fmt, delimiter=',')
            np.savetxt(outputdir + '/e_init_y.dat', Ei[1,], fmt=fmt, delimiter=',')
            np.savetxt(outputdir + '/e_init_z.dat', Ei[2,], fmt=fmt, delimiter=',')
            np.savetxt(outputdir + '/e_scat_x.dat', Es[0,], fmt=fmt, delimiter=',')
            np.savetxt(outputdir + '/e_scat_y.dat', Es[1,], fmt=fmt, delimiter=',')
            np.savetxt(outputdir + '/e_scat_z.dat', Es[2,], fmt=fmt, delimiter=',')
            if show_internal_field:
                np.savetxt(outputdir + '/e_int_x.dat', En[0,], fmt=fmt, delimiter=',')
                np.savetxt(outputdir + '/e_int_y.dat', En[1,], fmt=fmt, delimiter=',')
                np.savetxt(outputdir + '/e_int_z.dat', En[2,], fmt=fmt, delimiter=',')
        else:
            raise ValueError('Currently, only hdf5 or ascii output data formats are available')


def show_scattered_far_field(simulation, show_plots=True, show_opts=[{'label':'scattered_far_field'}],
                             save_plots=False, save_opts=None,
                             save_data=False, data_format='hdf5', outputdir='.',
                             flip_downward=True, split=True, log_scale=False,
                             polar_angles='default', azimuthal_angles='default', angular_resolution=None):
    """Display and export the scattered far field.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object
        show_plots (bool):                          Display plots if True
        show_opts (dict list):  List of dictionaries containing options to be passed to pcolormesh for plotting.
                                If save_plots=True, a 1:1 correspondence between show_opts and save_opts dictionaries
                                is assumed. For simplicity, one can also provide a single show_opts entry that will
                                be applied to all save_opts.
                                The following keys are available (see matplotlib.pyplot.pcolormesh documentation):
                                'alpha'     (None)
                                'cmap'      ('inferno')
                                'norm'      (None), is set to matplotlib.colors.LogNorm() if log_scale is True
                                'vmin'      (None), applies only to 2D plots
                                'vmax'      (None), applies only to 2D plots
                                'shading'   ('nearest'), applies only to 2D plots. 'gouraud' is also available
                                'linewidth' (None), applies only to 1D plots
                                'linestyle' (None), applies only to 1D plots
                                'marker'    (None), applies only to 1D plots
                                An optional extra key called 'label' of type string is shown in the plot title
                                and appended to the associated file if save_plots is True
        save_plots (bool):      If True, plots are exported to file.
        save_opts (dict list):  List of dictionaries containing options to be passed to savefig.
                                A 1:1 correspondence between save_opts and show_opts dictionaries is assumed. For
                                simplicity, one can also provide a single save_opts entry that will be applied to
                                all show_opts.
                                The following keys are made available (see matplotlib.pyplot.savefig documentation):
                                'dpi'           (None)
                                'orientation'   (None)
                                'format'        ('png'), also available: eps, jpeg, jpg, pdf, ps, svg, tif, tiff ...
                                'transparent'   (False)
                                'bbox_inches'   ('tight')
                                'pad_inches'    (0.1)
        save_data (bool):       If True, raw data are exported to file
        data_format (str):      Output data format string, 'hdf5' and 'ascii' formats are available
        outputdir (str):                        Path to the directory where files are to be saved
        flip_downward (bool):                   If True, represent downward directions as 0-90 deg instead of 90-180
        split (bool):                           If True, show two different plots for upward and downward directions
        log_scale (bool):                       If True, set a logarithmic scale
        polar_angles (numpy.ndarray or str):    Polar angles values (radian).
                                                If 'default', use smuthi.fields.default_polar_angles
        azimuthal_angles (numpy.ndarray or str):Azimuthal angle values (radian).
                                                If 'default', use smuthi.fields.default_azimuthal_angles
        angular_resolution (float):             If provided, angular arrays are generated with this angular resolution
                                                over the default angular range
    """

    infld = simulation.initial_field
    plst = simulation.particle_list
    lsys = simulation.layer_system
    far_field = ff.scattered_far_field(vacuum_wavelength=infld.vacuum_wavelength,
                                       particle_list=plst,
                                       layer_system=lsys,
                                       polar_angles=polar_angles,
                                       azimuthal_angles=azimuthal_angles,
                                       angular_resolution=angular_resolution)

    [d.setdefault('label','scattered_far_field') for d in show_opts]

    show_far_field(far_field=far_field, save_plots=save_plots, save_opts=save_opts, show_plots=show_plots,
                   show_opts=show_opts, save_data=save_data, data_format=data_format, outputdir=outputdir,
                   flip_downward=flip_downward, split=split, log_scale=log_scale)


def show_total_far_field(simulation, show_plots=True, show_opts=[{'label':'total_far_field'}],
                         save_plots=False, save_opts=None,
                         save_data=False, data_format='hdf5', outputdir='.',
                         flip_downward=True, split=True, log_scale=False,
                         polar_angles='default', azimuthal_angles='default', angular_resolution=None):
    """Display and export the total far field. This function cannot be used if the inital field is a plane wave.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object
        show_plots (bool):                          Display plots if True
        show_opts (dict list):  List of dictionaries containing options to be passed to pcolormesh for plotting.
                                If save_plots=True, a 1:1 correspondence between show_opts and save_opts dictionaries
                                is assumed. For simplicity, one can also provide a single show_opts entry that will
                                be applied to all save_opts.
                                The following keys are available (see matplotlib.pyplot.pcolormesh documentation):
                                'alpha'     (None)
                                'cmap'      ('inferno')
                                'norm'      (None), is set to matplotlib.colors.LogNorm() if log_scale is True
                                'vmin'      (None), applies only to 2D plots
                                'vmax'      (None), applies only to 2D plots
                                'shading'   ('nearest'), applies only to 2D plots. 'gouraud' is also available
                                'linewidth' (None), applies only to 1D plots
                                'linestyle' (None), applies only to 1D plots
                                'marker'    (None), applies only to 1D plots
                                An optional extra key called 'label' of type string is shown in the plot title
                                and appended to the associated file if save_plots is True
        save_plots (bool):      If True, plots are exported to file.
        save_opts (dict list):  List of dictionaries containing options to be passed to savefig.
                                A 1:1 correspondence between save_opts and show_opts dictionaries is assumed. For
                                simplicity, one can also provide a single save_opts entry that will be applied to
                                all show_opts.
                                The following keys are made available (see matplotlib.pyplot.savefig documentation):
                                'dpi'           (None)
                                'orientation'   (None)
                                'format'        ('png'), also available: eps, jpeg, jpg, pdf, ps, svg, tif, tiff ...
                                'transparent'   (False)
                                'bbox_inches'   ('tight')
                                'pad_inches'    (0.1)
        save_data (bool):       If True, raw data are exported to file
        data_format (str):      Output data format string, 'hdf5' and 'ascii' formats are available
        outputdir (str):                        Path to the directory where files are to be saved
        flip_downward (bool):                   If True, represent downward directions as 0-90 deg instead of 90-180
        split (bool):                           If True, show two different plots for upward and downward directions
        log_scale (bool):                       If True, set a logarithmic scale
        polar_angles (numpy.ndarray or str):    Polar angles values (radian).
                                                If 'default', use smuthi.fields.default_polar_angles
        azimuthal_angles (numpy.ndarray or str):Azimuthal angle values (radian).
                                                If 'default', use smuthi.fields.default_azimuthal_angles
        angular_resolution (float):             If provided, angular arrays are generated with this angular resolution
                                                over the default angular range
    """
    infld = simulation.initial_field
    plst = simulation.particle_list
    lsys = simulation.layer_system
    far_field, _, _ = ff.total_far_field(initial_field=infld,
                                         particle_list=plst,
                                         layer_system=lsys,
                                         polar_angles=polar_angles,
                                         azimuthal_angles=azimuthal_angles,
                                         angular_resolution=angular_resolution)

    [d.setdefault('label','total_far_field') for d in show_opts]

    show_far_field(far_field=far_field, save_plots=save_plots, save_opts=save_opts, show_plots=show_plots,
                   show_opts=show_opts, save_data=save_data, data_format=data_format, outputdir=outputdir,
                   flip_downward=flip_downward, split=split, log_scale=log_scale)


def show_scattering_cross_section(simulation, show_plots=True, show_opts=[{'label':'scattering_cross_section'}],
                                  save_plots=False, save_opts=None,
                                  save_data=False, data_format='hdf5', outputdir='.',
                                  flip_downward=True, split=True, log_scale=False,
                                  polar_angles='default', azimuthal_angles='default', angular_resolution=None):
    """Display and export the differential scattering cross section.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object
        show_plots (bool):                          Display plots if True
        show_opts (dict list):  List of dictionaries containing options to be passed to pcolormesh for plotting.
                                If save_plots=True, a 1:1 correspondence between show_opts and save_opts dictionaries
                                is assumed. For simplicity, one can also provide a single show_opts entry that will
                                be applied to all save_opts.
                                The following keys are available (see matplotlib.pyplot.pcolormesh documentation):
                                'alpha'     (None)
                                'cmap'      ('inferno')
                                'norm'      (None), is set to matplotlib.colors.LogNorm() if log_scale is True
                                'vmin'      (None), applies only to 2D plots
                                'vmax'      (None), applies only to 2D plots
                                'shading'   ('nearest'), applies only to 2D plots. 'gouraud' is also available
                                'linewidth' (None), applies only to 1D plots
                                'linestyle' (None), applies only to 1D plots
                                'marker'    (None), applies only to 1D plots
                                An optional extra key called 'label' of type string is shown in the plot title
                                and appended to the associated file if save_plots is True
        save_plots (bool):      If True, plots are exported to file.
        save_opts (dict list):  List of dictionaries containing options to be passed to savefig.
                                A 1:1 correspondence between save_opts and show_opts dictionaries is assumed. For
                                simplicity, one can also provide a single save_opts entry that will be applied to
                                all show_opts.
                                The following keys are made available (see matplotlib.pyplot.savefig documentation):
                                'dpi'           (None)
                                'orientation'   (None)
                                'format'        ('png'), also available: eps, jpeg, jpg, pdf, ps, svg, tif, tiff ...
                                'transparent'   (False)
                                'bbox_inches'   ('tight')
                                'pad_inches'    (0.1)
        save_data (bool):       If True, raw data are exported to file
        data_format (str):      Output data format string, 'hdf5' and 'ascii' formats are available
        outputdir (str):                        Path to the directory where files are to be saved
        flip_downward (bool):                   If True, represent downward directions as 0-90 deg instead of 90-180
        split (bool):                           If True, show two different plots for upward and downward directions
        log_scale (bool):                       If True, set a logarithmic scale
        polar_angles (numpy.ndarray or str):    Polar angles values (radian).
                                                If 'default', use smuthi.fields.default_polar_angles
        azimuthal_angles (numpy.ndarray or str):Azimuthal angle values (radian).
                                                If 'default', use smuthi.fields.default_azimuthal_angles
        angular_resolution (float):             If provided, angular arrays are generated with this angular resolution
                                                over the default angular range
    """

    infld = simulation.initial_field
    plst = simulation.particle_list
    lsys = simulation.layer_system
    far_field = ff.scattering_cross_section(initial_field=infld,
                                            particle_list=plst,
                                            layer_system=lsys,
                                            polar_angles=polar_angles,
                                            azimuthal_angles=azimuthal_angles,
                                            angular_resolution=angular_resolution)

    [d.setdefault('label','scattering_cross_section') for d in show_opts]

    show_far_field(far_field=far_field, save_plots=save_plots, save_opts=save_opts, show_plots=show_plots,
                   show_opts=show_opts, save_data=save_data, data_format=data_format, outputdir=outputdir,
                   flip_downward=flip_downward, split=split, log_scale=log_scale)


def show_far_field(far_field, show_plots=True, show_opts=[{'label':'far_field'}], save_plots=False, save_opts=None,
                   save_data=False, data_format='hdf5', outputdir='.',
                   flip_downward=True, split=True, log_scale=False):
    """Display and export the far field.

    Args:
        far_field (smuthi.field_expansion.FarField):    Far-field object to show and export
        show_plots (bool):      Display plots if True
        show_opts (dict list):  List of dictionaries containing options to be passed to pcolormesh for plotting.
                                If save_plots=True, a 1:1 correspondence between show_opts and save_opts dictionaries
                                is assumed. For simplicity, one can also provide a single show_opts entry that will
                                be applied to all save_opts.
                                The following keys are available (see matplotlib.pyplot.pcolormesh documentation):
                                'alpha'     (None)
                                'cmap'      ('inferno')
                                'norm'      (None), is set to matplotlib.colors.LogNorm() if log_scale is True
                                'vmin'      (None), applies only to 2D plots
                                'vmax'      (None), applies only to 2D plots
                                'shading'   ('nearest'), applies only to 2D plots. 'gouraud' is also available
                                'linewidth' (None), applies only to 1D plots
                                'linestyle' (None), applies only to 1D plots
                                'marker'    (None), applies only to 1D plots
                                An optional extra key called 'label' of type string is shown in the plot title
                                and appended to the associated file if save_plots is True
                                Finally, an optional 'figsize' key is available to set the width and height of
                                the figure window (see matplotlib.pyplot.figure documentation)
        save_plots (bool):      If True, plots are exported to file.
        save_opts (dict list):  List of dictionaries containing options to be passed to savefig.
                                A 1:1 correspondence between save_opts and show_opts dictionaries is assumed. For
                                simplicity, one can also provide a single save_opts entry that will be applied to
                                all show_opts.
                                The following keys are made available (see matplotlib.pyplot.savefig documentation):
                                'dpi'           (None)
                                'orientation'   (None)
                                'format'        ('png'), also available: eps, jpeg, jpg, pdf, ps, svg, tif, tiff ...
                                'transparent'   (False)
                                'bbox_inches'   ('tight')
                                'pad_inches'    (0.1)
        save_data (bool):       If True, raw data are exported to file
        data_format (str):      Output data format string, 'hdf5' and 'ascii' formats are available
        outputdir (str):        Path to the directory where files are to be saved
        flip_downward (bool):   If True, represent downward directions as 0-90 deg instead of 90-180
        split (bool):           If True, show two different plots for upward and downward directions
        log_scale (bool):       If True, set a logarithmic scale
    """

    if split and any(far_field.polar_angles < np.pi/2) and any(far_field.polar_angles > np.pi/2):
        for d in show_opts:
            d['label'] = d.get('label') + '_top'
        show_far_field(far_field.top(), show_plots, show_opts, save_plots, save_opts,
                       save_data, data_format, outputdir, True, False, log_scale)
        for d in show_opts:
            d['label'] = d.get('label').rstrip('_top') + '_bottom'
        show_far_field(far_field.bottom(), show_plots, show_opts, save_plots, save_opts,
                       save_data, data_format, outputdir, True, False, log_scale)
        return

    if (not os.path.exists(outputdir)) and (save_plots or save_data):
        os.makedirs(outputdir)

    if save_opts is None:
        save_opts = [{'format':'png'}]

    if save_data:
        pa = far_field.polar_angles
        aa = far_field.azimuthal_angles
        label = show_opts[0].get('label') # WARNING: filenames are based only on the label in show_opts[0]
        if data_format.lower() == 'hdf5':
            import h5py

            with h5py.File(outputdir + '/data.hdf5', 'a') as f:
                g = f.require_group('far_field')
                g.require_dataset('polar_angles', data=pa, shape=np.shape(pa), dtype=pa.dtype)
                g.require_dataset('azimuthal_angles', data=aa, shape=np.shape(aa), dtype=aa.dtype)
                g.create_dataset(label, data=far_field.signal, compression="gzip")
                g.create_dataset(label + '_polar', data=far_field.azimuthal_integral(), compression="gzip")
        elif data_format.lower() == 'ascii':
            np.savetxt(outputdir + '/' + label + '_TE.dat', far_field.signal[0, :, :],
                       header='Each line corresponds to a polar angle, each column corresponds to an azimuthal angle.')
            np.savetxt(outputdir + '/' + label + '_TM.dat', far_field.signal[1, :, :],
                       header='Each line corresponds to a polar angle, each column corresponds to an azimuthal angle.')
            np.savetxt(outputdir + '/' + label + '_polar_TE.dat', far_field.azimuthal_integral()[0, :],
                       header='Each line corresponds to a polar angle, each column corresponds to an azimuthal angle.')
            np.savetxt(outputdir + '/' + label + '_polar_TM.dat', far_field.azimuthal_integral()[1, :],
                       header='Each line corresponds to a polar angle, each column corresponds to an azimuthal angle.')
            np.savetxt(outputdir + '/polar_angles.dat', pa,
                       header='Polar angles of the far field in radians.')
            np.savetxt(outputdir + '/azimuthal_angles.dat', aa,
                       header='Azimuthal angles of the far field in radians.')
        else:
            raise ValueError('Currently, only hdf5 or ascii output data formats are available')

    alpha_grid = far_field.alpha_grid()
    beta_grid = far_field.beta_grid() * 180 / np.pi
    polar_array = far_field.polar_angles * 180 / np.pi
    if flip_downward and all(far_field.polar_angles >= np.pi / 2):
        beta_grid = 180 - beta_grid
        polar_array = 180 - polar_array

    if not show_plots:
        import matplotlib
        default_backend = matplotlib.get_backend()
        matplotlib.use('Agg') # a non-gui backend to avoid opening figure windows

    for show_opt, save_opt in zip(cycle(show_opts), save_opts) if len(show_opts) < len(save_opts) else zip(show_opts, cycle(save_opts)):
        # 2D polar plot of far field
        fig = plt.figure(figsize=show_opt.get('figsize',[6.4, 4.8]))
        ax = fig.add_subplot(111, polar=True)

        if log_scale: # in either case, this will be overridden if a 'norm' is also passed to show_opts
            color_norm = show_opt.get('norm', LogNorm(vmin=show_opt.get('vmin'), vmax=show_opt.get('vmax')))
        else:
            color_norm = show_opt.get('norm', Normalize(vmin=show_opt.get('vmin'), vmax=show_opt.get('vmax')))

        pcm = ax.pcolormesh(alpha_grid, beta_grid, (far_field.signal[0, :, :] + far_field.signal[1, :, :]),
                            alpha=show_opt.get('alpha'), norm=color_norm,
                            cmap=show_opt.get('cmap','inferno'), shading=show_opt.get('shading','nearest'))

        plt.colorbar(pcm, ax=ax)
        plt.title(show_opt.get('label').replace('_',' '))
        if save_plots:
            export_filename = outputdir + '/' + show_opt.get('label') + '.' + save_opt.get('format','png') # png unless specified
            plt.savefig(export_filename, dpi=save_opt.get('dpi'), orientation=save_opt.get('orientation','portrait'),
                        transparent=save_opt.get('transparent',False), bbox_inches=save_opt.get('bbox_inches','tight'),
                        pad_inches=save_opt.get('pad_inches',0.1))

        if not show_plots:
            plt.close(fig)

        # 1D polar plot of far field
        fig = plt.figure(figsize=show_opt.get('figsize',[6.4, 4.8]))

        plt.plot(polar_array, np.sum(far_field.azimuthal_integral(), axis=0) * np.pi / 180, alpha=show_opt.get('alpha'),
                 lw=show_opt.get('linewidth'), ls=show_opt.get('linestyle'), marker=show_opt.get('marker'))

        plt.xlabel('polar angle (degree)') # TODO: use radians instead?
        if far_field.signal_type == 'differential scattering cross section':
            plt.ylabel('d_CS/d_cos(beta)') # TODO: add units based on simulation.length_unit?
        elif far_field.signal_type == 'intensity':
            plt.ylabel('d_P/d_cos(beta)') # TODO: add units based on simulation.length_unit?

        if isinstance(color_norm, LogNorm):
            plt.yscale('log')
        elif isinstance(color_norm, SymLogNorm):
            linscale = color_norm._linscale_adj*(1.0 - 1/color_norm._base)
            plt.yscale('symlog', linthresh=color_norm.linthresh,
                                 base=color_norm._base, linscale=linscale)

        # the following line would apply the same vmin and vmax of 2D maps to 1D polar plots
        # plt.ylim([color_norm.vmin, color_norm.vmax])
        plt.grid(True)
        plt.title(show_opt.get('label').replace('_',' '))

        if save_plots:
            export_filename = outputdir + '/' + show_opt.get('label') + '_polar.' + save_opt.get('format','png')
            plt.savefig(export_filename, dpi=save_opt.get('dpi'), orientation=save_opt.get('orientation','portrait'),
                        transparent=save_opt.get('transparent',False), bbox_inches=save_opt.get('bbox_inches','tight'),
                        pad_inches=save_opt.get('pad_inches',0.1))
        if not show_plots:
            plt.close(fig)

    if not show_plots:
        matplotlib.use(default_backend) # now we can restore the original backend
    else:
        plt.show(block=False)

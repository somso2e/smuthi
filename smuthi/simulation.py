# -*- coding: utf-8 -*-
"""Provide class to manage a simulation."""

import smuthi.linearsystem.linear_system as lsys
import smuthi.fields as flds
import sys
import os
import datetime
import time
import shutil
import pickle
import numpy as np
import smuthi
import warnings
import smuthi.utility.logging as log


class Simulation:
    """Central class to manage a simulation.

    Args:
        layer_system (smuthi.layers.LayerSystem):               stratified medium
        particle_list (list):                                   list of smuthi.particles.Particle objects
        initial_field (smuthi.initial_field.InitialField):      initial field object
        k_parallel (numpy.ndarray or str):              in-plane wavenumber for Sommerfeld integrals and field
                                                        expansions. if 'default', use
                                                        smuthi.fields.default_Sommerfeld_k_parallel_array
        angular_resolution (float):                     If provided, angular arrays are generated with this angular
                                                        resolution over the default angular range
        neff_waypoints (list or ndarray):               Used to set default k_parallel arrays.
                                                        Corner points through which the contour runs
                                                        This quantity is dimensionless (effective
                                                        refractive index, will be multiplied by vacuum
                                                        wavenumber)
                                                        :ref:`MultipoleCutOffAnchor`
                                                        If not provided, reasonable waypoints are estimated.
        neff_imag (float):                              Used to set default k_parallel arrays.
                                                        Extent of the contour into the negative imaginary direction
                                                        (in terms of effective refractive index, n_eff=kappa/omega).
                                                        Only needed when no neff_waypoints are provided
        neff_max (float):                               Used to set default k_parallel arrays.
                                                        Truncation value of contour (in terms of effective refractive
                                                        index). Only needed when no neff_waypoints are
                                                        provided
        neff_max_offset (float):                        Used to set default k_parallel arrays.
                                                        Use the last estimated singularity location plus this value
                                                        (in terms of effective refractive index). Default=1
                                                        Only needed when no `neff_waypoints` are provided
                                                        and if no value for `neff_max` is specified.
        neff_resolution(float):                         Used to set default k_parallel arrays.
                                                        Resolution of contour, again in terms of effective refractive
                                                        index
        neff_minimal_branchpoint_distance (float):      Used to set default k_parallel arrays.
                                                        Minimal distance that contour points shall have from
                                                        branchpoint singularities (in terms of effective
                                                        refractive index). This is only relevant if not deflected
                                                        into imaginary. Default: One fifth of neff_resolution
        overwrite_default_contours (bool):              If true (default), the default contours are written even if
                                                        they have already been defined before
        solver_type (str):                      What solver type to use?
                                                Options: 'LU' for LU factorization, 'gmres' for GMRES iterative solver
        coupling_matrix_lookup_resolution (float or None): If type float, compute particle coupling by interpolation of
                                                           a lookup table with that spacial resolution. If None
                                                           (default), don't use a lookup table but compute the coupling
                                                           directly. This is more suitable for a small particle number.
        coupling_matrix_interpolator_kind (str):  Set to 'linear' (default) or 'cubic' interpolation of the lookup table.
        store_coupling_matrix (bool):           If True (default), the coupling matrix is stored. Otherwise it is
                                                recomputed on the fly during each iteration of the solver.
        length_unit (str):      what is the physical length unit? has no influence on the computations
        input_file (str):       path and filename of input file (for logging purposes)
        output_dir (str):       path to folder where to export data
        save_after_run(bool):   if true, the simulation object is exported to disc when over
        log_to_file(bool):      if true, the simulation log will be written to a log file
        log_to_terminal(bool):  if true, the simulation progress will be displayed in the terminal
        check_circumscribing_spheres(bool):  if true, check all particles for overlapping circumscribing spheres
                                             and print a warning if detected
        identical_particles (bool):          set this flag to true, if all particles have the same T-matrix (identical
                                             particles, located in the same background medium). Then, the T-matrix is
                                             computed only once for all particles.
         do_sanity_check (bool):             if true (default), check numerical input for some flaws. Warning: A passing
                                             sanity check does not guarantee correct numerical settings. For many
                                             particles, the sanity check might take some time and/or occupy large memory.
    """
    def __init__(self,
                 layer_system=None,
                 particle_list=None,
                 initial_field=None,
                 k_parallel='default',
                 angular_resolution=np.pi/360,
                 neff_waypoints=None,
                 neff_imag=1e-2,
                 neff_max=None,
                 neff_max_offset=1,
                 neff_resolution=1e-2,
                 neff_minimal_branchpoint_distance=None,
                 overwrite_default_contours=True,
                 solver_type='LU',
                 solver_tolerance=1e-4,
                 store_coupling_matrix=True,
                 coupling_matrix_lookup_resolution=None,
                 coupling_matrix_interpolator_kind='linear',
                 length_unit='length unit',
                 input_file=None,
                 output_dir='smuthi_output',
                 save_after_run=False,
                 log_to_file=False,
                 log_to_terminal=True,
                 check_circumscribing_spheres=True,
                 identical_particles=False,
                 do_sanity_check=True):

        # initialize attributes
        self.layer_system = layer_system
        self.particle_list = particle_list
        self.initial_field = initial_field
        self.k_parallel = k_parallel
        self.angular_resolution = angular_resolution
        self.neff_waypoints = neff_waypoints
        self.neff_imag = neff_imag
        self.neff_max = neff_max
        self.neff_max_offset = neff_max_offset
        self.neff_resolution = neff_resolution
        self.neff_minimal_branchpoint_distance = neff_minimal_branchpoint_distance
        self.overwrite_default_contours = overwrite_default_contours
        self.solver_type = solver_type
        self.solver_tolerance = solver_tolerance
        self.store_coupling_matrix = store_coupling_matrix
        self.coupling_matrix_lookup_resolution = coupling_matrix_lookup_resolution
        self.coupling_matrix_interpolator_kind = coupling_matrix_interpolator_kind
        self.post_processing = None
        self.length_unit = length_unit
        self.save_after_run = save_after_run
        self.check_circumscribing_spheres = check_circumscribing_spheres
        self.identical_particles = identical_particles
        self.do_sanity_check = do_sanity_check

        # output
        timestamp = '{:%Y%m%d%H%M%S}'.format(datetime.datetime.now())
        self.output_dir = output_dir + '/' + timestamp
        self.log_to_terminal = log_to_terminal
        self.log_to_file = log_to_file
        self.log_filename = self.output_dir + '/' + 'smuthi.log'
        self.set_logging()

        if input_file is not None and log_to_file:
            shutil.copyfile(input_file, self.output_dir + '/input.dat')

    def set_logging(self, log_to_terminal=None, log_to_file=None, log_filename=None):
        """Update logging behavior.

        Args:
            log_to_terminal (logical):  If true, print output to console.
            log_to_file (logical):      If true, print output to file
            log_filename (char):        If `log_to_file` is true, print output to a file with that name in the output
                                        directory. If the file already exists, it will be appended.
        """
        if log_to_terminal is not None:
            self.log_to_terminal = log_to_terminal
        if log_to_file is not None:
            self.log_to_file = log_to_file
        if log_filename is not None:
            self.log_filename = log_filename

        if not os.path.exists(self.output_dir) and self.log_to_file:
            os.makedirs(self.output_dir)
        sys.stdout = log.Logger(log_filename=self.log_filename,
                            log_to_file=self.log_to_file,
                            log_to_terminal=self.log_to_terminal)

    def __getstate__(self):
        """Return state values to be pickled."""
        return (self.layer_system, self.particle_list, self.initial_field, self.k_parallel, self.solver_type,
                self.solver_tolerance, self.store_coupling_matrix, self.coupling_matrix_lookup_resolution,
                self.coupling_matrix_interpolator_kind, self.post_processing, self.length_unit, self.save_after_run,
                flds.default_Sommerfeld_k_parallel_array, flds.default_initial_field_k_parallel_array,
                flds.default_polar_angles, flds.default_azimuthal_angles)

    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        (self.layer_system, self.particle_list, self.initial_field, self.k_parallel, self.solver_type,
         self.solver_tolerance, self.store_coupling_matrix, self.coupling_matrix_lookup_resolution,
         self.coupling_matrix_interpolator_kind, self.post_processing, self.length_unit, self.save_after_run,
         flds.default_Sommerfeld_k_parallel_array, flds.default_initial_field_k_parallel_array,
         flds.default_polar_angles, flds.default_azimuthal_angles) = state

    def print_simulation_header(self):
        smuthi.print_smuthi_header()
        sys.stdout.write("Starting simulation.\n")
        sys.stdout.flush()

    def save(self, filename=None):
        """Export simulation object to disc.

        Args:
            filename (str): path and file name where to store data
        """
        if filename is None:
            filename = self.output_dir + '/simulation.p'
        with open(filename, 'wb') as fn:
            pickle.dump(self, fn, -1)

    def sanity_check(self):
        """Check contour parameters for obvious problems"""
        maxrho = self.largest_lateral_distance()
        k = flds.angular_frequency(self.initial_field.vacuum_wavelength)
        thresh = 2 / (maxrho * k)

        # check neff_resol
        if self.neff_resolution is not None and self.k_parallel == "default":
            if self.neff_resolution > thresh:
                sys.stdout.write("Warning: Sanity check for neff_resolution failed!\n")
                sys.stdout.write("         neff_resolution = %.3e, max_rho = %.3e. k = %.3e.\n"%(self.neff_resolution,maxrho,k))
                sys.stdout.write("         It is recommended to select neff_resolution < 2 (k max_rho) = %.3e\n"%thresh)
                sys.stdout.flush()

        # check neff_imag
        if self.neff_imag is not None and self.k_parallel == "default" and self.neff_waypoints is None:
            if self.neff_imag > thresh:
                sys.stdout.write("Warning: Sanity check for neff_imag failed!\n")
                sys.stdout.write("         neff_imag = %.3e, max_rho = %.3e. k = %.3e.\n"%(self.neff_imag,maxrho,k))
                sys.stdout.write("         It is recommended to select neff_imag < 2 (k max_rho) = %.3e\n"%thresh)
                sys.stdout.flush()

    def largest_lateral_distance(self):
        """Compute the largest lateral distance between any two particles"""
        particle_x_array = np.array([particle.position[0] for particle in self.particle_list])
        particle_y_array = np.array([particle.position[1] for particle in self.particle_list])
        rho_squared = ((particle_x_array[:, None] - particle_x_array[None, :]) ** 2
                       + (particle_y_array[:, None] - particle_y_array[None, :]) ** 2)
        return np.sqrt(np.max(rho_squared))

    def initialize_linear_system(self):
        self.linear_system = lsys.LinearSystem(particle_list=self.particle_list,
                                               initial_field=self.initial_field,
                                               layer_system=self.layer_system,
                                               k_parallel=self.k_parallel,
                                               solver_type=self.solver_type,
                                               solver_tolerance=self.solver_tolerance,
                                               store_coupling_matrix=self.store_coupling_matrix,
                                               coupling_matrix_lookup_resolution=self.coupling_matrix_lookup_resolution,
                                               interpolator_kind=self.coupling_matrix_interpolator_kind,
                                               identical_particles=self.identical_particles)

    def circumscribing_spheres_disjoint(self):
        """Check if all circumscribing spheres are disjoint"""
        from scipy.spatial.distance import pdist, cdist, squareform
        pos = np.array([p.position for p in self.particle_list])
        csr = np.array([p.circumscribing_sphere_radius() for p in self.particle_list])
        try: # fully vectorized but memory expensive, finds all overlapping pairs
            distmatrix = squareform(pdist(pos))
            r1, r2 = np.meshgrid(csr, csr)
            overlap = np.triu(distmatrix < r1+r2, k=1)
            pidx = np.where(overlap)
            s = 's' if np.count_nonzero(overlap) > 1 else '' # pluralize message string
            msg = 'The circumscribing sphere of particle' + s
            for idx in pidx[0]:
                msg += " %i"%idx
            msg += ' overlaps with that of particle'+ s
            for idx in pidx[1]:
                msg += " %i"%idx
            msg += '.\n'
        except MemoryError: # less vectorized, more time consuming, stops at first overlap detected
            for i in range(len(self.particle_list)):
                dists = cdist(np.expand_dims(pos[i,:],0), pos[i+1:,])
                overlap = np.squeeze(dists < csr[i] + csr[i+1:])
                if overlap.any():
                    pidx = np.where(overlap)
                    msg = 'Found overlap between circumscribing sphere of particle %i and particle %i.\n'%(i,i+1+pidx[0][0])
                    break
        if overlap.any():
            sys.stdout.write(msg)
            sys.stdout.flush()
        return not overlap.any()

    def set_default_Sommerfeld_contour(self):
        """Set the default Sommerfeld k_parallel array"""

        flds.default_Sommerfeld_k_parallel_array = flds.reasonable_Sommerfeld_kpar_contour(
            vacuum_wavelength=self.initial_field.vacuum_wavelength,
            neff_waypoints=self.neff_waypoints,
            layer_refractive_indices=self.layer_system.refractive_indices,
            neff_imag=self.neff_imag,
            neff_max=self.neff_max,
            neff_max_offset=self.neff_max_offset,
            neff_resolution=self.neff_resolution,
            neff_minimal_branchpoint_distance=self.neff_minimal_branchpoint_distance)

    def set_default_angles(self):
        """Set the default polar and azimuthal angular arrays for pre-processing (i.e., initial field expansion)"""

        flds.default_azimuthal_angles, flds.default_polar_angles = \
        flds.angular_arrays(angular_resolution=self.angular_resolution)

    def set_default_initial_field_contour(self):
        """Set the default initial field k_parallel array"""

        if type(self.initial_field).__name__ == 'GaussianBeam':
            # in that case use only wavenumbers that propagate in the originating layer
            if self.initial_field.polar_angle <= np.pi / 2:
                neff_max = self.layer_system.refractive_indices[0].real
            else:
                neff_max = self.layer_system.refractive_indices[-1].real

            flds.default_initial_field_k_parallel_array = flds.reasonable_Sommerfeld_kpar_contour(
                vacuum_wavelength=self.initial_field.vacuum_wavelength,
                layer_refractive_indices=self.layer_system.refractive_indices,
                neff_imag=0,
                neff_max=neff_max,
                neff_resolution=self.neff_resolution,
                neff_minimal_branchpoint_distance=self.neff_minimal_branchpoint_distance)
        else:
            # case of dipoles etc ...
            # use a similar contour as for Sommerfeld integrals
            flds.default_initial_field_k_parallel_array = flds.reasonable_Sommerfeld_kpar_contour(
                vacuum_wavelength=self.initial_field.vacuum_wavelength,
                neff_waypoints=self.neff_waypoints,
                layer_refractive_indices=self.layer_system.refractive_indices,
                neff_imag=self.neff_imag,
                neff_max=self.neff_max,
                neff_max_offset=self.neff_max_offset,
                neff_resolution=self.neff_resolution,
                neff_minimal_branchpoint_distance=self.neff_minimal_branchpoint_distance)

    def set_default_contours(self):
        """Set the default initial field k_parallel array and the default Sommerfeld k_parallel array"""
        self.set_default_initial_field_contour()
        self.set_default_Sommerfeld_contour()

    def run(self):
        """Start the simulation."""
        self.print_simulation_header()

        # check for circumscribing sphere collisions.
        if self.check_circumscribing_spheres and len(self.particle_list) > 1 \
                and not self.circumscribing_spheres_disjoint():
            warnings.warn("Particles with circumscribing spheres detected.")

        # run sanity check
        if self.do_sanity_check and len(self.particle_list) > 1:
            self.sanity_check()

        # set default angular arrays
        if flds.default_azimuthal_angles is None or flds.default_polar_angles is None:
            self.set_default_angles()

        # check if default contours exists, otherwise set them
        if flds.default_Sommerfeld_k_parallel_array is None or self.overwrite_default_contours:
            self.set_default_Sommerfeld_contour()

        if flds.default_initial_field_k_parallel_array is None or self.overwrite_default_contours:
            self.set_default_initial_field_contour()

        # prepare and solve linear system
        start = time.time()
        self.initialize_linear_system()
        self.linear_system.prepare()
        end = time.time()
        preparation_time = end - start

        start = time.time()
        self.linear_system.solve()
        end = time.time()
        solution_time = end - start

        # post processing
        start = time.time()
        if self.post_processing:
            self.post_processing.run(self)
        end = time.time()
        postprocessing_time = end - start

        if self.save_after_run:
            self.save(self.output_dir + '/simulation.p')

        sys.stdout.write('\n')
        sys.stdout.flush()

        return preparation_time, solution_time, postprocessing_time

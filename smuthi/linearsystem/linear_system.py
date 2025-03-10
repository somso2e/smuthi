"""This package contains classes and functions to represent the system of
linear equations that needs to be solved in order to solve the scattering
problem, see section 3.7 of [Egel 2018 dissertation].

Symbolically, the linear system can be written like

.. math::
    (1 - TW)b = Ta,

where :math:`T` is the transition matrices of the particles, :math:`W` is the
particle coupling matrix, :math:`b` are the (unknown) coefficients of the
scattered field in terms of an outgoing spherical wave expansion and :math:`a`
are the coefficients of the initial field in terms of a regular spherical wave
expansion."""

import sys
import time
import inspect
import numpy as np
import scipy.special
import scipy.sparse.linalg
import scipy.linalg
import scipy.interpolate
from tqdm import tqdm
import warnings
import smuthi.fields as flds
import smuthi.linearsystem.tmatrix.t_matrix as tmt
import smuthi.fields.expansions as fldex
import smuthi.utility.cuda as cu
import smuthi.linearsystem.particlecoupling.direct_coupling as dircoup
import smuthi.linearsystem.particlecoupling.layer_mediated_coupling as laycoup
import smuthi.linearsystem.particlecoupling.prepare_lookup as look
import smuthi.linearsystem.linear_system_cuda as cusrc
from numba import config, set_num_threads
import smuthi.periodicboundaries as pb
import smuthi.periodicboundaries.coupling_helper as pbch


iter_num = 0


class LinearSystem:
    """Manage the assembly and solution of the linear system of equations.

    Args:
        particle_list (list):   List of smuthi.particles.Particle objects
        initial_field (smuthi.initial_field.InitialField):   Initial field object
        layer_system (smuthi.layers.LayerSystem):   Stratified medium
        k_parallel (numpy.ndarray or str): in-plane wavenumber.
                                           If 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
        solver_type (str):  What solver to use? Options: 'LU' for LU factorization, 'gmres' for GMRES iterative solver
        store_coupling_matrix (bool):   If True (default), the coupling matrix is stored. Otherwise it is recomputed on
                                        the fly during each iteration of the solver.
        coupling_matrix_lookup_resolution (float or None): If type float, compute particle coupling by interpolation of
                                                           a lookup table with that spacial resolution. A smaller number
                                                           implies higher accuracy and memory footprint.
                                                           If None (default), don't use a lookup table but compute the
                                                           coupling directly. This is more suitable for a small particle
                                                           number.
        interpolator_kind (str): interpolation order to be used, e.g. 'linear' or 'cubic'. This argument is ignored if
                                 coupling_matrix_lookup_resolution is None. In general, cubic interpolation is more
                                 accurate but a bit slower than linear.
        periodicity (tuple):                 tuple (a1, a2) specifying two 3-dimensional lattice vectors in Carthesian coordinates
                                             with a1, a2 (numpy.ndarrays)
        ewald_sum_separation_parameter (float):     Used to separate the real and reciprocal lattice sums to evaluate
                                                    particle coupling in periodic lattices.
        number_of_threads_periodic (int or str):    sets the number of threats used in a simulation with periodic particle arrangements
                                                    if 'default', all available CPU cores are used 
                                                    if negative, all but number_of_threads_periodic are used
        use_pvwf_coupling (bool):              If set to True, plane wave coupling is used to calculate
                                               the direct. Currently only possible in combination with direct solver strategy.
        pvwf_coupling_k_parallel (array):      k-parallel for PVWF coupling
    """

    def __init__(self,
                 particle_list,
                 initial_field,
                 layer_system,
                 k_parallel='default',
                 solver_type='LU',
                 solver_tolerance=1e-4,
                 store_coupling_matrix=True,
                 coupling_matrix_lookup_resolution=None,
                 interpolator_kind='cubic',
                 cuda_blocksize=None,
                 periodicity=None,
                 ewald_sum_separation_parameter='default',
                 number_of_threads_periodic='default',
                 use_pvwf_coupling=False,
                 pvwf_coupling_k_parallel=None):

        if cuda_blocksize is None:
            cuda_blocksize = cu.default_blocksize

        self.particle_list = particle_list
        self.initial_field = initial_field
        self.layer_system = layer_system
        self.k_parallel = k_parallel
        self.solver_type = solver_type
        self.solver_tolerance = solver_tolerance
        self.store_coupling_matrix = store_coupling_matrix
        self.coupling_matrix_lookup_resolution = coupling_matrix_lookup_resolution
        self.interpolator_kind = interpolator_kind
        self.cuda_blocksize = cuda_blocksize
        self.periodicity = periodicity
        self.ewald_sum_separation_parameter = ewald_sum_separation_parameter
        self.number_of_threads_periodic = number_of_threads_periodic
        self.use_pvwf_coupling = use_pvwf_coupling
        self.pvwf_coupling_k_parallel = pvwf_coupling_k_parallel

        if use_pvwf_coupling and (coupling_matrix_lookup_resolution is not None):
            warnings.warn(
                "The use_pvwf_coupling flag was defined in combination with a  "
                "lookup. This combination is not yet available. Ignoring the "
                "use_pvwf_coupling flag.")

        dummy_matrix = SystemMatrix(self.particle_list)
        sys.stdout.write('Number of unknowns: %i\n' % dummy_matrix.shape[0])

    def prepare(self):
        self.compute_initial_field_coefficients()
        self.compute_t_matrix()
        self.compute_coupling_matrix()
        self.master_matrix = MasterMatrix(t_matrix=self.t_matrix,
                                          coupling_matrix=self.coupling_matrix)

    def compute_initial_field_coefficients(self):
        """Evaluate initial field coefficients."""
        for particle in tqdm(self.particle_list,
                             desc='Initial field coefficients',
                             file=sys.stdout,
                             bar_format='{l_bar}{bar}| elapsed: {elapsed} remaining: {remaining}'):
            particle.initial_field = self.initial_field.spherical_wave_expansion(particle, self.layer_system)

    def compute_t_matrix(self):
        """Initialize T-matrix object."""
        for particle in tqdm(self.particle_list,
                             desc='T-matrices                ',
                             file=sys.stdout,
                             bar_format='{l_bar}{bar}| elapsed: {elapsed} remaining: {remaining}'):
            iS = self.layer_system.layer_number(particle.position[2])
            niS = self.layer_system.refractive_indices[iS]
            particle.t_matrix = particle.compute_t_matrix(self.initial_field.vacuum_wavelength, niS)
            if not particle.euler_angles == [0, 0, 0]:
                particle.t_matrix = tmt.rotate_t_matrix(particle.t_matrix, particle.l_max, particle.m_max,
                                                        particle.euler_angles, wdsympy=False)
        self.t_matrix = TMatrix(particle_list=self.particle_list)

    def compute_coupling_matrix(self):
        """Initialize coupling matrix object."""

        # LU factorization only possible with explicit matrix.
        if self.solver_type == 'LU':
            self.store_coupling_matrix = True

        if self.coupling_matrix_lookup_resolution is not None:
            z_list = [particle.position[2] for particle in self.particle_list]
            is_list = [self.layer_system.layer_number(z) for z in z_list]
            if not is_list.count(is_list[0]) == len(is_list):  # all particles in same layer?
                warnings.warn("Particles are not all in same layer. "
                              "Fall back to direct coupling matrix computation (no lookup).")
                self.coupling_matrix_lookup_resolution = None
            if self.periodicity is not None:
                warnings.warn("Periodic particle arrangement defined. "
                              "Disabling lookup.")
                self.coupling_matrix_lookup_resolution = None
            if self.store_coupling_matrix:
                warnings.warn("Explicit matrix computation using lookup currently not implemented. "
                              "Disabling lookup.")
                self.coupling_matrix_lookup_resolution = None
            else:  # use lookup
                if not self.interpolator_kind in ('linear', 'cubic'):
                    warnings.warn(self.interpolator_kind + ' interpolation not implemented. '
                                                           'Use "linear" instead')
                    self.interpolator_kind = 'linear'

                z_list = [particle.position[2] for particle in self.particle_list]
                if z_list.count(z_list[0]) == len(z_list):  # all particles at same height: use radial lookup
                    if cu.use_gpu:
                        sys.stdout.write('Coupling matrix computation by ' + self.interpolator_kind
                                         + ' interpolation of radial lookup on GPU.\n')
                        sys.stdout.flush()
                        self.coupling_matrix = CouplingMatrixRadialLookupCUDA(
                            vacuum_wavelength=self.initial_field.vacuum_wavelength,
                            particle_list=self.particle_list,
                            layer_system=self.layer_system,
                            k_parallel=self.k_parallel,
                            resolution=self.coupling_matrix_lookup_resolution,
                            cuda_blocksize=self.cuda_blocksize,
                            interpolator_kind=self.interpolator_kind)
                    else:
                        sys.stdout.write('Coupling matrix computation by ' + self.interpolator_kind
                                         + ' interpolation of radial lookup on CPU.\n')
                        sys.stdout.flush()
                        self.coupling_matrix = CouplingMatrixRadialLookupCPU(
                            vacuum_wavelength=self.initial_field.vacuum_wavelength,
                            particle_list=self.particle_list,
                            layer_system=self.layer_system,
                            k_parallel=self.k_parallel,
                            resolution=self.coupling_matrix_lookup_resolution,
                            interpolator_kind=self.interpolator_kind)
                else:  # not all particles at same height: use volume lookup
                    if cu.use_gpu:
                        sys.stdout.write('Coupling matrix computation by ' + self.interpolator_kind
                                         + ' interpolation of 3D lookup on GPU.\n')
                        sys.stdout.flush()
                        self.coupling_matrix = CouplingMatrixVolumeLookupCUDA(
                            vacuum_wavelength=self.initial_field.vacuum_wavelength,
                            particle_list=self.particle_list,
                            layer_system=self.layer_system,
                            k_parallel=self.k_parallel,
                            resolution=self.coupling_matrix_lookup_resolution,
                            interpolator_kind=self.interpolator_kind)
                    else:
                        sys.stdout.write('Coupling matrix computation by ' + self.interpolator_kind
                                         + ' interpolation of 3D lookup on CPU.\n')

                        sys.stdout.flush()
                        self.coupling_matrix = CouplingMatrixVolumeLookupCPU(
                            vacuum_wavelength=self.initial_field.vacuum_wavelength,
                            particle_list=self.particle_list,
                            layer_system=self.layer_system,
                            k_parallel=self.k_parallel,
                            resolution=self.coupling_matrix_lookup_resolution,
                            interpolator_kind=self.interpolator_kind)

        if self.coupling_matrix_lookup_resolution is None:
            if not self.store_coupling_matrix:
                warnings.warn("With lookup disabled, coupling matrix needs to be stored.")
                self.store_coupling_matrix = True
                
            if self.periodicity:
                sys.stdout.write('Explicit, periodic coupling matrix computation on CPU.\n')
                sys.stdout.flush()
                self.coupling_matrix = CouplingMatrixPeriodicGridNumba(
                    initial_field=self.initial_field,
                    particle_list=self.particle_list,
                    layer_system=self.layer_system,
                    periodicity=self.periodicity,
                    ewald_sum_separation_parameter=self.ewald_sum_separation_parameter,
                    num_threads=self.number_of_threads_periodic)
            
            else:
                sys.stdout.write('Explicit coupling matrix computation on CPU.\n')
                sys.stdout.flush()
                self.coupling_matrix = CouplingMatrixExplicit(
                    vacuum_wavelength=self.initial_field.vacuum_wavelength,
                    particle_list=self.particle_list,
                    layer_system=self.layer_system,
                    k_parallel=self.k_parallel,
                    use_pvwf_coupling=self.use_pvwf_coupling,
                    pvwf_coupling_k_parallel=self.pvwf_coupling_k_parallel)

    def solve(self):
        """Compute scattered field coefficients and store them
        in the particles' spherical wave expansion objects."""
        if len(self.particle_list) > 0:
            if self.solver_type.lower() == 'lu':
                sys.stdout.write('Solve (LU decomposition)  : ...')
                if not hasattr(self.master_matrix.linear_operator, 'A'):
                    raise ValueError('LU factorization only possible '
                                     'with the option "store coupling matrix".')
                if not hasattr(self.master_matrix, 'LU_piv'):
                    lu, piv = scipy.linalg.lu_factor(self.master_matrix.linear_operator.A,
                                                     overwrite_a=False)
                    self.master_matrix.LU_piv = (lu, piv)
                b = scipy.linalg.lu_solve(self.master_matrix.LU_piv,
                                          self.t_matrix.right_hand_side())
                sys.stdout.write(' done\n')
                sys.stdout.flush()
            elif self.solver_type.lower() == 'gmres':
                rhs = self.t_matrix.right_hand_side()
                start_time = time.time()
                iter_num = 0

                def status_msg(rk):
                    nonlocal iter_num
                    iter_msg = ('Solve (GMRES)             : Iter ' + str(iter_num)
                                + ' | Rel. residual: '
                                + "{:.2e}".format(np.linalg.norm(rk))
                                + ' | elapsed: ' + str(int(time.time() - start_time)) + 's')
                    sys.stdout.write('\r' + iter_msg)
                    iter_num += 1

                b, info = scipy.sparse.linalg.gmres(self.master_matrix.linear_operator, rhs, rhs,
                                                    tol=self.solver_tolerance, callback=status_msg)
            #                sys.stdout.write('\n')
            elif self.solver_type.lower() == 'lgmres':
                rhs = self.t_matrix.right_hand_side()
                start_time = time.time()
                iter_num = 0

                def status_msg(xk):
                    frame = inspect.currentframe().f_back
                    # for some reason, lgmres passes xk to callback instead of rk ...
                    # then, we could inspect r_norm from the currentframe, but it's not there ... (?)
                    # therefore we have to re-calculate it again, which is a total waste of time ...
                    rk = frame.f_locals['nrm2'](frame.f_locals['r_outer'])/frame.f_locals['b_norm']
                    nonlocal iter_num
                    iter_msg = ('Solve (LGMRES)            : Iter ' + str(iter_num)
                                + ' | Rel. residual: '
                                + "{:.2e}".format(np.linalg.norm(rk))
                                + ' | elapsed: ' + str(int(time.time() - start_time)) + 's')
                    sys.stdout.write('\r' + iter_msg)
                    iter_num += 1

                b, info = scipy.sparse.linalg.lgmres(self.master_matrix.linear_operator, rhs, rhs,
                                                    tol=self.solver_tolerance, callback=status_msg)
            #                sys.stdout.write('\n')
            elif self.solver_type.lower() == 'gcrotmk':
                rhs = self.t_matrix.right_hand_side()
                start_time = time.time()
                iter_num = 0

                def status_msg(xk):
                    frame = inspect.currentframe().f_back
                    # also gcrotmk passes xk to callback instead of rk ...
                    # let's calculate residuals on the fly also in this case ...
                    rk = frame.f_locals['nrm2'](frame.f_locals['r'])/frame.f_locals['b_norm']
                    nonlocal iter_num
                    iter_msg = ('Solve (GCROTMK)           : Iter ' + str(iter_num)
                                + ' | Rel. residual: '
                                + "{:.2e}".format(np.linalg.norm(rk))
                                + ' | elapsed: ' + str(int(time.time() - start_time)) + 's')
                    sys.stdout.write('\r' + iter_msg)
                    iter_num += 1

                b, info = scipy.sparse.linalg.gcrotmk(self.master_matrix.linear_operator, rhs, rhs,
                                                    tol=self.solver_tolerance, callback=status_msg)
            #                sys.stdout.write('\n')
            else:
                raise ValueError('This solver type is currently not implemented.')

        for iS, particle in enumerate(self.particle_list):
            i_iS = self.layer_system.layer_number(particle.position[2])
            n_iS = self.layer_system.refractive_indices[i_iS]
            k = flds.angular_frequency(self.initial_field.vacuum_wavelength) * n_iS
            loz, upz = self.layer_system.lower_zlimit(i_iS), self.layer_system.upper_zlimit(i_iS)
            particle.scattered_field = fldex.SphericalWaveExpansion(k=k, l_max=particle.l_max, m_max=particle.m_max,
                                                                    kind='outgoing', reference_point=particle.position,
                                                                    lower_z=loz, upper_z=upz)
            particle.scattered_field.coefficients = b[self.master_matrix.index_block(iS)]
            particle.scattered_field.validity_conditions.append(particle.is_outside)
            self.initial_field.validity_conditions.append(particle.is_outside)


class SystemMatrix:
    r"""A system matrix is an abstract linear operator that operates on a
    system coefficient vector, i.e. a vector :math:`c = c_{\tau,l,m}^i`, where
    :math:`(\tau, l, m)` are the multipole indices and :math:`i` indicates the
    particle number. In other words, if we have a spherical wave expansion for
    each particle, and write all the expansion coefficients of these expansions
    into one (long) array, what we get is a system vector.
    """

    def __init__(self, particle_list):
        self.particle_list = particle_list
        blocksizes = [flds.blocksize(particle.l_max, particle.m_max) for particle in self.particle_list]
        self.shape = (sum(blocksizes), sum(blocksizes))

    def index_block(self, i):
        """
        Args:
            i (int): number of particle

        Returns:
            indices that correspond to the coefficients for that particle
        """
        blocksizes = [flds.blocksize(particle.l_max, particle.m_max) for particle in self.particle_list[:(i + 1)]]
        return range(sum(blocksizes[:i]), sum(blocksizes))

    def index(self, i, tau, l, m):
        r"""
        Args:
            i (int):    particle number
            tau (int):    spherical polarization index
            l (int):    multipole degree
            m (int):    multipole order

        Returns:
            Position in a system vector that corresponds to the :math:`(\tau, l, m)` coefficient of the i-th particle.
        """
        blocksizes = [flds.blocksize(particle.l_max, particle.m_max) for particle in self.particle_list[:i]]
        return sum(blocksizes) + flds.multi_to_single_index(tau, l, m, self.particle_list[i].l_max,
                                                             self.particle_list[i].m_max)


class MasterMatrix(SystemMatrix):
    r"""Represent the master matrix :math:`M = 1 - TW` as a linear operator.

    Args:
        t_matrix (SystemMatrix):          System T-matrix
        coupling_matrix (SystemMatrix):   System coupling matrix
    """

    def __init__(self, t_matrix, coupling_matrix):
        SystemMatrix.__init__(self, t_matrix.particle_list)
        if type(coupling_matrix.linear_operator).__name__ == 'MatrixLinearOperator':
            M = (np.eye(self.shape[0], dtype=complex)
                 - t_matrix.linear_operator.matmat(coupling_matrix.linear_operator.A))
            self.linear_operator = scipy.sparse.linalg.aslinearoperator(M)
        else:
            def apply_master_matrix(vector):
                return vector - t_matrix.linear_operator.dot(coupling_matrix.linear_operator.matvec(vector))

            self.linear_operator = scipy.sparse.linalg.LinearOperator(shape=self.shape, matvec=apply_master_matrix,
                                                                      dtype=complex)


class TMatrix(SystemMatrix):
    """Collect the particle T-matrices in a global lienear operator.

    Args:
        particle_list (list):   List of smuthi.particles.Particle objects containing a t_matrix attribute.
    """

    def __init__(self, particle_list):
        SystemMatrix.__init__(self, particle_list)

        def apply_t_matrix(vector):
            tv = np.zeros(vector.shape, dtype=complex)
            for i_s, particle in enumerate(particle_list):
                tv[self.index_block(i_s)] = particle.t_matrix.dot(vector[self.index_block(i_s)])
            return tv

        self.linear_operator = scipy.sparse.linalg.LinearOperator(shape=self.shape, matvec=apply_t_matrix,
                                                                  matmat=apply_t_matrix, dtype=complex)

    def right_hand_side(self):
        r"""The right hand side of the linear system is given by :math:`\sum_{\tau l m} T^i_{\tau l m} a^i_{\tau l m }`

        Returns:
            right hand side as a complex numpy.ndarray
        """
        tai = np.zeros(self.shape[0], dtype=complex)
        for i_s, particle in enumerate(self.particle_list):
            tai[self.index_block(i_s)] = particle.t_matrix.dot(particle.initial_field.coefficients)
        return tai


class CouplingMatrixExplicit(SystemMatrix):
    """Class for an explicit representation of the coupling matrix. Recommended for small particle numbers.

    Args:
        vacuum_wavelength (float):  Vacuum wavelength in length units
        particle_list (list):   List of smuthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem):   Stratified medium
        k_parallell (numpy.ndarray or str): In-plane wavenumber.
                                            If 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
    """

    def __init__(self, vacuum_wavelength, particle_list, layer_system, k_parallel='default', use_pvwf_coupling=False, pvwf_coupling_k_parallel=None):

        SystemMatrix.__init__(self, particle_list)
        coup_mat = np.zeros(self.shape, dtype=complex)
        sys.stdout.write('Coupling matrix memory footprint: ' + look.size_format(coup_mat.nbytes) + '\n')
        sys.stdout.flush()
        
        if layer_system.is_degenerate():
            sys.stdout.write('Omitting calculations of layer mediated coupling since no difference in layers was detected. \n')
            sys.stdout.flush()
        
        for s1, particle1 in enumerate(tqdm(particle_list, desc='Particle coupling matrix  ', file=sys.stdout,
                                            bar_format='{l_bar}{bar}| elapsed: {elapsed} ' 'remaining: {remaining}')):
            idx1 = np.array(self.index_block(s1))[:, None]
            for s2, particle2 in enumerate(particle_list):
                idx2 = self.index_block(s2)

                # direct contribution
                if use_pvwf_coupling:
                    coup_mat[idx1, idx2] = dircoup.direct_coupling_block_pvwf_mediated(vacuum_wavelength,
                                                                                       particle1,
                                                                                       particle2,
                                                                                       layer_system,
                                                                                       pvwf_coupling_k_parallel)
                else:
                    coup_mat[idx1, idx2] = dircoup.direct_coupling_block(vacuum_wavelength,
                                                                         particle1,
                                                                         particle2,
                                                                         layer_system)

                # layer mediated contribution
                if not layer_system.is_degenerate():
                    coup_mat[idx1, idx2] += laycoup.layer_mediated_coupling_block(vacuum_wavelength, particle1, particle2,
                                                                                  layer_system, k_parallel)
        self.linear_operator = scipy.sparse.linalg.aslinearoperator(coup_mat)


class CouplingMatrixPeriodicGridNumba(SystemMatrix):
    """ Class for an explicit representation of the coupling matrix of periodic particle arrangements.
        Computation supports Numba.
    Args:
        initial_field (smuthi.initial_field.PlaneWave):     initial plane wave object
        particle_list (list):                               list of smuthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem):           stratified medium
        periodicity (tuple):                                (a1, a2) lattice vector 1 and 2 in carthesian coordinates
        ewald_sum_separation_parameter (float):             Ewald sum separation parameter
        num_threads (int or str):                           if 'default' all available CPU cores are used
                                                            if negative, all but num_threads are used
    """
    
    def __init__(self, initial_field, particle_list, layer_system, periodicity,
                 ewald_sum_separation_parameter, num_threads='default'):
        
        max_num_threads = config.NUMBA_DEFAULT_NUM_THREADS
        if type(num_threads).__name__ == 'int':
            if 0 < num_threads <= max_num_threads:
                set_num_threads(num_threads)
                sys.stdout.write('Parallel Numba execution limited to %d of %d available CPU cores. \n' % (num_threads, max_num_threads))
            elif num_threads <= 0:
                set_num_threads(max_num_threads + num_threads)
                sys.stdout.write('Parallel Numba execution limited to %d of %d available CPU cores. \n' % (max_num_threads + num_threads, max_num_threads))
            else:
                sys.stdout.write('Parallel Numba execution uses all %d available CPU cores. \n' % (max_num_threads))
        else:
            sys.stdout.write('Parallel Numba execution uses all %d available CPU cores. \n' % (max_num_threads))
        sys.stdout.flush()
        
        a1, a2 = periodicity
        
        if ewald_sum_separation_parameter == 'default':
            eta = pb.default_Ewald_sum_separation
        else:
            eta = ewald_sum_separation_parameter
        # prepare a5b5 lookup table
        lmax_global = np.max([particle.l_max for particle in particle_list])
        mmax_global = np.max([particle.m_max for particle in particle_list])
        a5b5_mat = pbch.a5b5_lookup(lmax_global, lmax_global, mmax_global, mmax_global)
        
        i_sca = layer_system.layer_number(particle_list[0].position[2]) # all particles are within one layer
        if initial_field.polar_angle < np.pi:
            pwe_exc = initial_field.plane_wave_expansion(layer_system, i_sca)[0]
        else:
            pwe_exc = initial_field.plane_wave_expansion(layer_system, i_sca)[1]
        k0t = np.array([pwe_exc.k_parallel[0] * np.cos(pwe_exc.azimuthal_angles)[0],
                        pwe_exc.k_parallel[0] * np.sin(pwe_exc.azimuthal_angles)[0]]) 
        
        num_particles = len(particle_list)
        positions = np.zeros([num_particles, 3], np.float64)
        radii = np.zeros([num_particles], np.float64)
        lmax_array = np.zeros([num_particles], np.int32)
        mmax_array = np.zeros([num_particles], np.int32)
        for idx, particle in enumerate(particle_list):
            positions[idx] = particle.position
            lmax_array[idx] = particle.l_max
            mmax_array[idx] = particle.m_max
            if type(particle).__name__ == 'Sphere':
                radii[idx] = particle.radius
            elif type(particle).__name__ == 'FiniteCylinder':
                radii[idx] = np.sqrt(particle.cylinder_radius ** 2 + (particle.cylinder_height / 2) ** 2)

        import smuthi.periodicboundaries.particle_coupling as pbcoup
        coup_mat = pbcoup.periodic_coupling_matrix(initial_field.vacuum_wavelength, k0t,
                            pwe_exc.azimuthal_angles[0], np.array(layer_system.thicknesses, np.float64),
                            np.array(layer_system.refractive_indices, np.complex128),
                            i_sca, positions, radii, lmax_array, mmax_array, a1, a2, eta, a5b5_mat, mmax_global)
        sys.stdout.write('Coupling matrix computation finished. \n')
        sys.stdout.flush()
        
        self.linear_operator = scipy.sparse.linalg.aslinearoperator(coup_mat)


class CouplingMatrixVolumeLookup(SystemMatrix):
    """Base class for 3D lookup based coupling matrix either on CPU or on GPU (CUDA).

    Args:
        vacuum_wavelength (float): vacuum wavelength in length units
        particle_list (list): list of sumthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem): stratified medium
        k_parallel (numpy.ndarray or str): in-plane wavenumber.
                                           If 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
        resolution (float or None): spatial resolution of the lookup in the radial direction
    """

    def __init__(self, vacuum_wavelength, particle_list, layer_system, k_parallel='default', resolution=None):
        z_list = [particle.position[2] for particle in particle_list]
        is_list = [layer_system.layer_number(z) for z in z_list]
        assert is_list.count(is_list[0]) == len(is_list)  # all particles in same layer?

        SystemMatrix.__init__(self, particle_list)

        self.l_max = max([particle.l_max for particle in particle_list])
        self.m_max = max([particle.m_max for particle in particle_list])
        self.blocksize = flds.blocksize(self.l_max, self.m_max)
        self.resolution = resolution
        lkup = look.volumetric_coupling_lookup_table(vacuum_wavelength=vacuum_wavelength, particle_list=particle_list,
                                                     layer_system=layer_system, k_parallel=k_parallel,
                                                     resolution=resolution)
        self.lookup_table_plus, self.lookup_table_minus = lkup[0], lkup[1]
        self.rho_array, self.sum_z_array, self.diff_z_array = lkup[2], lkup[3], lkup[4]


class CouplingMatrixVolumeLookupCPU(CouplingMatrixVolumeLookup):
    """Class for 3D lookup based coupling matrix running on CPU. This is used when no suitable GPU device is detected
    or when PyCuda is not installed.

    Args:
        vacuum_wavelength (float): vacuum wavelength in length units
        particle_list (list): list of sumthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem): stratified medium
        k_parallel (numpy.ndarray or str): in-plane wavenumber.
                                           If 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
        resolution (float or None): spatial resolution of the lookup in the radial direction
        interpolator_kind (str): 'linear' or 'cubic' interpolation
    """

    def __init__(self, vacuum_wavelength, particle_list, layer_system, k_parallel='default', resolution=None,
                 interpolator_kind='cubic'):

        if interpolator_kind == 'cubic':
            interpolation_order = 3
        else:
            interpolation_order = 1

        CouplingMatrixVolumeLookup.__init__(self, vacuum_wavelength, particle_list, layer_system, k_parallel,
                                            resolution)

        x_array = np.array([particle.position[0] for particle in particle_list])
        y_array = np.array([particle.position[1] for particle in particle_list])
        z_array = np.array([particle.position[2] for particle in particle_list])

        self.particle_rho_array = np.sqrt((x_array[:, None] - x_array[None, :]) ** 2
                                          + (y_array[:, None] - y_array[None, :]) ** 2)
        self.particle_phi_array = np.arctan2(y_array[:, None] - y_array[None, :], x_array[:, None] - x_array[None, :])
        self.particle_sz_array = z_array[:, None] + z_array[None, :]
        self.particle_dz_array = z_array[:, None] - z_array[None, :]

        # contains for each n all positions in the large system arrays that correspond to n:
        self.system_vector_index_list = [[] for i in range(self.blocksize)]

        # same size as system_vector_index_list, contains the according particle numbers:
        self.particle_number_list = [[] for i in range(self.blocksize)]
        self.m_list = [None for i in range(self.blocksize)]
        for i, particle in enumerate(particle_list):
            for m in range(-particle.m_max, particle.m_max + 1):
                for l in range(max(1, abs(m)), particle.l_max + 1):
                    for tau in range(2):
                        n_lookup = flds.multi_to_single_index(tau=tau, l=l, m=m, l_max=self.l_max, m_max=self.m_max)
                        self.system_vector_index_list[n_lookup].append(self.index(i, tau, l, m))
                        self.particle_number_list[n_lookup].append(i)
                        self.m_list[n_lookup] = m
        for n in range(self.blocksize):
            self.system_vector_index_list[n] = np.array(self.system_vector_index_list[n])
            self.particle_number_list[n] = np.array(self.particle_number_list[n])

        self.lookup_plus_real = [[None for i in range(self.blocksize)] for i2 in range(self.blocksize)]
        self.lookup_plus_imag = [[None for i in range(self.blocksize)] for i2 in range(self.blocksize)]
        self.lookup_minus_real = [[None for i in range(self.blocksize)] for i2 in range(self.blocksize)]
        self.lookup_minus_imag = [[None for i in range(self.blocksize)] for i2 in range(self.blocksize)]
        for n1 in range(self.blocksize):
            for n2 in range(self.blocksize):
                self.lookup_plus_real[n1][n2] = scipy.interpolate.RectBivariateSpline(
                    x=self.rho_array, y=self.sum_z_array, z=self.lookup_table_plus[:, :, n1, n2].real,
                    kx=interpolation_order, ky=interpolation_order)
                self.lookup_plus_imag[n1][n2] = scipy.interpolate.RectBivariateSpline(
                    x=self.rho_array, y=self.sum_z_array, z=self.lookup_table_plus[:, :, n1, n2].imag,
                    kx=interpolation_order, ky=interpolation_order)
                self.lookup_minus_real[n1][n2] = scipy.interpolate.RectBivariateSpline(
                    x=self.rho_array, y=self.diff_z_array, z=self.lookup_table_minus[:, :, n1, n2].real,
                    kx=interpolation_order, ky=interpolation_order)
                self.lookup_minus_imag[n1][n2] = scipy.interpolate.RectBivariateSpline(
                    x=self.rho_array, y=self.diff_z_array, z=self.lookup_table_minus[:, :, n1, n2].imag,
                    kx=interpolation_order, ky=interpolation_order)

        def matvec(in_vec):
            out_vec = np.zeros(shape=in_vec.shape, dtype=complex)
            for n1 in range(self.blocksize):
                i1 = self.particle_number_list[n1]
                idx1 = self.system_vector_index_list[n1]
                m1 = self.m_list[n1]
                for n2 in range(self.blocksize):
                    i2 = self.particle_number_list[n2]
                    idx2 = self.system_vector_index_list[n2]
                    m2 = self.m_list[n2]
                    rho = self.particle_rho_array[i1[:, None], i2[None, :]]
                    phi = self.particle_phi_array[i1[:, None], i2[None, :]]
                    sz = self.particle_sz_array[i1[:, None], i2[None, :]]
                    dz = self.particle_dz_array[i1[:, None], i2[None, :]]
                    Mpl = self.lookup_plus_real[n1][n2].ev(rho, sz) + 1j * self.lookup_plus_imag[n1][n2].ev(rho, sz)
                    Mmn = self.lookup_minus_real[n1][n2].ev(rho, dz) + 1j * self.lookup_minus_imag[n1][n2].ev(rho, dz)
                    M = (Mpl + Mmn) * np.exp(1j * (m2 - m1) * phi)
                    out_vec[idx1] += M.dot(in_vec[idx2])
            return out_vec

        self.linear_operator = scipy.sparse.linalg.LinearOperator(shape=self.shape, matvec=matvec, dtype=complex)


class CouplingMatrixVolumeLookupCUDA(CouplingMatrixVolumeLookup):
    """Class for 3D lookup based coupling matrix running on GPU.

    Args:
        vacuum_wavelength (float): vacuum wavelength in length units
        particle_list (list): list of sumthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem): stratified medium
        k_parallel (numpy.ndarray or str): in-plane wavenumber.
                                           If 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
        resolution (float or None): spatial resolution of the lookup in the radial direction
        cuda_blocksize (int): threads per block for cuda call
        interpolator_kind (str): 'linear' (default) or 'cubic' interpolation
    """

    def __init__(self, vacuum_wavelength, particle_list, layer_system, k_parallel='default', resolution=None,
                 cuda_blocksize=None, interpolator_kind='linear'):

        if cuda_blocksize is None:
            cuda_blocksize = cu.default_blocksize

        CouplingMatrixVolumeLookup.__init__(self, vacuum_wavelength, particle_list, layer_system, k_parallel,
                                            resolution)

        sys.stdout.write('Prepare CUDA kernel and device lookup data ... ')
        sys.stdout.flush()
        start_time = time.time()

        if interpolator_kind == 'linear':
            coupling_source = cusrc.linear_volume_lookup_source % (self.blocksize, self.shape[0], len(self.sum_z_array),
                                                                min(self.rho_array), min(self.sum_z_array),
                                                                min(self.diff_z_array), self.resolution)
        elif interpolator_kind == 'cubic':
            coupling_source = cusrc.cubic_volume_lookup_source % (self.blocksize, self.shape[0], len(self.sum_z_array),
                                                               min(self.rho_array), min(self.sum_z_array),
                                                               min(self.diff_z_array), self.resolution)

        coupling_function = cu.SourceModule(coupling_source).get_function("coupling_kernel")

        n_lookup_array = np.zeros(self.shape[0], dtype=np.uint32)
        m_particle_array = np.zeros(self.shape[0], dtype=np.float32)
        x_array = np.zeros(self.shape[0], dtype=np.float32)
        y_array = np.zeros(self.shape[0], dtype=np.float32)
        z_array = np.zeros(self.shape[0], dtype=np.float32)

        i_particle = 0
        for i, particle in enumerate(particle_list):
            for m in range(-particle.m_max, particle.m_max + 1):
                for l in range(max(1, abs(m)), particle.l_max + 1):
                    for tau in range(2):
                        i_taulm = flds.multi_to_single_index(tau, l, m, particle.l_max, particle.m_max)
                        idx = i_particle + i_taulm

                        n_lookup_array[idx] = flds.multi_to_single_index(tau, l, m, self.l_max, self.m_max)
                        m_particle_array[idx] = m

                        # scale the x and y position to the lookup resolution:
                        x_array[idx] = particle.position[0]
                        y_array[idx] = particle.position[1]
                        z_array[idx] = particle.position[2]

            i_particle += flds.blocksize(particle.l_max, particle.m_max)

        re_lookup_pl = self.lookup_table_plus.real.astype(dtype=np.float32)
        im_lookup_pl = self.lookup_table_plus.imag.astype(dtype=np.float32)
        re_lookup_mn = self.lookup_table_minus.real.astype(dtype=np.float32)
        im_lookup_mn = self.lookup_table_minus.imag.astype(dtype=np.float32)

        # transfer data to gpu
        n_lookup_array_d = cu.gpuarray.to_gpu(n_lookup_array)
        m_particle_array_d = cu.gpuarray.to_gpu(m_particle_array)
        x_array_d = cu.gpuarray.to_gpu(x_array)
        y_array_d = cu.gpuarray.to_gpu(y_array)
        z_array_d = cu.gpuarray.to_gpu(z_array)
        re_lookup_pl_d = cu.gpuarray.to_gpu(re_lookup_pl)
        im_lookup_pl_d = cu.gpuarray.to_gpu(im_lookup_pl)
        re_lookup_mn_d = cu.gpuarray.to_gpu(re_lookup_mn)
        im_lookup_mn_d = cu.gpuarray.to_gpu(im_lookup_mn)

        sys.stdout.write('done | elapsed: ' + str(int(time.time() - start_time)) + 's\n')
        sys.stdout.flush()

        cuda_gridsize = (self.shape[0] + cuda_blocksize - 1) // cuda_blocksize

        def matvec(in_vec):
            re_in_vec_d = cu.gpuarray.to_gpu(np.float32(in_vec.real))
            im_in_vec_d = cu.gpuarray.to_gpu(np.float32(in_vec.imag))
            re_result_d = cu.gpuarray.zeros(in_vec.shape, dtype=np.float32)
            im_result_d = cu.gpuarray.zeros(in_vec.shape, dtype=np.float32)
            coupling_function(n_lookup_array_d.gpudata, m_particle_array_d.gpudata, x_array_d.gpudata,
                              y_array_d.gpudata, z_array_d.gpudata, re_lookup_pl_d.gpudata, im_lookup_pl_d.gpudata,
                              re_lookup_mn_d.gpudata, im_lookup_mn_d.gpudata, re_in_vec_d.gpudata, im_in_vec_d.gpudata,
                              re_result_d.gpudata, im_result_d.gpudata, block=(cuda_blocksize, 1, 1),
                              grid=(cuda_gridsize, 1))
            return re_result_d.get() + 1j * im_result_d.get()

        self.linear_operator = scipy.sparse.linalg.LinearOperator(shape=self.shape, matvec=matvec, dtype=complex)


class CouplingMatrixRadialLookup(SystemMatrix):
    """Base class for radial lookup based coupling matrix either on CPU or on GPU (CUDA).

    Args:
        vacuum_wavelength (float): vacuum wavelength in length units
        particle_list (list): list of sumthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem): stratified medium
        k_parallel (numpy.ndarray or str): in-plane wavenumber.
                                           If 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
        resolution (float or None): spatial resolution of the lookup in the radial direction
    """

    def __init__(self, vacuum_wavelength, particle_list, layer_system, k_parallel='default', resolution=None):
        z_list = [particle.position[2] for particle in particle_list]
        assert z_list.count(z_list[0]) == len(z_list)

        SystemMatrix.__init__(self, particle_list)

        self.l_max = max([particle.l_max for particle in particle_list])
        self.m_max = max([particle.m_max for particle in particle_list])
        self.blocksize = flds.blocksize(self.l_max, self.m_max)
        self.resolution = resolution
        self.lookup_table, self.radial_distance_array = look.radial_coupling_lookup_table(
            vacuum_wavelength=vacuum_wavelength, particle_list=particle_list, layer_system=layer_system,
            k_parallel=k_parallel, resolution=resolution)


class CouplingMatrixRadialLookupCUDA(CouplingMatrixRadialLookup):
    """Radial lookup based coupling matrix either on GPU (CUDA).

    Args:
        vacuum_wavelength (float): vacuum wavelength in length units
        particle_list (list): list of sumthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem): stratified medium
        k_parallel (numpy.ndarray or str): in-plane wavenumber.
                                           If 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
        resolution (float or None): spatial resolution of the lookup in the radial direction
        cuda_blocksize (int): threads per block when calling CUDA kernel
    """

    def __init__(self, vacuum_wavelength, particle_list, layer_system, k_parallel='default', resolution=None,
                 cuda_blocksize=None, interpolator_kind='linear'):

        if cuda_blocksize is None:
            cuda_blocksize = cu.default_blocksize

        CouplingMatrixRadialLookup.__init__(self, vacuum_wavelength, particle_list, layer_system, k_parallel,
                                            resolution)

        sys.stdout.write('Prepare CUDA kernel and device lookup data ... ')
        sys.stdout.flush()
        start_time = time.time()

        if interpolator_kind == 'linear':
            coupling_source = cusrc.linear_radial_lookup_source % (self.blocksize, self.shape[0],
                                                                self.radial_distance_array.min(), resolution)
        elif interpolator_kind == 'cubic':
            coupling_source = cusrc.cubic_radial_lookup_source % (self.blocksize, self.shape[0],
                                                               self.radial_distance_array.min(), resolution)

        coupling_function = cu.SourceModule(coupling_source).get_function("coupling_kernel")

        n_lookup_array = np.zeros(self.shape[0], dtype=np.uint32)
        m_particle_array = np.zeros(self.shape[0], dtype=np.float32)
        x_array = np.zeros(self.shape[0], dtype=np.float32)
        y_array = np.zeros(self.shape[0], dtype=np.float32)

        i_particle = 0
        for i, particle in enumerate(particle_list):
            for m in range(-particle.m_max, particle.m_max + 1):
                for l in range(max(1, abs(m)), particle.l_max + 1):
                    for tau in range(2):
                        # idx = self.index(i, tau, l, m)
                        i_taulm = flds.multi_to_single_index(tau, l, m, particle.l_max, particle.m_max)
                        idx = i_particle + i_taulm

                        n_lookup_array[idx] = flds.multi_to_single_index(tau, l, m, self.l_max, self.m_max)
                        m_particle_array[idx] = m

                        # scale the x and y position to the lookup resolution:
                        x_array[idx] = particle.position[0]
                        y_array[idx] = particle.position[1]

            i_particle += flds.blocksize(particle.l_max, particle.m_max)

        # lookup as numpy array in required shape
        re_lookup = self.lookup_table.real.astype(np.float32)
        im_lookup = self.lookup_table.imag.astype(np.float32)

        # transfer data to gpu
        n_lookup_array_d = cu.gpuarray.to_gpu(n_lookup_array)
        m_particle_array_d = cu.gpuarray.to_gpu(m_particle_array)
        x_array_d = cu.gpuarray.to_gpu(x_array)
        y_array_d = cu.gpuarray.to_gpu(y_array)
        re_lookup_d = cu.gpuarray.to_gpu(re_lookup)
        im_lookup_d = cu.gpuarray.to_gpu(im_lookup)

        sys.stdout.write('done | elapsed: ' + str(int(time.time() - start_time)) + 's\n')
        sys.stdout.flush()

        cuda_gridsize = (self.shape[0] + cuda_blocksize - 1) // cuda_blocksize

        def matvec(in_vec):
            re_in_vec_d = cu.gpuarray.to_gpu(np.float32(in_vec.real))
            im_in_vec_d = cu.gpuarray.to_gpu(np.float32(in_vec.imag))
            re_result_d = cu.gpuarray.zeros(in_vec.shape, dtype=np.float32)
            im_result_d = cu.gpuarray.zeros(in_vec.shape, dtype=np.float32)
            coupling_function(n_lookup_array_d.gpudata, m_particle_array_d.gpudata, x_array_d.gpudata,
                              y_array_d.gpudata,
                              re_lookup_d.gpudata, im_lookup_d.gpudata, re_in_vec_d.gpudata, im_in_vec_d.gpudata,
                              re_result_d.gpudata, im_result_d.gpudata, block=(cuda_blocksize, 1, 1),
                              grid=(cuda_gridsize, 1))
            return re_result_d.get() + 1j * im_result_d.get()

        self.linear_operator = scipy.sparse.linalg.LinearOperator(shape=self.shape, matvec=matvec, dtype=complex)


class CouplingMatrixRadialLookupCPU(CouplingMatrixRadialLookup):
    """Class for radial lookup based coupling matrix running on CPU. This is used when no suitable GPU device is detected
    or when PyCuda is not installed.

    Args:
        vacuum_wavelength (float): vacuum wavelength in length units
        particle_list (list): list of sumthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem): stratified medium
        k_parallel (numpy.ndarray or str): in-plane wavenumber.
                                           If 'default', use smuthi.fields.default_Sommerfeld_k_parallel_array
        resolution (float or None): spatial resolution of the lookup in the radial direction
        kind (str): interpolation order, e.g. 'linear' or 'cubic'
    """

    def __init__(self, vacuum_wavelength, particle_list, layer_system, k_parallel='default', resolution=None,
                 interpolator_kind='linear'):

        z_list = [particle.position[2] for particle in particle_list]
        assert z_list.count(z_list[0]) == len(z_list)

        CouplingMatrixRadialLookup.__init__(self, vacuum_wavelength, particle_list, layer_system, k_parallel,
                                            resolution)

        x_array = np.array([particle.position[0] for particle in particle_list])
        y_array = np.array([particle.position[1] for particle in particle_list])

        self.particle_rho_array = np.sqrt((x_array[:, None] - x_array[None, :]) ** 2
                                          + (y_array[:, None] - y_array[None, :]) ** 2)
        self.particle_phi_array = np.arctan2(y_array[:, None] - y_array[None, :], x_array[:, None] - x_array[None, :])

        # contains for each n all positions in the large system arrays that correspond to n:
        self.system_vector_index_list = [[] for i in range(self.blocksize)]

        # same size as system_vector_index_list, contains the according particle numbers:
        self.particle_number_list = [[] for i in range(self.blocksize)]
        self.m_list = [None for i in range(self.blocksize)]
        for i, particle in enumerate(particle_list):
            for m in range(-particle.m_max, particle.m_max + 1):
                for l in range(max(1, abs(m)), particle.l_max + 1):
                    for tau in range(2):
                        n_lookup = flds.multi_to_single_index(tau=tau, l=l, m=m, l_max=self.l_max, m_max=self.m_max)
                        self.system_vector_index_list[n_lookup].append(self.index(i, tau, l, m))
                        self.particle_number_list[n_lookup].append(i)
                        self.m_list[n_lookup] = m
        for n in range(self.blocksize):
            self.system_vector_index_list[n] = np.array(self.system_vector_index_list[n])
            self.particle_number_list[n] = np.array(self.particle_number_list[n])

        lookup = [[None for i in range(self.blocksize)] for i2 in range(self.blocksize)]
        for n1 in range(self.blocksize):
            for n2 in range(self.blocksize):
                lookup[n1][n2] = scipy.interpolate.interp1d(x=self.radial_distance_array,
                                                            y=self.lookup_table[:, n1, n2],
                                                            kind=interpolator_kind, axis=-1, assume_sorted=True)

        def matvec(in_vec):
            out_vec = np.zeros(shape=in_vec.shape, dtype=complex)
            for n1 in range(self.blocksize):
                i1 = self.particle_number_list[n1]
                idx1 = self.system_vector_index_list[n1]
                m1 = self.m_list[n1]
                for n2 in range(self.blocksize):
                    i2 = self.particle_number_list[n2]
                    idx2 = self.system_vector_index_list[n2]
                    m2 = self.m_list[n2]
                    rho = self.particle_rho_array[i1[:, None], i2[None, :]]
                    phi = self.particle_phi_array[i1[:, None], i2[None, :]]
                    M = lookup[n1][n2](rho)
                    M = M * np.exp(1j * (m2 - m1) * phi)
                    out_vec[idx1] += M.dot(in_vec[idx2])
            return out_vec

        self.linear_operator = scipy.sparse.linalg.LinearOperator(shape=self.shape, matvec=matvec, dtype=complex)

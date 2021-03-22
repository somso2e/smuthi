# -*- coding: utf-8 -*-
"""Classes to manage the expansion of the electric field in plane wave and 
spherical wave basis sets."""

import numpy as np
import smuthi.fields
import smuthi.fields.vector_wave_functions as vwf
import smuthi.fields.expansions_cuda as cu_src
import smuthi.utility.cuda as cu
import smuthi.utility.numba_helpers as nh
import multiprocessing as mp
import smuthi.utility.multiprocessing_helpers as mp_helpers
from functools import partial
from enum import Enum
from psutil import virtual_memory
import platform
import copy
import math


class FieldExpansion:
    """Base class for field expansions."""

    def __init__(self):
        self.validity_conditions = []

    def valid(self, x, y, z):
        """Test if points are in definition range of the expansion. 
        Abstract method to be overwritten in child classes.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
         
        Returns:
            numpy.ndarray of bool datatype indicating if points are inside 
            definition domain.
        """
        ret = np.ones(x.shape, dtype=bool)
        for check in self.validity_conditions:
            ret = np.logical_and(ret, check(x, y, z))
        return ret

    def diverging(self, x, y, z):
        """Test if points are in domain where expansion could diverge. Virtual 
        method to be overwritten in child 
        classes.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
         
        Returns:
            numpy.ndarray of bool datatype indicating if points are inside 
            divergence domain.
        """
        pass

    def electric_field(self, x, y, z):
        """Evaluate electric field. Virtual method to be overwritten in child 
        classes.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
         
        Returns:
            Tuple of (E_x, E_y, E_z) numpy.ndarray objects with the Cartesian 
            coordinates of complex electric field.
        """
        pass
    
    def magnetic_field(self, x, y, z, vacuum_wavelength):
        """Evaluate magnetic field. Virtual method to be overwritten in child 
        classes.
        
        Args:
            x (numpy.ndarray):          x-coordinates of query points
            y (numpy.ndarray):          y-coordinates of query points
            z (numpy.ndarray):          z-coordinates of query points
            vacuum_wavelength (float):  Vacuum wavelength in length units
         
        Returns:
            Tuple of (H_x, H_y, H_z) numpy.ndarray objects with the Cartesian 
            coordinates of complex magnetic field.
        """
        pass


class PiecewiseFieldExpansion(FieldExpansion):
    r"""Manage a field that is expanded in different ways for different 
    domains, i.e., an expansion of the kind
    
    .. math::
        \mathbf{E}(\mathbf{r}) = \sum_{i} \mathbf{E}_i(\mathbf{r}),
    
    where
    
    .. math::
        \mathbf{E}_i(\mathbf{r}) = \begin{cases} \tilde{\mathbf{E}}_i(\mathbf{r}) & \text{ if }\mathbf{r}\in D_i \\ 0 & \text{ else} \end{cases}
    
    and :math:`\tilde{\mathbf{E_i}}(\mathbf{r})` is either a plane wave 
    expansion or a spherical wave expansion, and 
    :math:`D_i` is its domain of validity.
    """
    def __init__(self):
        FieldExpansion.__init__(self)
        self.expansion_list = []

    def valid(self, x, y, z):
        """Test if points are in definition range of the expansion.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
         
        Returns:
            numpy.ndarray of bool datatype indicating if points are inside 
            definition domain.
        """
        vld = np.zeros(x.shape, dtype=bool)
        for fex in self.expansion_list:
            vld = np.logical_or(vld, fex.valid(x, y, z))

        vld = np.logical_and(vld, FieldExpansion.valid(self, x, y, z))

        return vld

    def diverging(self, x, y, z):
        """Test if points are in domain where expansion could diverge.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
         
        Returns:
            numpy.ndarray of bool datatype indicating if points are inside 
            divergence domain.
        """
        dvg = np.zeros(x.shape, dtype=bool)
        for fex in self.expansion_list:
            dvg = np.logical_and(dvg, fex.diverging(x, y, z))
        return dvg
    
    def electric_field(self, x, y, z):
        """Evaluate electric field.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
         
        Returns:
            Tuple of (E_x, E_y, E_z) numpy.ndarray objects with the Cartesian 
            coordinates of complex electric field.
        """
        x, y, z = np.array(x, ndmin=1), np.array(y, ndmin=1), np.array(z, ndmin=1)
        ex = np.zeros(x.shape, dtype=complex)
        ey = np.zeros(x.shape, dtype=complex)
        ez = np.zeros(x.shape, dtype=complex)
        vld = self.valid(x, y, z)
        for fex in self.expansion_list:
            dex, dey, dez = fex.electric_field(x, y, z)
            ex[vld], ey[vld], ez[vld] = ex[vld] + dex[vld], ey[vld] + dey[vld], ez[vld] + dez[vld]
        return ex, ey, ez
    
    def magnetic_field(self, x, y, z, vacuum_wavelength):
        """Evaluate magnetic field. Virtual method to be overwritten in child 
        classes.
        
        Args:
            x (numpy.ndarray):          x-coordinates of query points
            y (numpy.ndarray):          y-coordinates of query points
            z (numpy.ndarray):          z-coordinates of query points
            vacuum_wavelength (float):  Vacuum wavelength in length units
         
        Returns:
            Tuple of (H_x, H_y, H_z) numpy.ndarray objects with the Cartesian 
            coordinates of complex magnetic field.
        """
        x, y, z = np.array(x, ndmin=1), np.array(y, ndmin=1), np.array(z, ndmin=1)
        hx = np.zeros(x.shape, dtype=complex)
        hy = np.zeros(x.shape, dtype=complex)
        hz = np.zeros(x.shape, dtype=complex)
        vld = self.valid(x, y, z)
        for fex in self.expansion_list:
            dhx, dhy, dhz = fex.magnetic_field(x, y, z, vacuum_wavelength)
            hx[vld], hy[vld], hz[vld] = hx[vld] + dhx[vld], hy[vld] + dhy[vld], hz[vld] + dhz[vld]
        return hx, hy, hz

    def compatible(self, other):
        """Returns always true, because any field expansion can be added to a 
        piecewise field expansion."""
        return True

    def __add__(self, other):
        """Addition of expansion objects.

        Args:
            other (FieldExpansion):  expansion object to add to this object

        Returns:
            PiecewiseFieldExpansion object as the sum of this expansion and the 
            other
        """
        # todo: testing
        pfe_sum = PiecewiseFieldExpansion()

        if type(other).__name__ == "PiecewiseFieldExpansion":
            added = [False for other_fex in other.expansion_list]
            for self_fex in self.expansion_list:
                fex = copy.deepcopy(self_fex)
                for i, other_fex in enumerate(other.expansion_list):
                    if (not added[i]) and self_fex.compatible(other_fex):
                        fex = fex + other_fex
                        added[i] = True
                pfe_sum.expansion_list.append(fex)
            for i, other_fex in enumerate(other.expansion_list):
                if not added[i]:
                    pfe_sum.expansion_list.append(other_fex)
        else:
            added = False
            for self_fex in self.expansion_list:
                fex = copy.deepcopy(self_fex)
                if (not added) and fex.compatible(other):
                    pfe_sum.expansion_list.append(fex + other)
                    added = True
                else:
                    pfe_sum.expansion_list.append(fex)
            if not added:
                pfe_sum.expansion_list.append(other)

        return pfe_sum


class SphericalWaveExpansion(FieldExpansion):
    r"""A class to manage spherical wave expansions of the form

    .. math::
        \mathbf{E}(\mathbf{r}) = \sum_{\tau=1}^2 \sum_{l=1}^\infty \sum_{m=-l}^l a_{\tau l m} 
        \mathbf{\Psi}^{(\nu)}_{\tau l m}(\mathbf{r} - \mathbf{r}_i)

    for :math:`\mathbf{r}` located in a layer defined by 
    :math:`z\in [z_{min}, z_{max}]`
    and where :math:`\mathbf{\Psi}^{(\nu)}_{\tau l m}` are the SVWFs, see 
    :meth:`smuthi.vector_wave_functions.spherical_vector_wave_function`.

    Internally, the expansion coefficients :math:`a_{\tau l m}` are stored as a 
    1-dimensional array running over a multi index :math:`n` subsumming over 
    the SVWF indices :math:`(\tau,l,m)`. The 
    mapping from the SVWF indices to the multi
    index is organized by the function :meth:`multi_to_single_index`.
    
    Args:
        k (float):    wavenumber in layer where expansion is valid
        l_max (int):  maximal multipole degree :math:`l_\mathrm{max}\geq 1` 
        where to truncate the expansion. m_max (int):  maximal multipole order 
        :math:`0 \leq m_\mathrm{max} \leq l_\mathrm{max}` where to truncate the 
        expansion.
        kind (str):   'regular' for :math:`\nu=1` or 'outgoing' for :math:`\nu=3`
        reference_point (list or tuple):  [x, y, z]-coordinates of point relative 
                                          to which the spherical waves are 
                                          considered (e.g., particle center).
        lower_z (float):   the expansion is valid on and above that z-coordinate
        upper_z (float):   the expansion is valid below that z-coordinate
        inner_r (float):   radius inside which the expansion diverges 
                           (e.g. circumscribing sphere of particle)
        outer_r (float):   radius outside which the expansion diverges

    Attributes:
        coefficients (numpy ndarray): expansion coefficients 
        :math:`a_{\tau l m}` ordered by multi index n
    """

    def __init__(self, k, l_max, m_max=None, kind=None, reference_point=None, lower_z=-np.inf, upper_z=np.inf,
                 inner_r=0, outer_r=np.inf):
        FieldExpansion.__init__(self)
        self.k = k
        self.l_max = l_max
        if m_max is not None:
            self.m_max = m_max
        else:
            self.m_max = l_max
        self.coefficients = np.zeros(smuthi.fields.blocksize(self.l_max, self.m_max), dtype=complex)
        self.kind = kind  # 'regular' or 'outgoing'
        self.reference_point = reference_point
        self.lower_z = lower_z
        self.upper_z = upper_z
        self.inner_r = inner_r
        self.outer_r = outer_r

    def valid(self, x, y, z):
        """Test if points are in definition range of the expansion.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
         
        Returns:
            numpy.ndarray of bool datatype indicating if points are inside 
            definition domain.
        """
        vld = np.logical_and(z >= self.lower_z, z < self.upper_z)
        return np.logical_and(vld, FieldExpansion.valid(self, x, y, z))

    def diverging(self, x, y, z):
        """Test if points are in domain where expansion could diverge.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
         
        Returns:
            numpy.ndarray of bool datatype indicating if points are inside 
            divergence domain.
        """
        r = np.sqrt((x - self.reference_point[0])**2 + (y - self.reference_point[1])**2
                    + (z - self.reference_point[2])**2)
        if self.kind == 'regular':
            return r >= self.outer_r
        if self.kind == 'outgoing':
            return r <= self.inner_r
        else:
            return None

    def coefficients_tlm(self, tau, l, m):
        """SWE coefficient for given (tau, l, m)

        Args:
            tau (int):  SVWF polarization (0 for spherical TE, 1 for spherical TM)
            l (int):    SVWF degree
            m (int):    SVWF order

        Returns:
            SWE coefficient
        """
        n = smuthi.fields.multi_to_single_index(tau, l, m, self.l_max, self.m_max)
        return self.coefficients[n]
    
    def electric_field(self, x, y, z):
        """Evaluate electric field.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
         
        Returns:
            Tuple of (E_x, E_y, E_z) numpy.ndarray objects with the Cartesian 
            coordinates of complex electric field.
        """
        x = np.array(x, ndmin=1)
        y = np.array(y, ndmin=1)
        z = np.array(z, ndmin=1)

        xr = x[self.valid(x, y, z)] - self.reference_point[0]
        yr = y[self.valid(x, y, z)] - self.reference_point[1]
        zr = z[self.valid(x, y, z)] - self.reference_point[2]
        ex = np.zeros(x.shape, dtype=complex)
        ey = np.zeros(x.shape, dtype=complex)
        ez = np.zeros(x.shape, dtype=complex)
        for tau in range(2):
            for m in range(-self.m_max, self.m_max + 1):
                for l in range(max(1, abs(m)), self.l_max + 1):
                    b = self.coefficients_tlm(tau, l, m)
                    if self.kind == 'regular':
                        Nx, Ny, Nz = vwf.spherical_vector_wave_function(xr, yr, zr, self.k, 1, tau, l, m)
                    elif self.kind == 'outgoing':
                        Nx, Ny, Nz = vwf.spherical_vector_wave_function(xr, yr, zr, self.k, 3, tau, l, m)
                    ex[self.valid(x, y, z)] += b * Nx
                    ey[self.valid(x, y, z)] += b * Ny
                    ez[self.valid(x, y, z)] += b * Nz
        return ex, ey, ez
    
    def magnetic_field(self, x, y, z, vacuum_wavelength):
        """Evaluate magnetic field.
        
        Args:
            x (numpy.ndarray):          x-coordinates of query points
            y (numpy.ndarray):          y-coordinates of query points
            z (numpy.ndarray):          z-coordinates of query points
            vacuum_wavelength (float):  Vacuum wavelength in length units
         
        Returns:
            Tuple of (H_x, H_y, H_z) numpy.ndarray objects with the Cartesian
            coordinates of complex electric field.
        """
        omega = smuthi.fields.angular_frequency(vacuum_wavelength)
        
        x = np.array(x, ndmin=1)
        y = np.array(y, ndmin=1)
        z = np.array(z, ndmin=1)

        xr = x[self.valid(x, y, z)] - self.reference_point[0]
        yr = y[self.valid(x, y, z)] - self.reference_point[1]
        zr = z[self.valid(x, y, z)] - self.reference_point[2]
        hx = np.zeros(x.shape, dtype=complex)
        hy = np.zeros(x.shape, dtype=complex)
        hz = np.zeros(x.shape, dtype=complex)
        for tau in range(2):
            for m in range(-self.m_max, self.m_max + 1):
                for l in range(max(1, abs(m)), self.l_max + 1):
                    b = self.coefficients_tlm(1-tau, l, m)
                    if self.kind == 'regular':
                        Nx, Ny, Nz = vwf.spherical_vector_wave_function(xr, yr, zr, self.k, 1, tau, l, m)
                    elif self.kind == 'outgoing':
                        Nx, Ny, Nz = vwf.spherical_vector_wave_function(xr, yr, zr, self.k, 3, tau, l, m)
                    hx[self.valid(x, y, z)] += b * Nx
                    hy[self.valid(x, y, z)] += b * Ny
                    hz[self.valid(x, y, z)] += b * Nz
        
        hx = - 1j * self.k / omega * hx
        hy = - 1j * self.k / omega * hy
        hz = - 1j * self.k / omega * hz      
        
        return hx, hy, hz

    def compatible(self, other):
        """Check if two spherical wave expansions are compatible in the sense 
        that they can be added coefficient-wise

        Args:
            other (FieldExpansion):  expansion object to add to this object

        Returns:
            bool (true if compatible, false else)
        """
        return (type(other).__name__ == "SphericalWaveExpansion" 
                and self.k == other.k 
                and self.l_max == other.l_max
                and self.m_max == other.m_max 
                and self.kind == other.kind
                and self.reference_point == other.reference_point)

    def __add__(self, other):
        """Addition of expansion objects (by coefficients).
        
        Args:
            other (SphericalWaveExpansion):  expansion object to add to this object
        
        Returns:
            SphericalWaveExpansion object as the sum of this expansion and the other
        """
        # todo: allow different l_max
        if not self.compatible(other):
            raise ValueError('SphericalWaveExpansions are inconsistent.')
        swe_sum = SphericalWaveExpansion(k=self.k, l_max=self.l_max, m_max=self.m_max, kind=self.kind,
                                         reference_point=self.reference_point, inner_r=max(self.inner_r, other.inner_r),
                                         outer_r=min(self.outer_r, other.outer_r),
                                         lower_z=max(self.lower_z, other.lower_z),
                                         upper_z=min(self.upper_z, other.upper_z))
        swe_sum.coefficients = self.coefficients + other.coefficients
        return swe_sum


class PlaneWaveExpansion(FieldExpansion):
    r"""A class to manage plane wave expansions of the form

    .. math::
        \mathbf{E}(\mathbf{r}) = \sum_{j=1}^2 \iint \mathrm{d}^2\mathbf{k}_\parallel \, g_{j}(\kappa, \alpha)
        \mathbf{\Phi}^\pm_j(\kappa, \alpha; \mathbf{r} - \mathbf{r}_i)

    for :math:`\mathbf{r}` located in a layer defined by :math:`z\in [z_{min}, z_{max}]`
    and :math:`\mathrm{d}^2\mathbf{k}_\parallel = \kappa\,\mathrm{d}\alpha\,\mathrm{d}\kappa`. 
    
    The double integral runs over :math:`\alpha\in[0, 2\pi]` and :math:`\kappa\in[0,\kappa_\mathrm{max}]`. 
    Further, :math:`\mathbf{\Phi}^\pm_j` are the PVWFs, see :meth:`plane_vector_wave_function`.

    Internally, the expansion coefficients :math:`g_{ij}^\pm(\kappa, \alpha)` 
    are stored as a 3-dimensional numpy ndarray.
    
    If the attributes k_parallel and azimuthal_angles have only a single entry, 
    a discrete distribution is assumed:

    .. math::
        g_{j}^\pm(\kappa, \alpha) \sim \delta^2(\mathbf{k}_\parallel - \mathbf{k}_{\parallel, 0})

    .. todo: update attributes doc

    Args:
        k (float):                          wavenumber in layer where expansion is valid
        k_parallel (numpy ndarray):         array of in-plane wavenumbers (can be float or complex)
        azimuthal_angles (numpy ndarray):   :math:`\alpha`, from 0 to :math:`2\pi`
        kind (str):                         'upgoing' for :math:`g^+` and 'downgoing' for :math:`g^-` type
                                            expansions 
        reference_point (list or tuple):    [x, y, z]-coordinates of point relative to which the plane waves are 
                                            defined.
        lower_z (float):                    the expansion is valid on and above that z-coordinate
        upper_z (float):                    the expansion is valid below that z-coordinate
        

    Attributes:
        coefficients (numpy ndarray): coefficients[j, k, l] contains 
        :math:`g^\pm_{j}(\kappa_{k}, \alpha_{l})`
    """
    def __init__(self, k, k_parallel, azimuthal_angles, kind=None, reference_point=None, lower_z=-np.inf,
                 upper_z=np.inf):
        FieldExpansion.__init__(self)
        self.k = k
        self.k_parallel = np.array(k_parallel, ndmin=1)
        self.azimuthal_angles = np.array(azimuthal_angles, ndmin=1)
        self.kind = kind  # 'upgoing' or 'downgoing'
        self.reference_point = reference_point
        self.lower_z = lower_z
        self.upper_z = upper_z

        # The coefficients :math:`g^\pm_{j}(\kappa,\alpha) are represented as a 3-dimensional numpy.ndarray.
        # The indices are:
        # - polarization (0=TE, 1=TM)
        # - index of the kappa dimension
        # - index of the alpha dimension
        self.coefficients = np.zeros((2, len(self.k_parallel), len(self.azimuthal_angles)), dtype=complex)

    def valid(self, x, y, z):
        """Test if points are in definition range of the expansion.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
         
        Returns:
            numpy.ndarray of bool datatype indicating if points are inside 
            definition domain.
        """
        if self.upper_z == float('inf'):
            vld = np.logical_and(z >= self.lower_z, 1)
        else:
            vld = np.logical_and(z >= self.lower_z, z < self.upper_z)
        vld_custom = FieldExpansion.valid(self, x, y, z)
        return np.logical_and(vld, vld_custom)

    def diverging(self, x, y, z):
        """Test if points are in domain where expansion could diverge.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
         
        Returns:
            numpy.ndarray of bool datatype indicating if points are inside 
            divergence domain.
        """
        return np.zeros(x.shape,dtype=bool)

    def k_parallel_grid(self):
        """Meshgrid of n_effective with respect to azimuthal_angles"""
        kp_grid, _ = np.meshgrid(self.k_parallel, self.azimuthal_angles, indexing='ij')
        return kp_grid

    def azimuthal_angle_grid(self):
        """Meshgrid of azimuthal_angles with respect to n_effective"""
        _, a_grid = np.meshgrid(self.k_parallel, self.azimuthal_angles, indexing='ij')
        return a_grid

    def k_z(self):
        if self.kind == 'upgoing':
            kz = smuthi.fields.k_z(k_parallel=self.k_parallel, k=self.k)
        elif self.kind == 'downgoing':
            kz = -smuthi.fields.k_z(k_parallel=self.k_parallel, k=self.k)
        else:
            raise ValueError('pwe kind undefined')
        return kz

    def k_z_grid(self):
        if self.kind == 'upgoing':
            kz = smuthi.fields.k_z(k_parallel=self.k_parallel_grid(), k=self.k)
        elif self.kind == 'downgoing':
            kz = -smuthi.fields.k_z(k_parallel=self.k_parallel_grid(), k=self.k)
        else:
            raise ValueError('pwe type undefined')
        return kz

    def compatible(self, other):
        """Check if two plane wave expansions are compatible in the sense that 
        they can be added coefficient-wise

        Args:
            other (FieldExpansion):  expansion object to add to this object

        Returns:
            bool (true if compatible, false else)
        """
        return (type(other).__name__=="PlaneWaveExpansion" and np.isclose(self.k, other.k)
                and all(np.isclose(self.k_parallel, other.k_parallel))
                and all(np.isclose(self.azimuthal_angles, other.azimuthal_angles)) and self.kind == other.kind
                and self.reference_point == other.reference_point)

    def __add__(self, other):
        if not self.compatible(other):
            raise ValueError('Plane wave expansion are inconsistent.')
        pwe_sum = PlaneWaveExpansion(k=self.k, k_parallel=self.k_parallel, azimuthal_angles=self.azimuthal_angles,
                                     kind=self.kind, reference_point=self.reference_point,
                                     lower_z=max(self.lower_z, other.lower_z),
                                     upper_z=min(self.upper_z, other.upper_z))
        pwe_sum.coefficients = self.coefficients + other.coefficients
        return pwe_sum
    
    def electric_field(self, x, y, z, max_chunksize=50, cpu_precision='single precision'):
        """Evaluate electric field.
        
        Args:
            x (numpy.ndarray):    x-coordinates of query points
            y (numpy.ndarray):    y-coordinates of query points
            z (numpy.ndarray):    z-coordinates of query points
            max_chunksize (int):  max number of field points that are simultaneously 
                                  evaluated when running in CPU mode. 
                                  In Windows/MacOS max_chunksize = chunksize, 
                                  in Linux it can be decreased considering available CPU cores.
            cpu_precision (string): set 'double precision' to use float64 and complex128 types
                                    instead of float32 and complex64
          
        Returns:
            Tuple of (E_x, E_y, E_z) numpy.ndarray objects with the Cartesian 
            coordinates of complex electric field.
        """
        ex = np.zeros(x.shape, dtype=complex)
        ey = np.zeros(x.shape, dtype=complex)
        ez = np.zeros(x.shape, dtype=complex)

        xr = x[self.valid(x, y, z)] - self.reference_point[0]
        yr = y[self.valid(x, y, z)] - self.reference_point[1]
        zr = z[self.valid(x, y, z)] - self.reference_point[2]
        
        if cu.use_gpu and xr.size and len(self.k_parallel) > 1:  # run calculations on gpu

            re_kp_d = cu.gpuarray.to_gpu(self.k_parallel.real.astype(np.float32))
            im_kp_d = cu.gpuarray.to_gpu(self.k_parallel.imag.astype(np.float32))

            re_kz_d = cu.gpuarray.to_gpu(self.k_z().real.astype(np.float32))
            im_kz_d = cu.gpuarray.to_gpu(self.k_z().imag.astype(np.float32))

            alpha_d = cu.gpuarray.to_gpu(self.azimuthal_angles.astype(np.float32))

            xr_d = cu.gpuarray.to_gpu(xr.astype(np.float32))
            yr_d = cu.gpuarray.to_gpu(yr.astype(np.float32))
            zr_d = cu.gpuarray.to_gpu(zr.astype(np.float32))

            re_g_te_d = cu.gpuarray.to_gpu(self.coefficients[0, :, :].real.astype(np.float32))
            im_g_te_d = cu.gpuarray.to_gpu(self.coefficients[0, :, :].imag.astype(np.float32))
            re_g_tm_d = cu.gpuarray.to_gpu(self.coefficients[1, :, :].real.astype(np.float32))
            im_g_tm_d = cu.gpuarray.to_gpu(self.coefficients[1, :, :].imag.astype(np.float32))

            re_e_x_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
            im_e_x_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
            re_e_y_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
            im_e_y_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
            re_e_z_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
            im_e_z_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))

            kernel_source = cu_src.pwe_electric_field_evaluation_code%(xr.size, len(self.k_parallel),
                                                                   len(self.azimuthal_angles), (1/self.k).real,
                                                                   (1/self.k).imag)

            kernel_function = cu.SourceModule(kernel_source).get_function("electric_field")
            cuda_blocksize = 128
            cuda_gridsize = (xr.size + cuda_blocksize - 1) // cuda_blocksize
            
            kernel_function(re_kp_d, im_kp_d, re_kz_d, im_kz_d, alpha_d, xr_d, yr_d, zr_d, re_g_te_d, im_g_te_d,
                            re_g_tm_d, im_g_tm_d, re_e_x_d, im_e_x_d, re_e_y_d, im_e_y_d, re_e_z_d, im_e_z_d,
                            block=(cuda_blocksize,1,1), grid=(cuda_gridsize,1))
            
            ex[self.valid(x, y, z)] = re_e_x_d.get() + 1j * im_e_x_d.get()
            ey[self.valid(x, y, z)] = re_e_y_d.get() + 1j * im_e_y_d.get()
            ez[self.valid(x, y, z)] = re_e_z_d.get() + 1j * im_e_z_d.get()
            
        else:  # run calculations on cpu
            ex[self.valid(x, y, z)], ey[self.valid(x, y, z)], ez[self.valid(x, y, z)] = \
                    self.__process_field_by_cpu(xr, yr, zr, max_chunksize, cpu_precision, 1,
                        self.__get_electric_field_integrands)

        return ex, ey, ez    

    def magnetic_field(self, x, y, z, vacuum_wavelength, max_chunksize=50, cpu_precision='single precision'):
        """Evaluate magnetic field.
        
        Args:
            x (numpy.ndarray):          x-coordinates of query points
            y (numpy.ndarray):          y-coordinates of query points
            z (numpy.ndarray):          z-coordinates of query points
            vacuum_wavelength (float):  Vacuum wavelength in length units
            chunksize (int):            number of field points that are simultaneously 
                                        evaluated when running in CPU mode
        Returns:
            Tuple of (H_x, H_y, H_z) numpy.ndarray objects with the Cartesian
            coordinates of complex magnetic field.
        """
        # todo: replace chunksize argument by automatic estimate (considering available RAM)
        hx = np.zeros(x.shape, dtype=complex)
        hy = np.zeros(x.shape, dtype=complex)
        hz = np.zeros(x.shape, dtype=complex)

        xr = x[self.valid(x, y, z)] - self.reference_point[0]
        yr = y[self.valid(x, y, z)] - self.reference_point[1]
        zr = z[self.valid(x, y, z)] - self.reference_point[2]
        
        omega = smuthi.fields.angular_frequency(vacuum_wavelength)
        
        if cu.use_gpu and xr.size and len(self.k_parallel) > 1:  # run calculations on gpu
            
            re_kp_d = cu.gpuarray.to_gpu(self.k_parallel.real.astype(np.float32))
            im_kp_d = cu.gpuarray.to_gpu(self.k_parallel.imag.astype(np.float32))
            
            re_kz_d = cu.gpuarray.to_gpu(self.k_z().real.astype(np.float32))
            im_kz_d = cu.gpuarray.to_gpu(self.k_z().imag.astype(np.float32))
            
            alpha_d = cu.gpuarray.to_gpu(self.azimuthal_angles.astype(np.float32))
            
            xr_d = cu.gpuarray.to_gpu(xr.astype(np.float32))
            yr_d = cu.gpuarray.to_gpu(yr.astype(np.float32))
            zr_d = cu.gpuarray.to_gpu(zr.astype(np.float32))
            
            re_g_te_d = cu.gpuarray.to_gpu(self.coefficients[0, :, :].real.astype(np.float32))
            im_g_te_d = cu.gpuarray.to_gpu(self.coefficients[0, :, :].imag.astype(np.float32))
            re_g_tm_d = cu.gpuarray.to_gpu(self.coefficients[1, :, :].real.astype(np.float32))
            im_g_tm_d = cu.gpuarray.to_gpu(self.coefficients[1, :, :].imag.astype(np.float32))
            
            re_h_x_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
            im_h_x_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
            re_h_y_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
            im_h_y_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
            re_h_z_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
            im_h_z_d = cu.gpuarray.to_gpu(np.zeros(xr.shape, dtype=np.float32))
                        
            kernel_source = cu_src.pwe_magnetic_field_evaluation_code%(xr.size, len(self.k_parallel), 
                                                                   len(self.azimuthal_angles), (self.k).real, 
                                                                   (self.k).imag)
            
            kernel_function = cu.SourceModule(kernel_source).get_function("magnetic_field") 
            cuda_blocksize = 128
            cuda_gridsize = (xr.size + cuda_blocksize - 1) // cuda_blocksize
            
            kernel_function(re_kp_d, im_kp_d, re_kz_d, im_kz_d, alpha_d, xr_d, yr_d, zr_d, re_g_te_d, im_g_te_d,
                            re_g_tm_d, im_g_tm_d, re_h_x_d, im_h_x_d, re_h_y_d, im_h_y_d, re_h_z_d, im_h_z_d,
                            block=(cuda_blocksize,1,1), grid=(cuda_gridsize,1))
          
            hx[self.valid(x, y, z)] = 1 / omega * (re_h_x_d.get() + 1j * im_h_x_d.get())
            hy[self.valid(x, y, z)] = 1 / omega * (re_h_y_d.get() + 1j * im_h_y_d.get())
            hz[self.valid(x, y, z)] = 1 / omega * (re_h_z_d.get() + 1j * im_h_z_d.get())
            
        else:  # run calculations on cpu
            hx[self.valid(x, y, z)], hy[self.valid(x, y, z)], hz[self.valid(x, y, z)] = \
                    self.__process_field_by_cpu(xr, yr, zr, max_chunksize, cpu_precision, omega,
                        self.__get_magnetic_field_integrands)

        return hx, hy, hz

    def __process_field_by_cpu(self, xr, yr, zr, max_chunksize, cpu_precision, omega, process_integrands):
        chunksize = self.__get_chunksize(max_chunksize, cpu_precision, xr.size)
                
        float_type = np.float32
        complex_type = np.complex64
        if cpu_precision == 'double precision':
            float_type = np.float64
            complex_type = np.complex128

        kpgrid = self.k_parallel_grid().astype(complex_type)
        agrid = self.azimuthal_angle_grid().astype(float_type)
        kx = kpgrid * np.cos(agrid)
        ky = kpgrid * np.sin(agrid)
        kz = self.k_z_grid().astype(complex_type)

        xr = xr.astype(float_type)
        yr = yr.astype(float_type)
        zr = zr.astype(float_type)

        f_x_flat = np.zeros(xr.size, dtype=complex_type)
        f_y_flat = np.zeros(xr.size, dtype=complex_type)
        f_z_flat = np.zeros(xr.size, dtype=complex_type)

        integrand_x, integrand_y, integrand_z = process_integrands(kz, agrid, kpgrid, complex_type)

        process_field_slice_method_with_context = partial(self.__process_field_slice_and_put_into_result, 
                        chunksize=chunksize,
                        xr=xr, yr=yr, zr=zr,
                        complex_type=complex_type,
                        integrand_x=integrand_x, integrand_y=integrand_y, integrand_z = integrand_z,
                        pwe=self,
                        kx = kx, ky = ky, kz = kz, omega = omega)

        results = []

        if self.__is_os_linux():
            # linux os.fork() works fine, so method "__process_field_slice_and_put_into_result"
            # can be paralleled from outside.
            # in Win and MasOS we have to parallel submethods in numba_helpers,
            # because otherwise it works incorrect.
            put_into_results = lambda results_to_be_filled, field_slices: results_to_be_filled.put(field_slices)
            results_q = mp.Queue()
            processes = []

            for i_chunk in range(math.ceil(xr.size / chunksize)):
                p = mp.Process(target=process_field_slice_method_with_context,
                            args = (i_chunk, results_q, put_into_results, self.OptimizationMethodsForLinux))
                processes.append(p)

            processes_clusters = mp_helpers.distribute_processes_into_clusters(processes, mp.cpu_count())

            for cluster in processes_clusters:
                cluster.execute()

                results = results + [results_q.get() for c in cluster.processes]

        else:
            put_into_results = lambda results_to_be_filled, field_slices: results_to_be_filled.append(field_slices)
            for i_chunk in range(math.ceil(xr.size / chunksize)):
                process_field_slice_method_with_context(i_chunk, results, put_into_results,
                                                        self.OptimizationMethodsFor_Not_Linux)

        self.__fill_flattened_arrays_by_field_slices(results, f_x_flat, f_y_flat, f_z_flat)

        return f_x_flat.reshape(xr.shape), f_y_flat.reshape(xr.shape), f_z_flat.reshape(xr.shape)

    def __get_chunksize(self, max_chunksize, cpu_precision, max_sensble_chunksize):
        reserved_ram_coefficient = 0.6
        ram_for_chunksize = virtual_memory().free * reserved_ram_coefficient

        memory_per_complex_value = 8
        if cpu_precision == 'double precision':
            memory_per_complex_value = 16

        available_parallel_processes = 1
        if self.__is_os_linux():
            available_parallel_processes = mp.cpu_count()

        kpgrid = self.k_parallel_grid()
        available_chunksize = int(ram_for_chunksize /
                        (4 * kpgrid.size * memory_per_complex_value * available_parallel_processes))
        # 4 because integrand_x_eikr, integrand_y_eikr, integrand_z_eikr and kr are the main memory consumers.       

        sensible_chunksize_per_process = int(max_sensble_chunksize / available_parallel_processes) + 1

        do_we_have_enough_memory = available_chunksize > max_chunksize
        do_we_have_much_stuff_to_process = sensible_chunksize_per_process > max_chunksize
        should_max_chunk_size_limit_everything = do_we_have_enough_memory and do_we_have_much_stuff_to_process

        if should_max_chunk_size_limit_everything:
            return max_chunksize

        should_we_break_processing_on_very_small_pieces = available_chunksize > sensible_chunksize_per_process
        if should_we_break_processing_on_very_small_pieces:
            return sensible_chunksize_per_process

        return available_chunksize

    def __get_electric_field_integrands(self, kz, agrid, kpgrid, complex_type):
        #pol=0
        integrand_x = (-np.sin(agrid) * self.coefficients[0, :, :]).astype(complex_type)[None, :, :]
        integrand_y = (np.cos(agrid) * self.coefficients[0, :, :]).astype(complex_type)[None, :, :]

        #pol=1
        integrand_x += (np.cos(agrid) * kz / self.k * self.coefficients[1, :, :]).astype(complex_type)[None, :, :]
        integrand_y += (np.sin(agrid) * kz / self.k * self.coefficients[1, :, :]).astype(complex_type)[None, :, :]
        integrand_z = (-kpgrid / self.k * self.coefficients[1, :, :]).astype(complex_type)[None, :, :]

        return integrand_x, integrand_y, integrand_z

    def __get_magnetic_field_integrands(self, kz, agrid, kpgrid, complex_type):
        #pol=0
        integrand_x = (-kz * np.cos(agrid) * self.coefficients[0, :, :]).astype(complex_type)[None, :, :]
        integrand_y = (-kz * np.sin(agrid) * self.coefficients[0, :, :]).astype(complex_type)[None, :, :]
        integrand_z = (kpgrid * self.coefficients[0, :, :]).astype(complex_type)[None, :, :]

        #pol=1
        integrand_x += (- np.sin(agrid) * self.k * self.coefficients[1, :, :]).astype(complex_type)[None, :, :]
        integrand_y += (np.cos(agrid) * self.k * self.coefficients[1, :, :]).astype(complex_type)[None, :, :]

        return integrand_x, integrand_y, integrand_z

    @staticmethod
    def __process_field_slice_and_put_into_result(i_chunk, results, put_into_results, optimization_methods,
                        chunksize, xr, yr, zr, complex_type,
                        integrand_x, integrand_y, integrand_z, pwe,
                        kx, ky, kz, omega = 1):
            chunk_idcs = range(i_chunk * chunksize, min((i_chunk + 1) * chunksize, xr.size))
            xr_chunk = xr.flatten()[chunk_idcs]
            yr_chunk = yr.flatten()[chunk_idcs]
            zr_chunk = zr.flatten()[chunk_idcs]

            numba_3tensordots_1dim_times_2dim = optimization_methods.numba_3tensordots_1dim_times_2dim
            kr = numba_3tensordots_1dim_times_2dim(xr_chunk, yr_chunk, zr_chunk, kx, ky, kz)

            evaluate_r_times_eikr = optimization_methods.evaluate_r_times_eikr
            integrand_x_eikr, integrand_y_eikr, integrand_z_eikr = evaluate_r_times_eikr(
                                                                    integrand_x, integrand_y, integrand_z, kr)

            res_x = np.zeros((chunk_idcs.stop - chunk_idcs.start), dtype=complex_type)
            res_y = np.zeros((chunk_idcs.stop - chunk_idcs.start), dtype=complex_type)
            res_z = np.zeros((chunk_idcs.stop - chunk_idcs.start), dtype=complex_type)

            if len(pwe.k_parallel) > 1:
                numba_trapz_3dim_array = optimization_methods.numba_trapz_3dim_array

                res_x = 1 / omega * np.trapz(
                    numba_trapz_3dim_array(integrand_x_eikr, pwe.azimuthal_angles)
                    * pwe.k_parallel, pwe.k_parallel)

                res_y = 1 / omega * np.trapz(
                    numba_trapz_3dim_array(integrand_y_eikr, pwe.azimuthal_angles)
                    * pwe.k_parallel, pwe.k_parallel)
                        
                res_z = 1 / omega * np.trapz(
                    numba_trapz_3dim_array(integrand_z_eikr, pwe.azimuthal_angles)
                    * pwe.k_parallel, pwe.k_parallel)
            else:
                res_x = 1 / omega * np.squeeze(integrand_x_eikr)
                res_y = 1 / omega * np.squeeze(integrand_y_eikr)
                res_z = 1 / omega * np.squeeze(integrand_z_eikr)

            put_into_results(results, (pwe.RawSliceOfField('x', chunk_idcs, res_x), \
                    pwe.RawSliceOfField('y', chunk_idcs, res_y), \
                    pwe.RawSliceOfField('z', chunk_idcs, res_z)))

    class RawSliceOfField:
        def __init__(self, axis, chunks, values):
            self.axis = axis
            self.chunks = chunks
            self.values = values

    @staticmethod
    def __fill_flattened_arrays_by_field_slices(results, f_x_flat, f_y_flat, f_z_flat):
        extracted_results = []

        for single_tuple in results:
            for parallel_result in single_tuple:
                extracted_results.append(parallel_result) # TODO Raplace by less ugly implementation.

        # get list where each element contains 'x'/'y'/'z' in 'axis' field
        results_x = list(filter(lambda x: 'x' in x.axis, extracted_results))
        results_y = list(filter(lambda x: 'y' in x.axis, extracted_results))
        results_z = list(filter(lambda x: 'z' in x.axis, extracted_results))

        for result in results_x:
            f_x_flat[result.chunks] = result.values
        for result in results_y:
            f_y_flat[result.chunks] = result.values
        for result in results_z:
            f_z_flat[result.chunks] = result.values

    def __is_os_linux(self):
        return 'Linux' in platform.system()

    class OptimizationMethodsForLinux(Enum):
        @staticmethod
        def __evaluate_r_times_eikr(integrand_x, integrand_y, integrand_z, kr):
            eikr = np.exp(1j * kr)

            integrand_x_eikr = integrand_x * eikr
            integrand_y_eikr = integrand_y * eikr
            integrand_z_eikr = integrand_z * eikr

            return integrand_x_eikr, integrand_y_eikr, integrand_z_eikr

        numba_3tensordots_1dim_times_2dim = nh.numba_3tensordots_1dim_times_2dim
        evaluate_r_times_eikr = __evaluate_r_times_eikr
        numba_trapz_3dim_array = nh.numba_trapz_3dim_array

    class OptimizationMethodsFor_Not_Linux(Enum):
        numba_3tensordots_1dim_times_2dim = nh.numba_3tensordots_1dim_times_2dim_parallel
        evaluate_r_times_eikr = nh.evaluate_r_times_eikr
        numba_trapz_3dim_array = nh.numba_trapz_3dim_array_parallel

    def _electric_field_cpu_legacy(self, x, y, z, chunksize=50):
        """This method is the previous method to compute the electric field of
        a plane wave expansion on the cpu. We still keep it for testing
        purpose."""
        ex = np.zeros(x.shape, dtype=complex)
        ey = np.zeros(x.shape, dtype=complex)
        ez = np.zeros(x.shape, dtype=complex)

        abc = self.valid(x, y, z)
        xr = x[self.valid(x, y, z)] - self.reference_point[0]
        yr = y[self.valid(x, y, z)] - self.reference_point[1]
        zr = z[self.valid(x, y, z)] - self.reference_point[2]

        kpgrid = self.k_parallel_grid()
        agrid = self.azimuthal_angle_grid()
        kx = kpgrid * np.cos(agrid)
        ky = kpgrid * np.sin(agrid)
        kz = self.k_z_grid()

        e_x_flat = np.zeros(xr.size, dtype=np.complex64)
        e_y_flat = np.zeros(xr.size, dtype=np.complex64)
        e_z_flat = np.zeros(xr.size, dtype=np.complex64)

        for i_chunk in range(math.ceil(xr.size / chunksize)):
            chunk_idcs = range(i_chunk * chunksize,
                               min((i_chunk + 1) * chunksize, xr.size))
            xr_chunk = xr.flatten()[chunk_idcs]
            yr_chunk = yr.flatten()[chunk_idcs]
            zr_chunk = zr.flatten()[chunk_idcs]

            kr = np.zeros((len(xr_chunk), len(self.k_parallel),
                           len(self.azimuthal_angles)), dtype=np.complex64)
            kr += np.tensordot(xr_chunk, kx, axes=0)
            kr += np.tensordot(yr_chunk, ky, axes=0)
            kr += np.tensordot(zr_chunk, kz, axes=0)

            eikr = np.exp(1j * kr)

            integrand_x = np.zeros((len(xr_chunk), len(self.k_parallel),
                                    len(self.azimuthal_angles)),
                                   dtype=np.complex64)
            integrand_y = np.zeros((len(yr_chunk), len(self.k_parallel),
                                    len(self.azimuthal_angles)),
                                   dtype=np.complex64)
            integrand_z = np.zeros((len(zr_chunk), len(self.k_parallel),
                                    len(self.azimuthal_angles)),
                                   dtype=np.complex64)

            # pol=0
            integrand_x += (-np.sin(agrid) * self.coefficients[0, :, :])[None,
                           :, :] * eikr
            integrand_y += (np.cos(agrid) * self.coefficients[0, :, :])[None, :,
                           :] * eikr
            # pol=1
            integrand_x += (np.cos(agrid) * kz / self.k * self.coefficients[1,
                                                          :, :])[None, :,
                           :] * eikr
            integrand_y += (np.sin(agrid) * kz / self.k * self.coefficients[1,
                                                          :, :])[None, :,
                           :] * eikr
            integrand_z += (-kpgrid / self.k * self.coefficients[1, :, :])[None,
                           :, :] * eikr

            if len(self.k_parallel) > 1:
                e_x_flat[chunk_idcs] = np.trapz(
                    np.trapz(integrand_x, self.azimuthal_angles)
                    * self.k_parallel, self.k_parallel)
                e_y_flat[chunk_idcs] = np.trapz(
                    np.trapz(integrand_y, self.azimuthal_angles)
                    * self.k_parallel, self.k_parallel)
                e_z_flat[chunk_idcs] = np.trapz(
                    np.trapz(integrand_z, self.azimuthal_angles)
                    * self.k_parallel, self.k_parallel)
            else:
                e_x_flat[chunk_idcs] = np.squeeze(integrand_x)
                e_y_flat[chunk_idcs] = np.squeeze(integrand_y)
                e_z_flat[chunk_idcs] = np.squeeze(integrand_z)

        ex[self.valid(x, y, z)] = e_x_flat.reshape(xr.shape)
        ey[self.valid(x, y, z)] = e_y_flat.reshape(xr.shape)
        ez[self.valid(x, y, z)] = e_z_flat.reshape(xr.shape)

        return ex, ey, ez

    def _magnetic_field_cpu_legacy(self, x, y, z, vacuum_wavelength, chunksize=50):
        """This method is the previous method to compute the magnetic field of
        a plane wave expansion on the cpu. We still keep it for testing
        purpose."""
        hx = np.zeros(x.shape, dtype=complex)
        hy = np.zeros(x.shape, dtype=complex)
        hz = np.zeros(x.shape, dtype=complex)

        xr = x[self.valid(x, y, z)] - self.reference_point[0]
        yr = y[self.valid(x, y, z)] - self.reference_point[1]
        zr = z[self.valid(x, y, z)] - self.reference_point[2]

        omega = smuthi.fields.angular_frequency(vacuum_wavelength)

        kpgrid = self.k_parallel_grid()
        agrid = self.azimuthal_angle_grid()
        kx = kpgrid * np.cos(agrid)
        ky = kpgrid * np.sin(agrid)
        kz = self.k_z_grid()

        h_x_flat = np.zeros(xr.size, dtype=np.complex64)
        h_y_flat = np.zeros(xr.size, dtype=np.complex64)
        h_z_flat = np.zeros(xr.size, dtype=np.complex64)

        for i_chunk in range(math.ceil(xr.size / chunksize)):
            chunk_idcs = range(i_chunk * chunksize,
                               min((i_chunk + 1) * chunksize, xr.size))
            xr_chunk = xr.flatten()[chunk_idcs]
            yr_chunk = yr.flatten()[chunk_idcs]
            zr_chunk = zr.flatten()[chunk_idcs]

            kr = np.zeros((len(xr_chunk), len(self.k_parallel),
                           len(self.azimuthal_angles)), dtype=np.complex64)
            kr += np.tensordot(xr_chunk, kx, axes=0)
            kr += np.tensordot(yr_chunk, ky, axes=0)
            kr += np.tensordot(zr_chunk, kz, axes=0)

            eikr = np.exp(1j * kr)

            integrand_x = np.zeros((len(xr_chunk), len(self.k_parallel),
                                    len(self.azimuthal_angles)),
                                   dtype=np.complex64)
            integrand_y = np.zeros((len(yr_chunk), len(self.k_parallel),
                                    len(self.azimuthal_angles)),
                                   dtype=np.complex64)
            integrand_z = np.zeros((len(zr_chunk), len(self.k_parallel),
                                    len(self.azimuthal_angles)),
                                   dtype=np.complex64)

            # pol=0
            integrand_x += (-kz * np.cos(agrid) * self.coefficients[0, :,
                                                  :])[None, :, :] * eikr
            integrand_y += (-kz * np.sin(agrid) * self.coefficients[0, :,
                                                  :])[None, :, :] * eikr
            integrand_z += (kpgrid * self.coefficients[0, :, :])[None, :,
                           :] * eikr
            # pol=1
            integrand_x += (- np.sin(agrid) * self.k * self.coefficients[1,
                                                       :, :])[None, :,
                           :] * eikr
            integrand_y += (np.cos(agrid) * self.k * self.coefficients[1, :,
                                                     :])[None, :, :] * eikr

            if len(self.k_parallel) > 1:
                h_x_flat[chunk_idcs] = 1 / omega * np.trapz(
                    np.trapz(integrand_x, self.azimuthal_angles)
                    * self.k_parallel, self.k_parallel)
                h_y_flat[chunk_idcs] = 1 / omega * np.trapz(
                    np.trapz(integrand_y, self.azimuthal_angles)
                    * self.k_parallel, self.k_parallel)
                h_z_flat[chunk_idcs] = 1 / omega * np.trapz(
                    np.trapz(integrand_z, self.azimuthal_angles)
                    * self.k_parallel, self.k_parallel)
            else:
                h_x_flat[chunk_idcs] = 1 / omega * np.squeeze(integrand_x)
                h_y_flat[chunk_idcs] = 1 / omega * np.squeeze(integrand_y)
                h_z_flat[chunk_idcs] = 1 / omega * np.squeeze(integrand_z)

        hx[self.valid(x, y, z)] = h_x_flat.reshape(xr.shape)
        hy[self.valid(x, y, z)] = h_y_flat.reshape(xr.shape)
        hz[self.valid(x, y, z)] = h_z_flat.reshape(xr.shape)

        return hx, hy, hz
        
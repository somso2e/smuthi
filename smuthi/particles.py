# -*- coding: utf-8 -*-
"""Classes for the representation of scattering particles."""
import smuthi.linearsystem.tmatrix.t_matrix as tmt
import smuthi.linearsystem.tmatrix.nfmds.indexconverter as nfic
import smuthi.linearsystem.tmatrix.nfmds.stlmanager as stlc
import smuthi.linearsystem.tmatrix.nfmds.nfmds as nfmds
import smuthi.utility.logging as log
import numpy as np
import tempfile
# import cylinder_t_matrix.t_matrix as t_alan   <- uncomment as soon as Alan's code is on PyPi


nfmds_logfile = None


class Particle:
    """Base class for scattering particles.

    Args:
        position (list):            Particle position in the format [x, y, z] (length unit)
        euler_angles (list):        Particle Euler angles in the format [alpha, beta, gamma]
        refractive_index (complex): Complex refractive index of particle
        l_max (int):                Maximal multipole degree used for the spherical wave expansion of incoming and
                                    scattered field
        m_max (int):                Maximal multipole order used for the spherical wave expansion of incoming and
                                    scattered field
    """

    def __init__(self, position=None, euler_angles=None, refractive_index=1 + 0j, l_max=None, m_max=None):

        if position is None:
            self.position = [0, 0, 0]
        else:
            self.position = position

        if euler_angles is None:
            self.euler_angles = [0, 0, 0]
        else:
            self.euler_angles = euler_angles

        self.refractive_index = refractive_index
        self.l_max = l_max
        if m_max is not None:
            self.m_max = m_max
        else:
            self.m_max = l_max
        self.initial_field = None
        self.scattered_field = None
        self.t_matrix = None

    def circumscribing_sphere_radius(self):
        """Virtual method to be overwritten"""
        pass

    def is_inside(self, x, y, z):
        """Virtual method to be overwritten.
        Until all child classes implement it: return False
        """
        return False

    def is_outside(self, x, y, z):
        """Virtual method to be overwritten.
        Until all child classes implement it: return True
        """
        return True

    def compute_t_matrix(self, vacuum_wavelength, n_medium):
        """Return the T-matrix of a particle.

        Args:
            vacuum_wavelength(float)
            n_medium(float or complex):             Refractive index of surrounding medium
            particle(smuthi.particles.Particle):    Particle object

        Returns:
            T-matrix as ndarray
        """
        raise ValueError('T-matrix for ' + type(self).__name__ + ' currently not implemented.')


class Sphere(Particle):
    """Particle subclass for spheres.

    Args:
        position (list):            Particle position in the format [x, y, z] (length unit)
        refractive_index (complex): Complex refractive index of particle
        radius (float):             Particle radius (length unit)
        l_max (int):                Maximal multipole degree used for the spherical wave expansion of incoming and
                                    scattered field
        m_max (int):                Maximal multipole order used for the spherical wave expansion of incoming and
                                    scattered field
    """

    def __init__(self, position=None, refractive_index=1 + 0j, radius=1, l_max=None, m_max=None):
        Particle.__init__(self, position=position, refractive_index=refractive_index, l_max=l_max, m_max=m_max)

        self.radius = radius

    def circumscribing_sphere_radius(self):
        return self.radius

    def is_inside(self, x, y, z):
        return (x - self.position[0]) ** 2 + (y - self.position[1]) ** 2 + (
                z - self.position[2]) ** 2 <= self.radius ** 2

    def is_outside(self, x, y, z):
        return (x - self.position[0]) ** 2 + (y - self.position[1]) ** 2 + (
                z - self.position[2]) ** 2 > self.radius ** 2

    def compute_t_matrix(self, vacuum_wavelength, n_medium):
        k_medium = 2 * np.pi / vacuum_wavelength * n_medium
        k_particle = 2 * np.pi / vacuum_wavelength * self.refractive_index
        radius = self.radius
        t = tmt.t_matrix_sphere(k_medium, k_particle, radius, self.l_max, self.m_max)
        return t


class AnisotropicSphere(Particle):
    """Particle subclass for anisotropic spheres.

    Args:
        position (list):              Particle position in the format [x, y, z] (length unit)
        euler_angles (list):          Euler angles [alpha, beta, gamma] in (zy'z''-convention) in radian.
                                      Alternatively, you can specify the polar and azimuthal angle of the axis of
                                      revolution.
        polar_angle (float):          Polar angle of axis of revolution.
        azimuthal_angle (float):      Azimuthal angle of axis of revolution.
        refractive_index (complex):   Complex refractive index of particle in x-y plane (if not rotated)
        refractive_index_z (complex): Complex refractive index of particle along z-axis (if not rotated)
        radius (float):               Sphere radius
        l_max (int):                  Maximal multipole degree used for the spherical wave expansion of incoming and
                                      scattered field
        m_max (int):                  Maximal multipole order used for the spherical wave expansion of incoming and
                                      scattered field
        n_rank (int):                 Maximal multipole order used for in NFMDS (default: l_max + 5)
    """

    def __init__(self, position=None, euler_angles=None, polar_angle=0, azimuthal_angle=0, refractive_index=1 + 0j,
                 radius=1, refractive_index_z=2 + 0j, l_max=None, m_max=None, n_rank=None):

        self.radius = radius
        self.refractive_index_z = refractive_index_z
        self.nrank = n_rank

        if euler_angles is None:
            euler_angles = [azimuthal_angle, polar_angle, 0]

        Particle.__init__(self, position=position, euler_angles=euler_angles, refractive_index=refractive_index,
                          l_max=l_max, m_max=m_max)

    def circumscribing_sphere_radius(self):
        return self.radius

    def compute_t_matrix(self, vacuum_wavelength, n_medium):
        return self._compute_t_matrix_nfmds(vacuum_wavelength, n_medium)

    def _compute_t_matrix_nfmds(self, vacuum_wavelength, n_medium):
        nrank = self.nrank if self.nrank is not None else self.l_max + 5
        Nmax = nrank * (2 + nrank)
        r = self.radius
        surf = np.array([r,r,r])
        
        with log.LoggerLowLevelMuted(filename=nfmds_logfile):
            tnfmds = nfmds.tnonaxsym(surf, Nmax, filegeom=0,
                                     wavelength=vacuum_wavelength, ind_refrel=self.refractive_index / n_medium + 0j,
                                     ind_refrelz=self.refractive_index_z / n_medium + 0j,
                                     nrank=nrank, mrank=nrank, ind_refmed=n_medium,
                                     anisotropic=1, typegeom=1, nparam=1)
        t = nfic.nfmds_to_smuthi_matrix(tnfmds, l_max=self.l_max, m_max=self.m_max)
        return t


class Spheroid(Particle):
    """Particle subclass for spheroids.

    Args:
        position (list):            Particle position in the format [x, y, z] (length unit)
        euler_angles (list):        Euler angles [alpha, beta, gamma] in (zy'z''-convention) in radian.
                                    Alternatively, you can specify the polar and azimuthal angle of the axis of 
                                    revolution.
        polar_angle (float):        Polar angle of axis of revolution. 
        azimuthal_angle (float):    Azimuthal angle of axis of revolution.
        refractive_index (complex): Complex refractive index of particle
        semi_axis_c (float):        Spheroid half axis in direction of axis of revolution (z-axis if not rotated)
        semi_axis_a (float):        Spheroid half axis in lateral direction (x- and y-axis if not rotated)
        l_max (int):                Maximal multipole degree used for the spherical wave expansion of incoming and
                                    scattered field
        m_max (int):                Maximal multipole order used for the spherical wave expansion of incoming and
                                    scattered field
        n_rank (int):               Maximal multipole order used in NFMDS (default: l_max + 5)
    """

    def __init__(self, position=None, euler_angles=None, polar_angle=0, azimuthal_angle=0, refractive_index=1 + 0j,
                 semi_axis_c=1, semi_axis_a=1, l_max=None, m_max=None, n_rank=None):

        self.semi_axis_c = semi_axis_c
        self.semi_axis_a = semi_axis_a
        self.nrank = n_rank
        if euler_angles is None:
            euler_angles = [azimuthal_angle, polar_angle, 0]

        Particle.__init__(self, position=position, euler_angles=euler_angles, refractive_index=refractive_index,
                          l_max=l_max, m_max=m_max)

    def circumscribing_sphere_radius(self):
        return max([self.semi_axis_a, self.semi_axis_c])

    def compute_t_matrix(self, vacuum_wavelength, n_medium):
        return self._compute_t_matrix_nfmds(vacuum_wavelength, n_medium)

    def _compute_t_matrix_nfmds(self, vacuum_wavelength, n_medium):
        nrank = self.nrank if self.nrank is not None else self.l_max + 5
        Nmax = nrank * (2 + nrank)
        surf = [self.semi_axis_c, self.semi_axis_a]
        with log.LoggerLowLevelMuted(filename=nfmds_logfile):
            tnfmds = nfmds.taxsym(surf, Nmax, typegeom=1, nparam=1,
                                  wavelength=vacuum_wavelength, ind_refrel=self.refractive_index / n_medium + 0j,
                                  nrank=nrank, ind_refmed=n_medium)
        t = nfic.nfmds_to_smuthi_matrix(tnfmds, l_max=self.l_max, m_max=self.m_max)
        return t


class LayeredSpheroid(Particle):
    """Particle subclass for layered spheroid.

    Args:
        position (list):            Particle position in the format [x, y, z] (length unit)
        euler_angles (list):        Euler angles [alpha, beta, gamma] in (zy'z''-convention) in radian.
                                    Alternatively, you can specify the polar and azimuthal angle of the axis of
                                    revolution.
        polar_angle (float):        Polar angle of axis of revolution.
        azimuthal_angle (float):    Azimuthal angle of axis of revolution.
        layer_refractive_indices (complex): Complex refractive index of particle
        layer_semi_axes_c (float):    Spheroid half axis in direction of axis of revolution (z-axis if not rotated)
        layer_semi_axes_a (float):    Spheroid half axis in lateral direction (x- and y-axis if not rotated)
        l_max (int):                Maximal multipole degree used for the spherical wave expansion of incoming and
                                    scattered field
        m_max (int):                Maximal multipole order used for the spherical wave expansion of incoming and
                                    scattered field
        n_rank (int):               Maximal multipole order used in NFMDS (default: l_max + 5)
    """

    def __init__(self, position=None, euler_angles=None, polar_angle=0, azimuthal_angle=0, layer_refractive_indices=1 + 0j,
                 layer_semi_axes_c=1, layer_semi_axes_a=1, l_max=None, m_max=None, n_rank=None):
        if euler_angles is None:
            euler_angles = [azimuthal_angle, polar_angle, 0]
        Particle.__init__(self, position=position, euler_angles=euler_angles, refractive_index=None, l_max=l_max, m_max=m_max)
        self.n_rank = n_rank
        self.layer_refractive_indices = layer_refractive_indices
        self.layer_semi_axes_c = layer_semi_axes_c
        self.layer_semi_axes_a = layer_semi_axes_a

    def compute_t_matrix(self, vacuum_wavelength, n_medium):
        nrank = self.nrank if self.nrank is not None else self.l_max + 5
        surf = np.hstack((self.layer_semi_axes_c[:,np.newaxis],self.layer_semi_axes_a[:,np.newaxis]))
        k = 2 * np.pi / vacuum_wavelength
        with log.LoggerLowLevelMuted(filename=nfmds_logfile):
            tnfmds = nfmds.tlay(k,self.layer_refractive_indices,surf,nrank)
        t = nfic.nfmds_to_smuthi_matrix(tnfmds, l_max=self.l_max, m_max=self.m_max)
        return t


class FiniteCylinder(Particle):
    """Particle subclass for finite cylinders.

    Args:
        position (list):            Particle position in the format [x, y, z] (length unit)
        euler_angles (list):        Euler angles [alpha, beta, gamma] in (zy'z''-convention) in radian.
                                    Alternatively, you can specify the polar and azimuthal angle of the axis of 
                                    revolution.
        polar_angle (float):        Polar angle of axis of revolution. 
        azimuthal_angle (float):    Azimuthal angle of axis of revolution.
        refractive_index (complex): Complex refractive index of particle
        cylinder_radius (float):    Radius of cylinder (length unit)
        cylinder_height (float):    Height of cylinder, in z-direction if not rotated (length unit)
        l_max (int):                Maximal multipole degree used for the spherical wave expansion of incoming and
                                    scattered field
        m_max (int):                Maximal multipole order used for the spherical wave expansion of incoming and
                                    scattered field
        n_rank (int):               Maximal multipole order used for in NFMDS (default: l_max + 5)
        use_python_tmatrix (bool):  If true, use Alan Zhan's Python code to compute the T-matrix rather than NFM-DS
        nint (int):                 Number of angles used in integral (only for python t-mnatrix)
    """

    def __init__(self, position=None, euler_angles=None, polar_angle=0, azimuthal_angle=0, refractive_index=1 + 0j,
                 cylinder_radius=1, cylinder_height=1, l_max=None, m_max=None, n_rank=None, use_python_tmatrix=False,
                 nint=100):

        self.cylinder_radius = cylinder_radius
        self.cylinder_height = cylinder_height
        self.nrank = n_rank
        self.python_tmt = use_python_tmatrix
        self.nint = nint
        if euler_angles is None:
            euler_angles = [azimuthal_angle, polar_angle, 0]

        Particle.__init__(self, position=position, euler_angles=euler_angles, refractive_index=refractive_index,
                          l_max=l_max, m_max=m_max)

    def circumscribing_sphere_radius(self):
        return np.sqrt((self.cylinder_height / 2) ** 2 + self.cylinder_radius ** 2)

    def compute_t_matrix(self, vacuum_wavelength, n_medium):
        if self.python_tmt:
            return self._compute_t_matrix_python(vacuum_wavelength, n_medium)
        else:
            return self._compute_t_matrix_nfmds(vacuum_wavelength, n_medium)

    def _compute_t_matrix_nfmds(self, vacuum_wavelength, n_medium):
        nrank = self.nrank if self.nrank is not None else self.l_max + 5
        Nmax = nrank * (2 + nrank)
        surf = [self.cylinder_height / 2, self.cylinder_radius]
        with log.LoggerLowLevelMuted(filename=nfmds_logfile):
            tnfmds = nfmds.taxsym(surf, Nmax, typegeom=2, nparam=3, wavelength=vacuum_wavelength,
                                  ind_refrel=self.refractive_index / n_medium + 0j, nrank=nrank, ind_refmed=n_medium)
        t = nfic.nfmds_to_smuthi_matrix(tnfmds, l_max=self.l_max, m_max=self.m_max)
        return t

    def _compute_t_matrix_python(self, vacuum_wavelength, n_medium):
        nrank = self.nrank if self.nrank is not None else self.l_max + 5
        Ntheta = self.nint
        geometric_params = [self.cylinder_radius, self.cylinder_height]
        t, dt = t_alan.compute_T(lmax=nrank, Ntheta=Ntheta, geometric_params=geometric_params,
                                 n0=n_medium, ns=self.refractive_index, wavelength=vacuum_wavelength,
                                 particle_type="cylinder")
        t = -nfic.python_to_smuthi_matrix(t, Nrank=nrank, l_max=self.l_max, m_max=self.m_max)
        return t


class CustomParticle(Particle):
    """Particle subclass for custom particle shapes defined via FEM file.

    Args:
        position (list):            Particle position in the format [x, y, z] (length unit)
        euler_angles (list):        Euler angles [alpha, beta, gamma] in (zy'z''-convention) in radian.
                                    Alternatively, you can specify the polar and azimuthal angle of the axis of
                                    revolution.
        polar_angle (float):        Polar angle of axis of revolution.
        azimuthal_angle (float):    Azimuthal angle of axis of revolution.
        geometry_filename (string):      Path to FEM file
        scale (float):              Scaling factor for particle dimensions (relative to provided dimensions)
        l_max (int):                Maximal multipole degree used for the spherical wave expansion of incoming and
                                    scattered field
        m_max (int):                Maximal multipole order used for the spherical wave expansion of incoming and
                                    scattered field
        n_rank (int):               Maximal multipole order used for in NFMDS (default: l_max + 5)
    """

    def __init__(self, position=None, euler_angles=None, polar_angle=0, azimuthal_angle=0, refractive_index=1 + 0j,
                 geometry_filename=None, scale=1, l_max=None, m_max=None, n_rank=None):
        if euler_angles is None:
            euler_angles = [azimuthal_angle, polar_angle, 0]
        Particle.__init__(self, position=position, euler_angles=euler_angles, refractive_index=refractive_index,
                l_max=l_max, m_max=m_max)
        self.nrank = n_rank
        self.geometry_filename = geometry_filename
        self.scale = scale

    def circumscribing_sphere_radius(self):
        return self.scale

    def compute_t_matrix(self, vacuum_wavelength, n_medium):
        if self.geometry_filename.endswith(".fem") or self.geometry_filename.endswith(".FEM"):
            return self._compute_t_matrix_nfmds(vacuum_wavelength, n_medium, self.geometry_filename)

        elif self.geometry_filename.endswith(".stl") or self.geometry_filename.endswith(".STL"):
            with tempfile.TemporaryDirectory() as tempdir:
                stlc.convert_stl_to_fem(stlname=self.geometry_filename,
                                        femname=tempdir + "/temp.fem")
                return self._compute_t_matrix_nfmds(vacuum_wavelength, n_medium, tempdir + "/temp.fem")
        else:
            raise Exception("Invalid geometry file extension.")

    def _compute_t_matrix_nfmds(self, vacuum_wavelength, n_medium, fem_file):
        """Private t-matrix method function"""
        nrank = self.nrank if self.nrank is not None else self.l_max + 5
        Nmax = nrank * (2 + nrank)
        surf = np.array([1, 1, 1])	
        with log.LoggerLowLevelMuted(filename=nfmds_logfile):
            tnfmds = nfmds.tnonaxsym(surf, Nmax, filefem=fem_file,
                                     wavelength=vacuum_wavelength / self.scale,
                                     ind_refrel=self.refractive_index / n_medium + 0j,
                                     nrank=nrank, mrank=nrank, ind_refmed=n_medium)
        t = nfic.nfmds_to_smuthi_matrix(tnfmds, l_max=self.l_max, m_max=self.m_max)
        return t

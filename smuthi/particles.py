# -*- coding: utf-8 -*-
"""Provide class for the representation of scattering particles."""
import smuthi.fields as flds
import smuthi.linearsystem.tmatrix.t_matrix as tmt
import numpy as np
import smuthi.linearsystem.tmatrix.nfmds.indexconverter as nfic
import smuthi.linearsystem.tmatrix.nfmds.nfmds as nfmds


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
    def __init__(self, position=None, euler_angles=None, refractive_index=1+0j, l_max=None, m_max=None):
        
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
        t_matrix_method (dict):     Dictionary containing the parameters for the algorithm to compute the T-matrix
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


class CustomParticle(Particle):
    """Particle subclass for custom particle shapes defined via FEM file.

    Args:
        position (list):            Particle position in the format [x, y, z] (length unit)
        euler_angles (list):        Euler angles [alpha, beta, gamma] in (zy'z''-convention) in radian.
                                    Alternatively, you can specify the polar and azimuthal angle of the axis of 
                                    revolution.
        polar_angle (float):        Polar angle of axis of revolution. 
        azimuthal_angle (float):    Azimuthal angle of axis of revolution.
        fem_filename (string):      Path to FEM file
        scale (float):              Scaling factor for particle dimensions (relative to provided dimensions)
        l_max (int):                Maximal multipole degree used for the spherical wave expansion of incoming and
                                    scattered field
        m_max (int):                Maximal multipole order used for the spherical wave expansion of incoming and
                                    scattered field
    """

    def __init__(self, position=None, euler_angles=None, polar_angle=0, azimuthal_angle=0, refractive_index=1 + 0j,
                 fem_filename=None, scale=1, l_max=None, m_max=None, t_matrix_method=None):
        if euler_angles is None:
            euler_angles = [azimuthal_angle, polar_angle, 0]
        if t_matrix_method is None:
            self.t_matrix_method = {}
        else:
            self.t_matrix_method = t_matrix_method
        Particle.__init__(self, position=position, euler_angles=euler_angles, refractive_index=refractive_index,
                          l_max=l_max, m_max=m_max)

        self.fem_filename = fem_filename
        self.scale = scale

    def circumscribing_sphere_radius(self):
        return self.scale

    def compute_t_matrix(self, vacuum_wavelength, n_medium):
        Nrank=self.t_matrix_method.get('nrank', self.l_max + 2)
        Mrank=self.t_matrix_method.get('mrank', self.l_max + 2)
        Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
        tnfmds = nfmds.tnonaxsym([1, 1, 1, 0, 0, 0, 0, 0, 0, 0], Nmax, filefem=self.fem_filename,
                                 wavelength=vacuum_wavelength / self.scale, ind_refrel=self.refractive_index/n_medium + 0j,
                                 nrank=Nrank, mrank=Mrank,ind_refmed=n_medium)
        t = nfic.nfmds_to_smuthi_matrix(tnfmds,l_max=self.l_max)
        return t


class AxisymmetricParticle(Particle):
    """Particle subclass for axisymmetric particles calculated with NFMDS.

    Args:
        position (list):            Particle position in the format [x, y, z] (length unit)
        euler_angles (list):        Euler angles [alpha, beta, gamma] in (zy'z''-convention) in radian.
                                    Alternatively, you can specify the polar and azimuthal angle of the axis of
                                    revolution.
        polar_angle (float):        Polar angle of axis of revolution.
        azimuthal_angle (float):    Azimuthal angle of axis of revolution.
        refractive_index (complex): Complex refractive index of particle
        geometry_type (int):        Particle shape in NFMDS convention (see table below)
        geometry_parameters (list): List of parameters specific to a given shape (see table below)
        l_max (int):                Maximal multipole degree used for the spherical wave expansion of incoming and
                                    scattered field
        m_max (int):                Maximal multipole order used for the spherical wave expansion of incoming and
                                    scattered field
!    Particle      TypeGeom   Nsurf   Nparam                surf                   !
!    spheroid         1         2       1         surf(1) - length of the semi-    !
!                                                           axis along the         !
!                                                           symmetry axis          !
!                                                 surf(2) - length of the second   !
!                                                           semi-axis              !
!    cylinder         2         2       3         surf(1) - half-length of         !
!                                                           the cylinder           !
!                                                 surf(2) - cylinder radius        !
!                                                                                  !
!    rounded          3         2       3         surf(1) - half-length of         !
!     oblate                                                the cylinder           ! 
!    cylinder                                     surf(2) - cylinder radius        !
!                                                           including the rounded  !
!                                                           part                   ! 
    """

    def __init__(self, position=None, euler_angles=None, polar_angle=0, azimuthal_angle=0, refractive_index=1 + 0j,
                 geometry_type = None, geometry_parameters = None, l_max=None, m_max=None, t_matrix_method=None):
        if euler_angles is None:
            euler_angles = [azimuthal_angle, polar_angle, 0]

        Particle.__init__(self, position=position, euler_angles=euler_angles, refractive_index=refractive_index,
                          l_max=l_max, m_max=m_max)

        if t_matrix_method is None:
            self.t_matrix_method={}

        self.geometry_type = geometry_type
        self.geometry_parameters = geometry_parameters
        if geometry_type == 1:
            self.nparam = 1
        else:
            self.nparam = 3
    def circumscribing_sphere_radius(self):
        return self.geometry_parameters[0]

    def compute_t_matrix(self, vacuum_wavelength, n_medium):
        Nrank=self.t_matrix_method.get('nrank', self.l_max + 2)
        Mrank=self.t_matrix_method.get('mrank', self.l_max + 2)
        Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
        nsurf=len(self.geometry_parameters)
        surf=list(self.geometry_parameters)+[0]*(10-nsurf)
        tnfmds = nfmds.taxsym(surf, Nmax, typegeom=self.geometry_type, nsurf=nsurf, nparam=self.nparam,
                             wavelength=vacuum_wavelength, ind_refrel=self.refractive_index/n_medium + 0j,
                             nrank=Nrank,ind_refmed=n_medium)
        t = nfic.nfmds_to_smuthi_matrix(tnfmds,l_max=self.l_max)
        return t

class Spheroid(AxisymmetricParticle):
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
        t_matrix_method (dict):     Dictionary containing the parameters for the algorithm to compute the T-matrix
    """

    def __init__(self, position=None, euler_angles=None, polar_angle=0, azimuthal_angle=0, refractive_index=1 + 0j,
                 semi_axis_c=1, semi_axis_a=1, l_max=None, m_max=None, t_matrix_method=None):

        if euler_angles is None:
            euler_angles = [azimuthal_angle, polar_angle, 0]

        if t_matrix_method is None:
            self.t_matrix_method = {}
        else:
            self.t_matrix_method = t_matrix_method

        self.semi_axis_c = semi_axis_c
        self.semi_axis_a = semi_axis_a
        AxisymmetricParticle.__init__(self, position=position, euler_angles=euler_angles, refractive_index=refractive_index,geometry_type = 1, geometry_parameters = [self.semi_axis_c, self.semi_axis_a], l_max=l_max, m_max=m_max, t_matrix_method=t_matrix_method)

    def circumscribing_sphere_radius(self):
        return max([self.semi_axis_a, self.semi_axis_c])


class FiniteCylinder(AxisymmetricParticle):
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
    """

    def __init__(self, position=None, euler_angles=None, polar_angle=0, azimuthal_angle=0, refractive_index=1 + 0j,
                 cylinder_radius=1, cylinder_height=1, l_max=None, m_max=None, t_matrix_method=None):

        if euler_angles is None:
            euler_angles = [azimuthal_angle, polar_angle, 0]

        if t_matrix_method is None:
            self.t_matrix_method = {}
        else:
            self.t_matrix_method = t_matrix_method

        self.cylinder_radius = cylinder_radius
        self.cylinder_height = cylinder_height

        AxisymmetricParticle.__init__(self, position=position, euler_angles=euler_angles, refractive_index=refractive_index,geometry_type = 2, geometry_parameters = [self.cylinder_height / 2, self.cylinder_radius], l_max=l_max, m_max=m_max, t_matrix_method=t_matrix_method)

    def circumscribing_sphere_radius(self):
        return np.sqrt((self.cylinder_height / 2) ** 2 + self.cylinder_radius ** 2)



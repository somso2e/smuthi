"""Manage post processing steps to evaluate the scattered far field"""

import os
import numpy as np
import smuthi.fields as flds
import smuthi.postprocessing.scattered_field as sf

class FarField:
    r"""Represent the far field intensity of an electromagnetic field.

    .. math::
        P = \sum_{j=1}^2 \iint \mathrm{d}^2 \Omega \, I_{\Omega,j}(\beta, \alpha),

    where :math:`P` is the radiative power, :math:`j` indicates the polarization and
    :math:`\mathrm{d}^2 \Omega = \mathrm{d}\alpha \sin\beta \mathrm{d}\beta` denotes the infinitesimal solid angle.

    Args:
        polar_angles (numpy.ndarray):       array of polar angles for plane wave expansions. If 'default', use
                                            smuthi.fields.default_polar_angles
        azimuthal_angles (ndarray or str):  array of azimuthal angles for plane wave expansions. If 'default', use
                                            smuthi.fields.default_azimuthal_angles
        angular_resolution (float):         If provided, angular arrays are generated with this angular resolution
                                            over the default angular range
        signal_type (str):                  Type of the signal (e.g., 'intensity' for power flux far fields).
    """

    def __init__(self, polar_angles='default', azimuthal_angles='default', angular_resolution=None, signal_type='intensity'):
        azimuthal_angles, polar_angles = _angular_arrays(azimuthal_angles, polar_angles, angular_resolution)
        self.polar_angles = polar_angles
        self.azimuthal_angles = azimuthal_angles

        # The far field signal is represented as a 3-dimensional numpy.ndarray.
        # The indices are:
        # - polarization (0=TE, 1=TM)
        # - index of the polar angle
        # - index of the azimuthal angle
        self.signal = np.zeros((2, len(polar_angles), len(azimuthal_angles)))
        self.signal.fill(np.nan)
        self.signal_type = signal_type

    def azimuthal_integral_times_sin_beta(self):
        r"""Far field as a function of polar angle only.

        .. math::
            P = \sum_{j=1}^2 \int \mathrm{d} \beta \, I_{\beta,j}(\beta),

        with

        .. math::
            I_{\beta,j}(\beta) = \int \mathrm{d} \alpha \, \sin\beta I_j(\beta, \alpha),

        Returns:
            :math:`I_{\beta,j}(\beta)` as numpy ndarray. First index is polarization, second is polar angle.
        """
        if len(self.azimuthal_angles) > 2:
            return np.trapz(self.signal, self.azimuthal_angles[None, None, :]) * np.sin(self.polar_angles[None, :])
        else:
            return None

    def azimuthal_integral(self):
        r"""Far field as a function of the polar angle cosine only.

        .. math::
            P = \sum_{j=1}^2 \int \mathrm{d} \cos\beta \, I_{\cos\beta,j}(\beta),

        with

        .. math::
            I_{\beta,j}(\beta) = \int \mathrm{d} \alpha \, I_j(\beta, \alpha),

        Returns:
            :math:`I_{\cos\beta,j}(\beta)` as numpy ndarray. First index is polarization, second is polar angle.
        """
        if len(self.azimuthal_angles) > 2:
            return np.trapz(self.signal, self.azimuthal_angles[None, None, :])
        else:
            return None

    def integral(self):
        r"""Integrate intensity to obtain total power :math:`P`.

        Returns:
            :math:`P_j` as numpy 1D-array with length 2, the index referring to polarization.
        """
        if len(self.azimuthal_angles) > 2:
            return np.trapz(self.azimuthal_integral_times_sin_beta(), self.polar_angles[None, :])
        else:
            return None

    def top(self):
        r"""Split far field into top and bottom part.

        Returns:
            FarField object with only the intensity for top hemisphere (:math:`\beta\leq\pi/2`)
        """
        if any(self.polar_angles <= np.pi / 2):
            ff = FarField(polar_angles=self.polar_angles[self.polar_angles <= np.pi / 2],
                          azimuthal_angles=self.azimuthal_angles, signal_type=self.signal_type)
            ff.signal = self.signal[:, self.polar_angles <= np.pi / 2, :]
            return ff
        else:
            return None

    def bottom(self):
        r"""Split far field into top and bottom part.

        Returns:
            FarField object with only the intensity for bottom hemisphere (:math:`\beta\geq\pi/2`)
        """
        if any(self.polar_angles >= np.pi / 2):
            ff = FarField(polar_angles=self.polar_angles[self.polar_angles >= np.pi / 2],
                          azimuthal_angles=self.azimuthal_angles, signal_type=self.signal_type)
            ff.signal = self.signal[:, self.polar_angles >= np.pi / 2, :]
            return ff
        else:
            return None

    def alpha_grid(self):
        r"""
        Returns:
            Meshgrid with :math:`\alpha` values.
        """
        agrid, _ = np.meshgrid(self.azimuthal_angles, self.polar_angles.real)
        return agrid

    def beta_grid(self):
        r"""
        Returns:
            Meshgrid with :math:`\beta` values.
        """
        _, bgrid = np.meshgrid(self.azimuthal_angles, self.polar_angles.real)
        return bgrid

    def append(self, other):
        """Combine two FarField objects with disjoint angular ranges. The other far field is appended to this one.

        Args:
            other (FarField): far field to append to this one.
        """

        if not all(self.azimuthal_angles == other.azimuthal_angles):
            raise ValueError('azimuthal angles not consistent')
        if not self.signal_type == other.signal_type:
            raise ValueError('signal type not consistent')
        if max(self.polar_angles) <= min(other.polar_angles):
            self.polar_angles = np.concatenate((self.polar_angles, other.polar_angles))
            self.signal = np.concatenate((self.signal, other.signal), 1)
        elif min(self.polar_angles) >= max(other.polar_angles):
            self.polar_angles = np.concatenate((other.polar_angles, self.polar_angles))
            self.signal = np.concatenate((other.signal, self.signal), 1)
        else:
            raise ValueError('far fields have overlapping polar angle domains')


def _angular_arrays(azimuthal_angles='default', polar_angles='default', angular_resolution=None):
    """ Utility function returning angular arrays based on user input """ # NOTE: if both arrays and a resolution are provided, the arrays are ignored!
    if angular_resolution is not None:
        azimuthal_angles, polar_angles = flds.angular_arrays(angular_resolution)
    else:
        if type(polar_angles) == str and polar_angles == 'default':
            polar_angles = flds.default_polar_angles
        if type(azimuthal_angles) == str and azimuthal_angles == 'default':
            azimuthal_angles = flds.default_azimuthal_angles
    return azimuthal_angles, polar_angles


def pwe_to_ff_conversion(vacuum_wavelength, plane_wave_expansion):
    """Compute the far field of a plane wave expansion object.

    Args:
        vacuum_wavelength (float):                 Vacuum wavelength in length units.
        plane_wave_expansion (PlaneWaveExpansion): Plane wave expansion to convert into far field object.

    Returns:
        A FarField object containing the far field intensity.
    """
    omega = flds.angular_frequency(vacuum_wavelength)
    k = plane_wave_expansion.k
    kp = plane_wave_expansion.k_parallel
    if plane_wave_expansion.kind == 'upgoing':
        polar_angles = np.arcsin(kp / k)
    elif plane_wave_expansion.kind == 'downgoing':
        polar_angles = np.pi - np.arcsin(kp / k)
    else:
        raise ValueError('PWE type not specified')
    if any(polar_angles.imag):
        raise ValueError('complex angles are not allowed')
    azimuthal_angles = plane_wave_expansion.azimuthal_angles
    kkz2 = flds.k_z(k_parallel=kp, k=k) ** 2 * k
    intens = (2 * np.pi ** 2 / omega * kkz2[np.newaxis, :, np.newaxis]
              * abs(plane_wave_expansion.coefficients) ** 2).real
    srt_idcs = np.argsort(polar_angles)  # reversing order in case of downgoing
    ff = FarField(polar_angles=polar_angles[srt_idcs], azimuthal_angles=azimuthal_angles)
    ff.signal = intens[:, srt_idcs, :]
    return ff


def total_far_field(initial_field, particle_list, layer_system,
                    polar_angles='default', azimuthal_angles='default', angular_resolution=None):
    """
    Evaluate the total far field, the initial far field and the scattered far field. Cannot be used if initial field
    is a plane wave.

    Args:
        initial_field (smuthi.initial_field.InitialField): represents the initial field
        particle_list (list):                       list of smuthi.Particle objects
        layer_system (smuthi.layers.LayerSystem):   represents the stratified medium
        polar_angles (numpy.ndarray or str):        polar angles values (radian).
                                                    if 'default', use smuthi.fields.default_polar_angles
        azimuthal_angles (numpy.ndarray or str):    azimuthal angle values (radian)
                                                    if 'default', use smuthi.fields.default_azimuthal_angles
        angular_resolution (float):                 If provided, angular arrays are generated with this angular
                                                    resolution over the default angular range

    Returns:
        A tuple of three smuthi.field_expansion.FarField objects for total, initial and scattered far field. Mind that the scattered far field
        has no physical meaning and is for illustration purposes only.
    """
    if not (type(initial_field).__name__ == 'GaussianBeam' or type(initial_field).__name__ == 'DipoleSource'
            or type(initial_field).__name__ == 'DipoleCollection'):
        raise ValueError('only for Gaussian beams and dipole sources')
    omega = initial_field.angular_frequency()
    vacuum_wavelength = initial_field.vacuum_wavelength
    azimuthal_angles, polar_angles = _angular_arrays(azimuthal_angles, polar_angles, angular_resolution)

    if any(polar_angles.imag):
        raise ValueError("complex angles not allowed in far field")

    i_top = layer_system.number_of_layers() - 1
    top_polar_angles = polar_angles[polar_angles < (np.pi / 2)]
    bottom_polar_angles = polar_angles[polar_angles > (np.pi / 2)]
    neff_top = np.sort(np.sin(top_polar_angles) * layer_system.refractive_indices[i_top])
    neff_bottom = np.sort(np.sin(bottom_polar_angles) * layer_system.refractive_indices[0])

    if len(top_polar_angles) > 1 and layer_system.refractive_indices[i_top].imag == 0:
        pwe_scat_top, _ = sf.scattered_field_pwe(vacuum_wavelength, particle_list, layer_system, i_top,
                                                 k_parallel=neff_top*omega, azimuthal_angles=azimuthal_angles,
                                                 include_direct=True, include_layer_response=True)
        pwe_in_top, _ = initial_field.plane_wave_expansion(layer_system, i_top, k_parallel_array=neff_top*omega,
                                                           azimuthal_angles_array=azimuthal_angles)
        pwe_top = pwe_scat_top + pwe_in_top
        top_far_field = pwe_to_ff_conversion(vacuum_wavelength=vacuum_wavelength, plane_wave_expansion=pwe_top)
        top_far_field_init = pwe_to_ff_conversion(vacuum_wavelength=vacuum_wavelength,
                                                  plane_wave_expansion=pwe_in_top)
        top_far_field_scat = pwe_to_ff_conversion(vacuum_wavelength=vacuum_wavelength,
                                                  plane_wave_expansion=pwe_scat_top)
    else:
        top_far_field = None
        top_far_field_init = None
        top_far_field_scat = None

    if len(bottom_polar_angles) > 1 and layer_system.refractive_indices[0].imag == 0:
        _, pwe_scat_bottom = sf.scattered_field_pwe(vacuum_wavelength, particle_list, layer_system, 0,
                                                   k_parallel=neff_bottom*omega, azimuthal_angles=azimuthal_angles,
                                                   include_direct=True, include_layer_response=True)
        _, pwe_in_bottom = initial_field.plane_wave_expansion(layer_system, 0, k_parallel_array=neff_bottom*omega,
                                                              azimuthal_angles_array=azimuthal_angles)
        pwe_bottom = pwe_scat_bottom + pwe_in_bottom
        bottom_far_field = pwe_to_ff_conversion(vacuum_wavelength=vacuum_wavelength,
                                                plane_wave_expansion=pwe_bottom)
        bottom_far_field_init = pwe_to_ff_conversion(vacuum_wavelength=vacuum_wavelength,
                                                     plane_wave_expansion=pwe_in_bottom)
        bottom_far_field_scat = pwe_to_ff_conversion(vacuum_wavelength=vacuum_wavelength,
                                                     plane_wave_expansion=pwe_scat_bottom)
    else:
        bottom_far_field = None
        bottom_far_field_init = None
        bottom_far_field_scat = None

    if top_far_field is not None:
        far_field = top_far_field
        far_field_init = top_far_field_init
        far_field_scat = top_far_field_scat
        if bottom_far_field is not None:
            far_field.append(bottom_far_field)
            far_field_init.append(bottom_far_field_init)
            far_field_scat.append(bottom_far_field_scat)
    else:
        far_field = bottom_far_field
        far_field_init = bottom_far_field_init
        far_field_scat = bottom_far_field_scat

    far_field.polar_angles = far_field.polar_angles.real
    far_field_init.polar_angles = far_field_init.polar_angles.real
    far_field_scat.polar_angles = far_field_scat.polar_angles.real
    return far_field, far_field_init, far_field_scat


def scattered_far_field(vacuum_wavelength, particle_list, layer_system,
                        polar_angles='default', azimuthal_angles='default', angular_resolution=None):
    """
    Evaluate the scattered far field.

    Args:
        vacuum_wavelength (float):                  in length units
        particle_list (list):                       list of smuthi.Particle objects
        layer_system (smuthi.layers.LayerSystem):   represents the stratified medium
        polar_angles (numpy.ndarray or str):        polar angles values (radian).
                                                    if 'default', use smuthi.fields.default_polar_angles
        azimuthal_angles (numpy.ndarray or str):    azimuthal angle values (radian)
                                                    if 'default', use smuthi.fields.default_azimuthal_angles
        angular_resolution (float):                 If provided, angular arrays are generated with this angular
                                                    resolution over the default angular range

    Returns:
        A smuthi.field_expansion.FarField object of the scattered field.
    """
    omega = flds.angular_frequency(vacuum_wavelength)
    azimuthal_angles, polar_angles = _angular_arrays(azimuthal_angles, polar_angles, angular_resolution)

    if any(polar_angles.imag):
        raise ValueError("complex angles not allowed in far field")

    i_top = layer_system.number_of_layers() - 1
    top_polar_angles = polar_angles[polar_angles < (np.pi / 2)]
    bottom_polar_angles = polar_angles[polar_angles > (np.pi / 2)]
    neff_top = np.sort(np.sin(top_polar_angles) * layer_system.refractive_indices[i_top])
    neff_bottom = np.sort(np.sin(bottom_polar_angles) * layer_system.refractive_indices[0])

    if len(top_polar_angles) > 1 and layer_system.refractive_indices[i_top].imag == 0:
        pwe_top, _ = sf.scattered_field_pwe(vacuum_wavelength, particle_list, layer_system, i_top,
                                            k_parallel=neff_top*omega, azimuthal_angles=azimuthal_angles,
                                            include_direct=True, include_layer_response=True)
        top_far_field = pwe_to_ff_conversion(vacuum_wavelength=vacuum_wavelength, plane_wave_expansion=pwe_top)
    else:
        top_far_field = None

    if len(bottom_polar_angles) > 1 and layer_system.refractive_indices[0].imag == 0:
        _, pwe_bottom = sf.scattered_field_pwe(vacuum_wavelength, particle_list, layer_system, 0,
                                               k_parallel=neff_bottom*omega, azimuthal_angles=azimuthal_angles,
                                               include_direct=True, include_layer_response=True)
        bottom_far_field = pwe_to_ff_conversion(vacuum_wavelength=vacuum_wavelength,
                                                plane_wave_expansion=pwe_bottom)
    else:
        bottom_far_field = None

    if top_far_field is not None:
        far_field = top_far_field
        if bottom_far_field is not None:
            far_field.append(bottom_far_field)
    else:
        far_field = bottom_far_field

    far_field.polar_angles = far_field.polar_angles.real
    return far_field


def scattering_cross_section(initial_field, particle_list, layer_system,
                             polar_angles='default', azimuthal_angles='default', angular_resolution=None):
    """Evaluate and display the differential scattering cross section as a function of solid angle.

    Args:
        initial_field (smuthi.initial.PlaneWave):   Initial Plane wave
        particle_list (list):                       scattering particles
        layer_system (smuthi.layers.LayerSystem):   stratified medium
        polar_angles (numpy.ndarray or str):        polar angles values (radian).
                                                    if 'default', use smuthi.fields.default_polar_angles
        azimuthal_angles (numpy.ndarray or str):    azimuthal angle values (radian)
                                                    if 'default', use smuthi.fields.default_azimuthal_angles
        angular_resolution (float):                 If provided, angular arrays are generated with this angular
                                                    resolution over the default angular range

    Returns:
        A smuthi.field_expansion.FarField object.
    """
    if not type(initial_field).__name__ == 'PlaneWave':
        raise ValueError('Cross section only defined for plane wave excitation.')

    i_top = layer_system.number_of_layers() - 1
    vacuum_wavelength = initial_field.vacuum_wavelength
    omega = flds.angular_frequency(vacuum_wavelength)
    k_bot = omega * layer_system.refractive_indices[0]
    k_top = omega * layer_system.refractive_indices[-1]

    # read plane wave parameters
    A_P = initial_field.amplitude
    beta_P = initial_field.polar_angle
    if np.cos(beta_P) > 0:
        i_P = 0
        n_P = layer_system.refractive_indices[i_P]
    else:
        i_P = i_top
        n_P = layer_system.refractive_indices[i_P]
    if n_P.imag:
        raise ValueError('plane wave from absorbing layer: cross section undefined')
    else:
        n_P = n_P.real

    initial_intensity = abs(A_P) ** 2 * n_P / 2

    dscs = scattered_far_field(vacuum_wavelength, particle_list, layer_system,
                               polar_angles, azimuthal_angles, angular_resolution)
    dscs.signal_type = 'differential scattering cross section'
    dscs.signal = dscs.signal / initial_intensity

    return dscs

def total_scattering_cross_section(simulation=None, initial_field=None, particle_list=None, layer_system=None,
                                   polar_angles='default', azimuthal_angles='default', angular_resolution=None):
    """Evaluate the total scattering cross section.

    Args:
        simulation (smuthi.Simulation.simulation):        Simulation object (optional)
        initial_field (smuthi.initial_field.PlaneWave):   Initial Plane wave (optional)
        particle_list (list):                       scattering particles (optional)
        layer_system (smuthi.layers.LayerSystem):   stratified medium (optional)
        polar_angles (numpy.ndarray or str):        polar angles values (radian, default None).
                                                    If None, use smuthi.fields.default_polar_angles
        azimuthal_angles (numpy.ndarray or str):    azimuthal angle values (radian, default None).
                                                    If None, use smuthi.fields.default_azimuthal_angles
        angular_resolution (float):                 If provided, angular arrays are generated with this angular
                                                    resolution over the default angular range

    Returns:
        A tuple of smuthi.field_expansion.FarField objects, one for forward scattering (i.e., into the top hemisphere) and one for backward
        scattering (bottom hemisphere).
    """

    if initial_field is None:
        initial_field = simulation.initial_field

    if particle_list is None:
        particle_list = simulation.particle_list

    if layer_system is None:
        layer_system = simulation.layer_system

    azimuthal_angles, polar_angles = _angular_arrays(azimuthal_angles, polar_angles, angular_resolution)

    dscs = scattering_cross_section(initial_field, particle_list, layer_system,
                                    polar_angles, azimuthal_angles, angular_resolution)
    scs = dscs.integral()
    return scs[0] + scs[1]


def extinction_cross_section(simulation=None, initial_field=None, particle_list=None, layer_system=None,
                             only_l=None, only_m=None, only_pol=None, only_tau=None,
                             extinction_direction='both'):
    """Evaluate the extinction cross section.

    Args:
        simulation (smuthi.Simulation.simulation):  Simulation object (optional)
        initial_field (smuthi.initial_field.PlaneWave): Plane wave object (optional)
        particle_list (list): List of smuthi.particles.Particle objects (optional)
        layer_system (smuthi.layers.LayerSystem): Representing the stratified medium
        only_pol (int):  if set to 0 or 1, only this plane wave polarization (0=TE, 1=TM) is considered
        only_tau (int):  if set to 0 or 1, only this spherical vector wave polarization (0 — magnetic, 1 — electric) is
                         considered
        only_l (int):    if set to positive number, only this multipole degree is considered
        only_m (int):    if set non-negative number, only this multipole order is considered
        extinction_direction (string): if set to 'both': return full excinction,
                                       if to 'reflection': extinction of reflected wave,
                                       if to 'transmission': extinction of transmitted wave.
                                       See section on :ref:`ExtinctionCrossSectionAnchor` for details.
    Returns:
        Extinction cross section.
    """
    if initial_field is None:
        initial_field = simulation.initial_field

    if particle_list is None:
        particle_list = simulation.particle_list

    if layer_system is None:
        layer_system = simulation.layer_system

    if not type(initial_field).__name__ == 'PlaneWave':
        raise ValueError('Cross section only defined for plane wave excitation.')

    i_top = layer_system.number_of_layers() - 1
    vacuum_wavelength = initial_field.vacuum_wavelength
    omega = flds.angular_frequency(vacuum_wavelength)
    k_bot = omega * layer_system.refractive_indices[0]
    k_top = omega * layer_system.refractive_indices[-1]

    # read plane wave parameters
    pol_P = initial_field.polarization
    beta_P = initial_field.polar_angle
    alpha_P = initial_field.azimuthal_angle

    if np.cos(beta_P) > 0:
        i_P = 0
        n_P = layer_system.refractive_indices[i_P]
        k_P = k_bot
    else:
        i_P = i_top
        n_P = layer_system.refractive_indices[i_P]
        k_P = k_top
    if n_P.imag:
        raise ValueError('plane wave from absorbing layer: cross section undefined')
    else:
        n_P = n_P.real

    # complex amplitude of initial wave (including phase factor for reference point)
    kappa_P = np.sin(beta_P) * k_P
    kx = np.cos(alpha_P) * kappa_P
    ky = np.sin(alpha_P) * kappa_P
    pm_kz_P = k_P * np.cos(beta_P)
    kvec_P = np.array([kx, ky, pm_kz_P])
    rvec_iP = np.array([0, 0, layer_system.reference_z(i_P)])
    rvec_0 = np.array(initial_field.reference_point)
    ejkriP = np.exp(1j * np.dot(kvec_P, rvec_iP - rvec_0))
    A_P = initial_field.amplitude * ejkriP

    initial_intensity = abs(A_P) ** 2 * n_P / 2

    pwe_scat_top, _ = sf.scattered_field_pwe(vacuum_wavelength, particle_list, layer_system, i_top, kappa_P, alpha_P,
                                             only_l=only_l, only_m=only_m, only_pol=only_pol, only_tau=only_tau)

    _, pwe_scat_bottom = sf.scattered_field_pwe(vacuum_wavelength, particle_list, layer_system, 0, kappa_P, alpha_P,
                                                only_l=only_l, only_m=only_m, only_pol=only_pol, only_tau=only_tau)

    # bottom extinction
    _, pwe_init_bottom = initial_field.plane_wave_expansion(layer_system, 0)
    kz_bot = flds.k_z(k_parallel=kappa_P, k=k_bot)
    gRPbot =  np.squeeze(pwe_init_bottom.coefficients[pol_P])
    g_scat_bottom = np.squeeze(pwe_scat_bottom.coefficients[pol_P])
    P_bot_ext = 4 * np.pi ** 2 * kz_bot / omega * (gRPbot * np.conj(g_scat_bottom)).real
    bottom_extinction_cs = - P_bot_ext / initial_intensity

    # top extinction
    pwe_init_top, _ = initial_field.plane_wave_expansion(layer_system, i_top)
    gRPtop = np.squeeze(pwe_init_top.coefficients[pol_P])
    kz_top = flds.k_z(k_parallel=kappa_P, k=k_top)
    g_scat_top = np.squeeze(pwe_scat_top.coefficients[pol_P])
    P_top_ext = 4 * np.pi ** 2 * kz_top / omega * (gRPtop * np.conj(g_scat_top)).real
    top_extinction_cs = - P_top_ext / initial_intensity

    if extinction_direction == 'both':
        return bottom_extinction_cs + top_extinction_cs
    elif extinction_direction == 'reflection':
        if i_P == 0:
            return bottom_extinction_cs
        else:
            return top_extinction_cs
    elif extinction_direction == 'transmission':
        if i_P == 0:
            return top_extinction_cs
        else:
            return bottom_extinction_cs
    else:
        raise ValueError('extinction_direction can be only \'both\', \
                                        \'reflection\' or \'transmission\'')

"""Functions that assist the user in the choice of suitable numerical simulation parameters."""

import sys
import numpy as np
import matplotlib.pyplot as plt
import smuthi.fields as flds
import smuthi.postprocessing.far_field as ff
from smuthi.fields import angular_frequency
import smuthi.utility.logging as log


def evaluate(simulation, detector):
    """Run a simulation and evaluate the detector.
    Args:
        simulation (smuthi.simulation.Simulation):    simulation object
        detector (method or str):                     Specify a method that accepts a simulation as input and returns
                                                      a float. Otherwise, type "extinction cross section" to use the
                                                      extinction cross section as a detector.

    Returns:
        The detector value (float)
    """
    if detector == "extinction cross section":
        def detector(sim):
            ecs = ff.extinction_cross_section(initial_field=sim.initial_field,
                                              particle_list=sim.particle_list,
                                              layer_system=sim.layer_system)
            return ecs["top"] + ecs["bottom"]

    with log.LoggerMuted():
        simulation.run()
        return detector(simulation)


def converge_l_max(particle,
                   simulation,
                   detector="extinction cross section",
                   tolerance=1e-3,
                   tolerance_steps=2,
                   max_iter=100,
                   start_from_1=True,
                   ax=None):
    """Find suitable multipole cutoff degree `l_max` for a given particle and simulation. The routine starts with the
    current `l_max` of the particle. The value of `l_max` is successively incremented in a loop until the resulting
    relative change in the detector value is smaller than the specified tolerance. The method updates the input
    particle object with the `l_max` value for which convergence has been achieved.

    Args:
        particle (smuthi.particles.Particle):       Particle for which the l_max is incremented
        simulation (smuthi.simulation.Simulation):  Simulation object containing the particle
        detector (function or string):              Function that accepts a simulation object and returns a detector
                                                    value the change of which is used to define convergence.
                                                    Alternatively, use "extinction cross section" (default) to have
                                                    the extinction cross section as the detector value.
        tolerance (float):                          Relative tolerance for the detector value change.
        tolerance_steps (int):                      Number of consecutive steps at which the tolerance must be met,
                                                    valid for neff selection and multipole truncation. Default: 2
        max_iter (int):                             Break convergence loop after that number of iterations, even if
                                                    no convergence has been achieved.
        start_from_1 (logical):                     If true (default), start from `l_max=1`. Otherwise, start from the
                                                    current particle `l_max`.
        ax (np.array of AxesSubplot):               Array of AxesSubplots where to live-plot convergence output

    Returns:
        Detector value of converged or break-off parameter settings.
      """
    print("")
    print("------------------------")
    log.write_blue("Searching suitable l_max")

    if start_from_1:
        particle.l_max = 1
        particle.m_max = 1

    print("Start value: l_max=%i" % particle.l_max)

    current_value = evaluate(simulation, detector)

    x = np.array([particle.l_max])
    y = np.array([current_value.real])
    r = np.array([])
    if ax is not None:
        line1, = ax[0].plot(x, y, '.-')
        line2, = ax[1].plot(x[1:], r, '.-')

    for _ in range(max_iter):
        old_l_max = particle.l_max
        particle.l_max = old_l_max + 1  # l_max increment
        particle.m_max = particle.l_max

        print("---------------------------------------")
        print("Try l_max = %i and m_max=%i" % (particle.l_max, particle.m_max))

        new_value = evaluate(simulation, detector)
        rel_diff = abs(new_value - current_value) / abs(current_value)
        print("Old detector value:", current_value)
        print("New detector value:", new_value)
        print("Relative difference:", rel_diff)
        print("Allowed tolerance:  ", tolerance)

        neff_max = simulation.neff_max
        x = np.append(x, particle.l_max)
        y = np.append(y, new_value.real)
        r = np.append(r, rel_diff)
        if ax is not None:
            line1.set_data(x, y)
            line1.set_label('$n_{eff}^{max} = %g$'%neff_max.real)
            line2.set_data(x[1:], r)
            [( ax.relim(), ax.autoscale_view()) for ax in ax]
            plt.draw()
            plt.pause(0.001)
            ax[0].legend()

        if np.all(r[-tolerance_steps:] < tolerance):  # in this case: discard l_max increment
            particle.l_max = old_l_max
            particle.m_max = old_l_max
            log.write_green("Relative difference smaller than tolerance. Keep l_max = %i" % particle.l_max)
            return current_value
        else:
            current_value = new_value

    log.write_red("No convergence achieved. Keep l_max = %i"%particle.l_max)
    return current_value


def converge_m_max(particle,
                   simulation,
                   detector="extinction cross section",
                   tolerance=1e-3,
                   tolerance_steps=2,
                   current_value=None,
                   ax=None):
    """Find suitable multipole cutoff order `m_max` for a given particle and simulation. The routine starts with the
    current `l_max` of the particle, i.e. with `m_max=l_max`. The value of `m_max` is successively decremented in a loop
    until the resulting relative change in the detector value is larger than the specified tolerance. The method updates
    the input particle object with the so determined `m_max`.

    Args:
        particle (smuthi.particles.Particle):       Particle for which suitable m_max is searched
        simulation (smuthi.simulation.Simulation):  Simulation object containing the particle
        detector (function or string):              Function that accepts a simulation object and returns a detector
                                                    value the change of which is used to define convergence.
                                                    Alternatively, use "extinction cross section" (default) to have
                                                    the extinction cross section as the detector value.
        tolerance (float):                          Relative tolerance for the detector value change.
        tolerance_steps (int):                      Number of consecutive steps at which the tolerance must be met,
                                                    valid for neff selection and multipole truncation. Default: 2
        max_iter (int):                             Break convergence loop after that number of iterations, even if
                                                    no convergence has been achieved.
        current_value (float):                      Start value of detector (for current settings). If not
                                                    specified the method starts with a simulation for the current
                                                    settings.
        ax (np.array of AxesSubplot):               Array of AxesSubplots where to live-plot convergence output

    Returns:
        Detector value of converged or break-off parameter settings.
    """
    print("")
    print("------------------------")
    log.write_blue("Searching suitable m_max")
    print("Start value: m_max=%i"%particle.m_max)

    if current_value is None:
        current_value = evaluate(simulation, detector)

    x = np.array([particle.m_max])
    y = np.array([current_value.real])
    r = np.array([])
    if ax is not None:
        line1, = ax[0].plot(x, y, '.-')
        line2, = ax[1].plot(x[1:], r, '.-')

    for m_max in range(particle.m_max, -1, -1):
        old_m_max = particle.m_max
        particle.m_max = m_max
        print("---------------------------------------")
        print("Try m_max=%i"%particle.m_max)
        new_value = evaluate(simulation, detector)
        rel_diff = abs(new_value - current_value) / abs(current_value)
        print("Old detector value:", current_value)
        print("New detector value:", new_value)
        print("Relative difference:", rel_diff)
        print("Allowed tolerance:  ", tolerance)

        x = np.append(x, particle.m_max)
        y = np.append(y, new_value.real)
        r = np.append(r, rel_diff)
        if ax is not None:
            line1.set_data(x, y)
            line1.set_label('$m_{max}$ for $l_{max} = %g$'%particle.l_max)
            line2.set_data(x[1:], r)
            [ax.relim() for ax in ax]
            [ax.autoscale_view() for ax in ax]
            plt.draw()
            plt.pause(0.001)
            ax[0].legend()

        if np.any(r > tolerance):  # in this case: discard m_max decrement
            particle.m_max = old_m_max
            log.write_green("Relative difference larger than tolerance. Keep m_max = %g"%particle.m_max)
            if ax is not None:
                titlestr = "relative diff < {:g}, keep $m_{{max}} = {:g}$".format(tolerance, particle.m_max)
                ax[1].title.set_text(titlestr)
                ax[1].title.set_color('g')
                plt.draw()
            return current_value
        else:
            titlestr = "no convergence achieved, keep $m_{{max}} = {:g}$".format(particle.m_max)
            current_value = new_value

    log.write_red("No convergence achieved. Keep m_max = %g" % particle.m_max)
    if ax is not None:
        ax[1].title.set_text(titlestr)
        ax[1].title.set_color('r')
        plt.draw()


def converge_multipole_cutoff(simulation,
                              detector="extinction cross section",
                              tolerance=1e-3,
                              tolerance_steps=2,
                              max_iter=100,
                              current_value=None,
                              converge_m=True,
                              ax=None):
    """Find suitable multipole cutoff degree `l_max` and order `m_max` for all particles in a given simulation object.
    The method updates the input simulation object with the so determined multipole truncation values.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object
        detector (function or string):              Function that accepts a simulation object and returns a detector
                                                    value the change of which is used to define convergence.
                                                    Alternatively, use "extinction cross section" (default) to have
                                                    the extinction cross section as the detector value.
        tolerance (float):                          Relative tolerance for the detector value change.
        tolerance_steps (int):                      Number of consecutive steps at which the tolerance must be met,
                                                    valid for neff selection and multipole truncation. Default: 2
        max_iter (int):                             Break convergence loops after that number of iterations, even if
                                                    no convergence has been achieved.
        current_value (float):                      Start value of detector (for current settings). If not
                                                    specified the method starts with a simulation for the current
                                                    settings.
        converge_m (logical):                       If false, only converge `l_max`, but keep `m_max=l_max`. Default
                                                    is true.
        ax (np.array of AxesSubplot):               Array of AxesSubplots where to live-plot convergence output

    Returns:
        Detector value of converged or break-off parameter settings.
    """
    print("")
    print("-----------------------------------")
    log.write_blue("Searching suitable multipole cutoff")

    for i, particle in enumerate(simulation.particle_list):
        print("")
        print("-------------------------------")
        log.write_blue("Checking particle number %i"%i)
        with log.LoggerIndented():
            current_value = converge_l_max(particle,
                                           simulation,
                                           detector,
                                           tolerance,
                                           tolerance_steps,
                                           max_iter,
                                           ax=ax)

            if converge_m:
                current_value = converge_m_max(particle,
                                               simulation,
                                               detector,
                                               tolerance,
                                               tolerance_steps,
                                               current_value,
                                               ax=ax)

    return current_value


def update_contour(simulation, neff_imag=5e-3, neff_max=None, neff_max_offset=0.5, neff_resolution=2e-3):
    """Update the default `k_parallel` arrays in smuthi.fields with a newly constructed Sommerfeld integral
    contours.

    Args:
        simulation (smuthi.simulation.Simulation):      Simulation object
        neff_imag (float):                              Extent of the contour into the negative imaginary direction
                                                        (in terms of effective refractive index, n_eff=kappa/omega).
        neff_max (float):                               Truncation value of contour (in terms of effective refractive
                                                        index).
        neff_max_offset (float):                        If no value for `neff_max` is specified, use the last estimated
                                                        singularity location plus this value (in terms of effective
                                                        refractive index).
        neff_resolution (float):                        Discretization of the contour (in terms of eff. refractive
                                                        index).
        update_default_contours (logical)               If true, overwrite the default contours in smuthi.fields module.
                                                        Otherwise, overwrite simulation.k_parallel array
    """
    simulation.neff_imag = neff_imag
    simulation.neff_max = neff_max
    simulation.neff_max_offset = neff_max_offset
    simulation.neff_resolution = neff_resolution
    simulation.set_default_initial_field_contour()
    simulation.set_default_Sommerfeld_contour()

    simulation.k_parallel = flds.default_Sommerfeld_k_parallel_array


def converge_neff_max(simulation,
                      detector="extinction cross section",
                      tolerance=1e-3,
                      tolerance_factor=0.1,
                      tolerance_steps=2,
                      max_iter=20,
                      neff_imag=1e-2,
                      neff_resolution=2e-3,
                      neff_max_increment=0.5,
                      converge_lm=True,
                      ax=None):
    """Find a suitable truncation value for the multiple scattering Sommerfeld integral contour and update the
    simulation object accordingly.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object
        detector (function or string):              Function that accepts a simulation object and returns a detector
                                                    value the change of which is used to define convergence.
                                                    Alternatively, use "extinction cross section" (default) to have
                                                    the extinction cross section as the detector value.
        tolerance (float):                          Relative tolerance for the detector value change.
        tolerance_factor (float):                   During neff selection, a smaller tolerance should be allowed to
                                                    avoid fluctuations of the order of ~tolerance which would
                                                    compromise convergence. Default: 0.1
        tolerance_steps (int):                      Number of consecutive steps at which the tolerance must be met,
                                                    valid for neff selection and multipole truncation. Default: 2
        max_iter (int):                             Break convergence loops after that number of iterations, even if
                                                    no convergence has been achieved.
        neff_imag (float):                          Extent of the contour into the negative imaginary direction
                                                    (in terms of effective refractive index, n_eff=kappa/omega).
        neff_resolution (float):                    Discretization of the contour (in terms of eff. refractive
                                                    index).
        neff_max_increment (float):                 Increment the neff_max parameter with that step size
        converge_lm (logical):                      If set to true, update multipole truncation during each step
                                                    (this takes longer time, but is necessary for critical use cases
                                                    like flat particles on a substrate)
        ax (np.array of AxesSubplot):               Array of AxesSubplots where to live-plot convergence output

    Returns:
        Detector value for converged settings.
    """

    print("")
    print("---------------------------")
    log.write_blue("Searching suitable neff_max")

    update_contour(simulation=simulation, neff_imag=neff_imag, neff_max_offset=0, neff_resolution=neff_resolution)

    simulation.k_parallel = flds.reasonable_Sommerfeld_kpar_contour(
        vacuum_wavelength=simulation.initial_field.vacuum_wavelength,
        layer_refractive_indices=simulation.layer_system.refractive_indices,
        neff_imag=neff_imag,
        neff_max_offset=0,
        neff_resolution=neff_resolution)

    neff_max = simulation.k_parallel[-1] / angular_frequency(simulation.initial_field.vacuum_wavelength)
    simulation.neff_max = neff_max
    print("Starting value: neff_max=%f"%neff_max.real)

    if converge_lm:
        with log.LoggerIndented():
            current_value = converge_multipole_cutoff(simulation=simulation,
                                                      detector=detector,
                                                      tolerance=tolerance*tolerance_factor,
                                                      tolerance_steps=2,
                                                      max_iter=max_iter,
                                                      converge_m=False,
                                                      ax=ax)
    else:
        current_value = evaluate(simulation, detector)

    x = np.array([neff_max.real])
    y = np.array([current_value.real])
    r = np.array([])
    if ax is not None:
        line1, = ax[-2].plot(x, y, 'k.-')
        line2, = ax[-1].plot(x[1:], r, 'k.-')
        ax[-1].set_xlabel('$n_{eff}^{max}$')
        ax[-2].set_title('$n_{eff}^{max}$ selection')

    for _ in range(max_iter):
        old_neff_max = neff_max
        neff_max = neff_max + neff_max_increment
        update_contour(simulation=simulation, neff_imag=neff_imag, neff_max=neff_max, neff_resolution=neff_resolution)

        print("---------------------------------------")
        print("Try neff_max = %f"%neff_max.real)

        if converge_lm:
            with log.LoggerIndented():
                new_value = converge_multipole_cutoff(simulation=simulation,
                                                      detector=detector,
                                                      tolerance=tolerance*tolerance_factor,
                                                      tolerance_steps=2,
                                                      max_iter=max_iter,
                                                      current_value=current_value,
                                                      converge_m=False,
                                                      ax=ax)
        else:
            new_value = evaluate(simulation, detector)

        rel_diff = abs(new_value - current_value) / abs(current_value)
        print("Old detector value:", current_value)
        print("New detector value:", new_value)
        print("Relative difference:", rel_diff)
        print("Allowed tolerance:  ", tolerance)

        x = np.append(x, neff_max.real)
        y = np.append(y, new_value.real)
        r = np.append(r, rel_diff)
        if ax is not None:
            line1.set_data(x, y)
            line1.set_label('$n_{eff}^{max}$')
            line2.set_data(x[1:], r)
            [( ax.relim(), ax.autoscale_view()) for ax in ax]
            plt.draw()
            plt.pause(0.001)
            [ax.legend() for ax in ax[::-2]]

        if np.all(r[-tolerance_steps:] < tolerance):  # in this case: discard l_max increment
            neff_max = old_neff_max
            update_contour(simulation=simulation, neff_imag=neff_imag, neff_max=neff_max, neff_resolution=neff_resolution)
            log.write_green("Relative difference smaller than tolerance. Keep neff_max = %g"%neff_max.real)
            if ax is not None:
                titlestr = "relative diff < {:g}, keep $n_{{eff}}^{{max}} = {:g}$".format(tolerance, neff_max.real)
                ax[-1].title.set_text(titlestr)
                ax[-1].title.set_color('g')
                plt.draw()
            return current_value
        else:
            titlestr = "no convergence achieved, keep $n_{{eff}}^{{max}} = {:g}$".format(neff_max.real)
            current_value = new_value

    log.write_red("No convergence achieved. Keep neff_max = %g" % neff_max.real)
    if ax is not None:
        ax[-1].title.set_text(titlestr)
        ax[-1].title.set_color('r')
        plt.draw()

    return None


def converge_neff_resolution(simulation,
                       detector="extinction cross section",
                       tolerance=1e-3,
                       max_iter=20,
                       neff_imag=1e-2,
                       neff_max=None,
                       neff_resolution=1e-2,
                       ax=None):
    """Find a suitable discretization step size for the multiple scattering Sommerfeld integral contour and update the
    simulation object accordingly.

    Args:
        simulation (smuthi.simulation.Simulation):    Simulation object
        detector (function or string):                Function that accepts a simulation object and returns a detector
                                                      value the change of which is used to define convergence.
                                                      Alternatively, use "extinction cross section" (default) to have
                                                      the extinction cross section as the detector value.
        tolerance (float):                            Relative tolerance for the detector value change.
        max_iter (int):                               Break convergence loops after that number of iterations, even if
                                                      no convergence has been achieved.
        neff_imag (float):                            Extent of the contour into the negative imaginary direction
                                                      (in terms of effective refractive index, n_eff=kappa/omega).
        neff_max (float):                             Truncation value of contour (in terms of effective refractive
                                                      index).
        neff_resolution (float):                      Discretization of the contour (in terms of eff. refractive
                                                      index) - start value for iteration
        ax (np.array of AxesSubplot):                 Array of AxesSubplots where to live-plot convergence output

    Returns:
        Detector value for converged settings.
    """
    print("")
    print("-----------------------")
    log.write_blue("Find suitable neff_resolution")

    update_contour(simulation=simulation, neff_imag=neff_imag, neff_max=neff_max, neff_resolution=neff_resolution)
    print("Starting value: neff_resolution=%f" % neff_resolution)
    current_value = evaluate(simulation, detector)

    if ax is not None:
        x = np.array([])
        y = np.array([])
        r = np.array([])
        init_x = neff_resolution
        init_y = current_value.real
        line1, = ax[0].plot(np.insert(x, 0, init_x), np.insert(y.real, 0, init_y),'.-')
        line2, = ax[1].plot(x,r,'.-')

    for _ in range(max_iter):
        old_neff_resolution = neff_resolution
        neff_resolution = neff_resolution / 2.0
        update_contour(simulation=simulation, neff_imag=neff_imag, neff_max=neff_max, neff_resolution=neff_resolution)

        print("---------------------------------------")
        print("Try neff_resolution = %f" % neff_resolution)

        new_value = evaluate(simulation, detector)

        rel_diff = abs(new_value - current_value) / abs(current_value)
        print("Old detector value:", current_value)
        print("New detector value:", new_value)
        print("Relative difference:", rel_diff)
        print("Allowed tolerance:  ", tolerance)

        if ax is not None:
            x = np.append(x, neff_resolution)
            y = np.append(y, new_value)
            r = np.append(r, rel_diff)
            line1.set_data(np.insert(x, 0, init_x), np.insert(y.real, 0, init_y.real))
            line1.set_label('$n_{{eff}}^{{max}} = {}$, $l_{{max}} = {}$, $m_{{max}} = {}$'.\
                            format(neff_max.real.round(decimals=2),
                                   simulation.particle_list[0].l_max,
                                   simulation.particle_list[0].m_max,
                                   neff_resolution))
            line2.set_data(x, r)
            [( ax.relim(), ax.autoscale_view()) for ax in ax]
            plt.draw()
            plt.pause(0.001)
            ax[0].legend()

        if rel_diff < tolerance:  # in this case: discard halfed neff_resolution
            neff_resolution = old_neff_resolution
            update_contour(simulation=simulation, neff_imag=neff_imag, neff_max=neff_max, neff_resolution=neff_resolution)
            log.write_green("Relative difference smaller than tolerance. Keep neff_resolution = %g" % neff_resolution)
            if ax is not None:
                titlestr = "relative diff < {:g}, keep $\delta n_{{eff}} = {:g}$".format(tolerance, neff_resolution)
                ax[1].title.set_text(titlestr)
                ax[1].title.set_color('g')
                plt.draw()
            return current_value
        else:
            titlestr = "no convergence achieved, keep $\delta n_{{eff}} = {:g}$".format(neff_max.real)
            current_value = new_value

    log.write_red("No convergence achieved. Keep neff_resolution = %g" % neff_resolution)
    if ax is not None:
        ax[1].title.set_text(titlestr)
        ax[1].title.set_color('r')
        plt.draw()

    return None


def select_numerical_parameters(simulation,
                                detector="extinction cross section",
                                tolerance=1e-3,
                                tolerance_factor=0.1,
                                tolerance_steps=2,
                                max_iter=20,
                                neff_imag=1e-2,
                                neff_resolution=1e-2,
                                select_neff_max=True,
                                neff_max_increment=0.5,
                                neff_max=None,
                                select_neff_resolution=True,
                                select_multipole_cutoff=True,
                                relative_convergence=True,
                                show_plot=True):
    """Trigger automatic selection routines for various numerical parameters.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object from which parameters are read and into which
                                                    results are stored.
        detector (function or string):              Function that accepts a simulation object and returns a detector
                                                    value the change of which is used to define convergence.
                                                    Alternatively, use "extinction cross section" (default) to have
                                                    the extinction cross section as the detector value.
        tolerance (float):                          Relative tolerance for the detector value change.
        tolerance_factor (float):                   During neff selection, a smaller tolerance should be allowed to
                                                    avoid fluctuations of the order of ~tolerance which would
                                                    compromise convergence. Default: 0.1
        tolerance_steps (int):                      Number of consecutive steps at which the tolerance must be met,
                                                    valid for neff selection and multipole truncation. Default: 2
        max_iter (int):                             Break convergence loops after that number of iterations, even if
                                                    no convergence has been achieved.
        neff_imag (float):                          Extent of the contour into the negative imaginary direction
                                                    (in terms of effective refractive index, n_eff=kappa/omega).
        neff_resolution (float):                    Discretization of the contour (in terms of eff. refractive
                                                    index) - if `select_neff_resolution` is true, this value will be
                                                    eventually overwritten. However, it is required in any case.
                                                    Default: 1e-2
        select_neff_max (logical):                  If set to true (default), the Sommerfeld integral truncation
                                                    parameter `neff_max` is determined automatically with the help
                                                    of a Cauchy convergence criterion.
        neff_max_increment (float):                 Only needed if `select_neff_max` is true.
                                                    Step size with which `neff_max` is incremented.
        neff_max (float):                           Only needed if `select_neff_max` is false.
                                                    Truncation value of contour (in terms of effective refractive
                                                    index).
        select_neff_resolution (logical):           If set to true (default), the Sommerfeld integral discretization
                                                    parameter `neff_resolution` is determined automatically with the
                                                    help of a Cauchy convergence criterion.
        select_multipole_cutoff (logical):          If set to true (default), the multipole expansion cutoff
                                                    parameters `l_max` and `m_max` are determined automatically with
                                                    the help of a Cauchy convergence criterion.
        relative_convergence (logical):             If set to true (default), the `neff_max` convergence and the
                                                    `l_max` and `m_max` convergence routine are performed in the
                                                    spirit of relative convergence, i.e., the multipole expansion
                                                    convergence is checked again for each value of the Sommerfeld
                                                    integral truncation. This takes more time, but is required at
                                                    least in the case of flat particles near interfaces.
    """

    def _init_fig(xlabel='x', ylabel='y', title=None, tol=None, allowedtol=None, cols=1):
        """
        Inner utility function returning figure and axes handles initialized
        following a common template for all parameter selection runs.

        Args:
            xlabel (string):    column-wise common x axis label
            ylabel (string):    y axis label for the detector quantity
            title  (string):    title string for the (lefmost, if cols=2) detector panel
            tol (float):        user-set tolerance level to be drawn in relative difference panel
            allowedtol (float): when specified, also draw the actual tolerance level for the current run
            cols (int):         number of subplot columns. select_neff_max() accepts 2
        """
        fig, ax_array = plt.subplots(2, cols, sharex='col', figsize=(cols*6.4,8), gridspec_kw={'height_ratios': [3, 1]})
        ax_array = ax_array.flatten(order='F')
        ax_array[1].set_xlabel(xlabel)
        ax_array[0].set_ylabel(ylabel)
        ax_array[1].set_ylabel('relative difference')
        ax_array[0].set_title(title)
        [ax.axhline(tol, color='grey', linestyle='dashed', label='tolerance') for ax in ax_array[::-2]]
        if allowedtol is not None:
            ax_array[1].axhline(allowedtol, color='grey', linestyle='dotted', label='allowed tolerance')
        [( ax.set_yscale('log'), ax.grid(), ax.legend() ) for ax in ax_array[::-2]]

        return fig, ax_array

    print("")
    print("----------------------------------------")
    print("Starting automatic paramateter selection")
    print("----------------------------------------")

    if tolerance_factor <= 0 or tolerance_factor > 1:
        raise ValueError('please provide 0 < tolerance_factor <= 1')

    if select_neff_max:
        if show_plot:
            plt.ion()
            _, ax_array = _init_fig('$l_{max}$', detector, 'multipole cutoff', tolerance, 0.1*tolerance, 2)
        else:
            ax_array = None
        converge_neff_max(simulation=simulation,
                          detector=detector,
                          tolerance=tolerance,
                          tolerance_factor=tolerance_factor,
                          tolerance_steps=tolerance_steps,
                          max_iter=max_iter,
                          neff_imag=neff_imag,
                          neff_resolution=neff_resolution,
                          neff_max_increment=neff_max_increment,
                          converge_lm=relative_convergence,
                          ax=ax_array)
        neff_max = simulation.neff_max

    if select_multipole_cutoff:
        if show_plot:
            plt.ion()
            _, ax_array = _init_fig('multipole order', detector, 'multipole cutoff selection', tolerance)
        else:
            ax_array = None
        converge_multipole_cutoff(simulation=simulation,
                                  detector=detector,
                                  tolerance=tolerance,
                                  max_iter=max_iter,
                                  ax=ax_array)

    if select_neff_resolution:
        if show_plot:
            plt.ion()
            _, ax_array = _init_fig('$\delta n_{eff}$', detector, '$\delta n_{eff}$ selection', tolerance)
        else:
            ax_array = None
        converge_neff_resolution(simulation=simulation,
                                 detector=detector,
                                 tolerance=tolerance,
                                 max_iter=max_iter,
                                 neff_imag=neff_imag,
                                 neff_max=neff_max,
                                 ax=ax_array)

    if show_plot:
        plt.ioff()

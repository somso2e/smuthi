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
        simulation (smuthi.simulation.Simulation):  simulation object
        detector (method or str):                   Specify a method that accepts a simulation as input and returns
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
            return ecs
    elif detector == "total scattering cross section":
        def detector(sim):
            scs = ff.total_scattering_cross_section(initial_field=sim.initial_field,
                                                    particle_list=sim.particle_list,
                                                    layer_system=sim.layer_system)
            return scs
    elif detector == "integrated scattered far field":
        def detector(sim):
            # _, _, scff = ff.total_far_field(initial_field=sim.initial_field,
            scff = ff.scattered_far_field(vacuum_wavelength=sim.initial_field.vacuum_wavelength,
                                          particle_list=sim.particle_list,
                                          layer_system=sim.layer_system)
            iscff = scff.integral()
            return iscff[0] + iscff[1]

    with log.LoggerMuted():
        simulation.run()
        return detector(simulation)

def update_lmax_mmax(simulation, l_max):
    """ Assign the same l_max and m_max = l_max to all particles in simulation """
    for particle in simulation.particle_list:
        particle.l_max = l_max
        particle.m_max = l_max

def converge_l_max(simulation,
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
        simulation (smuthi.simulation.Simulation):  Simulation object containing the particle
        detector (function or string):              Function that accepts a simulation object and returns a detector
                                                    value the change of which is used to define convergence.
                                                    Alternatively, use "extinction cross section" (default) to have
                                                    the extinction cross section as the detector value.
        tolerance (float):                          Relative tolerance for the detector value change.
        tolerance_steps (int):                      Number of consecutive steps at which the tolerance must be met
                                                    during multipole truncation convergence. Default: 2
        max_iter (int):                             Break convergence loop after that number of iterations, even if
                                                    no convergence has been achieved.
        start_from_1 (logical):                     If true (default), start from `l_max=1`. Otherwise, start from the
                                                    current particle `l_max`.
        ax (np.array of AxesSubplot):               Array of AxesSubplots where to live-plot convergence output

    Returns:
        A 3-tuple containing 
          - detector value of converged or break-off parameter settings.
          - series of lmax values
          - the detector values for the given lmax values
      """

    print("")
    print("------------------------")
    log.write_blue("Searching suitable l_max")

    if start_from_1:
        l_max = 1
        update_lmax_mmax(simulation, l_max)
    else:
        l_max = simulation.particle_list[0].l_max # assuming that all particles have the same l_max

    print("Start value: l_max=%i" % l_max)

    current_value = evaluate(simulation, detector)

    x = np.array([l_max])
    y = np.array([current_value])
    r = np.array([])
    if ax is not None:
        line1, = ax[0].plot(x, np.real(y), '.-')
        line2, = ax[1].plot(x[1:], r, '.-')

    for _ in range(max_iter):
        l_max += 1  # l_max increment
        update_lmax_mmax(simulation, l_max)

        print("---------------------------------------")
        print("Try l_max = %i and m_max=%i" % (l_max, l_max))

        new_value = evaluate(simulation, detector)
        rel_diff = abs(new_value - current_value) / abs(current_value)
        print("Old detector value:", current_value)
        print("New detector value:", new_value)
        print("Relative difference:", rel_diff)
        print("Allowed tolerance:  ", tolerance)

        neff_max = simulation.neff_max
        x = np.append(x, l_max)
        y = np.append(y, new_value)
        r = np.append(r, rel_diff)
        if ax is not None:
            line1.set_data(x, np.real(y))
            if neff_max is not None:
                line1.set_label('$n_{eff}^{max} = %g$'%neff_max.real)
            else:
                line1.set_label('$l_{max}$')
            line2.set_data(x[1:], r)
            for axi in ax:
                axi.relim()
                axi.autoscale_view()
            plt.draw()
            plt.pause(0.001)
            ax[0].legend()

        if np.all(r[-tolerance_steps:] < tolerance) and len(r) >= tolerance_steps: # then, discard l_max increment
            old_l_max = l_max - tolerance_steps
            update_lmax_mmax(simulation, old_l_max)
            log.write_green("Relative difference smaller than tolerance. Keep l_max = %i" % old_l_max)
            return y[-tolerance_steps-1], x, y
        else:
            current_value = new_value

    log.write_red("No convergence achieved. Keep l_max = %i"%l_max)
    return current_value, x, y


def converge_m_max(simulation,
                   detector="extinction cross section",
                   tolerance=1e-3,
                   target_value=None,
                   ax=None):
    """Find suitable multipole cutoff order `m_max` for a given particle and simulation. The routine starts with the
    current `l_max` of the particle, i.e. with `m_max=l_max`. The value of `m_max` is successively decremented in a loop
    until the resulting relative change in the detector value is larger than the specified tolerance. The method updates
    the input particle object with the so determined `m_max`.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object containing the particle
        detector (function or string):              Function that accepts a simulation object and returns a detector
                                                    value the change of which is used to define convergence.
                                                    Alternatively, use "extinction cross section" (default) to have
                                                    the extinction cross section as the detector value.
        tolerance (float):                          Relative tolerance for the detector value change.
        max_iter (int):                             Break convergence loop after that number of iterations, even if
                                                    no convergence has been achieved.
        target_value (float):                       If available (typically from preceding neff selection procedure),
                                                    use as target detector value
        ax (np.array of AxesSubplot):               Array of AxesSubplots where to live-plot convergence output

    Returns:
        Detector value of converged or break-off parameter settings.
    """

    def update_mmax(simulation, m_max):
        """ Assign the same m_max to all particles in simulation """
        for particle in simulation.particle_list:
            particle.m_max = m_max

    print("")
    print("------------------------")
    log.write_blue("Searching suitable m_max")
    m_max = simulation.particle_list[0].m_max # assuming that all particles have the same l_max
    print("Start value: m_max=%i"%m_max)

    current_value = evaluate(simulation, detector) if target_value is None else target_value

    x = np.array([])
    y = np.array([])
    r = np.array([])
    if ax is not None:
        ax[0].axhline(current_value, color='black', linestyle='dashed',
                      label='target value')
        line1, = ax[0].plot(x, y, '.-')
        line2, = ax[1].plot(x, r, '.-')

    for m_max in range(m_max, -1, -1):
        old_m_max = simulation.particle_list[0].m_max
        update_mmax(simulation, m_max)
        print("---------------------------------------")
        print("Try m_max=%i"%m_max)
        new_value = evaluate(simulation, detector)
        rel_diff = abs(new_value - current_value) / abs(current_value)
        print("Old detector value:", current_value)
        print("New detector value:", new_value)
        print("Relative difference:", rel_diff)
        print("Allowed tolerance:  ", tolerance)

        x = np.append(x, m_max)
        y = np.append(y, new_value.real)
        r = np.append(r, rel_diff)
        if ax is not None:
            line1.set_data(x, y)
            line1.set_label('$m_{max}$ for $l_{max} = %g$'%simulation.particle_list[0].l_max)
            line2.set_data(x, r)
            for axi in ax:
                axi.relim()
                axi.autoscale_view()
            plt.draw()
            plt.pause(0.001)
            ax[0].legend()

        if np.any(r > tolerance):  # in this case: discard m_max decrement
            update_mmax(simulation, old_m_max)
            log.write_green("Relative difference larger than tolerance. Keep m_max = %g"%old_m_max)
            if ax is not None:
                titlestr = "relative diff > {:g}, keep $m_{{max}} = {:g}$".format(tolerance, old_m_max)
                ax[1].title.set_text(titlestr)
                ax[1].title.set_color('g')
                plt.draw()
            return current_value
        else:
            titlestr = "no convergence achieved, keep $m_{{max}} = {:g}$".format(m_max)
            current_value = new_value if target_value is None else target_value

    log.write_red("No convergence achieved. Keep m_max = %g" % m_max)
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
                              l_max_list=None,
                              detector_value_list=None,
                              converge_m=True,
                              ax=None):
    """Find suitable multipole cutoff degree `l_max` and order `m_max` for all particles in a given simulation object.
    The method updates the input simulation object with the so determined multipole truncation values.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object
        detector (function or string):              Function that accepts a simulation object and returns a detector
                                                    value the change of which is used to define convergence
                                                    Alternatively, use "extinction cross section" (default) to have
                                                    the extinction cross section as the detector value
        tolerance (float):                          Relative tolerance for the detector value change
        tolerance_steps (int):                      Number of consecutive steps at which the tolerance must be met
                                                    during multipole truncation convergence. Default: 2
        max_iter (int):                             Break convergence loops after that number of iterations, even if
                                                    no convergence has been achieved
        current_value (float):                      If specified, skip l_max run and use this value for the
                                                    resulting detector value.
                                                    Otherwise, start with l_max run.
        l_max_list (list):                          If current_value was specified, the l_max run is skipped.
                                                    Then, this list is returned as the second item in the returned tuple.
        detector_value_list (list):                 If current_value was specified, the l_max run is skipped.
                                                    Then, this list is returned as the third item in the returned tuple.
        converge_m (logical):                       If false, only converge `l_max`, but keep `m_max=l_max`. Default
                                                    is true
        ax (np.array of AxesSubplot):               Array of AxesSubplots where to live-plot convergence output

    Returns:
        A 3-tuple containing
          - detector value of converged or break-off parameter settings.
          - series of lmax values
          - the detector values for the given lmax values
    """
    print("")
    print("-----------------------------------")
    log.write_blue("Searching suitable multipole cutoff")

    with log.LoggerIndented():
        if current_value is None:
            current_value, lmx, values = converge_l_max(simulation=simulation,
                                                        detector=detector,
                                                        tolerance=tolerance,
                                                        tolerance_steps=tolerance_steps,
                                                        max_iter=max_iter,
                                                        start_from_1=True, # start_from_1
                                                        ax=ax)
        else:
            log.write_green("Keep l_max = %i"%simulation.particle_list[0].l_max)
            lmx = l_max_list
            values = detector_value_list

        if converge_m:
            current_value = converge_m_max(simulation=simulation,
                                           detector=detector,
                                           tolerance=tolerance,
                                           target_value=current_value, # passing current_value as the target_value
                                           ax=ax)

    return current_value, lmx, values


def update_contour(simulation, neff_imag=5e-3, neff_max=None, neff_max_offset=0.5, neff_resolution=2e-3):
    """Update the default `k_parallel` arrays in smuthi.fields with a newly constructed Sommerfeld integral
    contour, and set the simulation object to use the default contour for particle coupling.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object
        neff_imag (float):                          Extent of the contour into the negative imaginary direction
                                                    (in terms of effective refractive index, n_eff=kappa/omega).
        neff_max (float):                           Truncation value of contour (in terms of effective refractive
                                                    index).
        neff_max_offset (float):                    If no value for `neff_max` is specified, use the last estimated
                                                    singularity location plus this value (in terms of effective
                                                    refractive index).
        neff_resolution (float):                    Discretization of the contour (in terms of effective refractive
                                                    index).
    """
    simulation.neff_imag = neff_imag
    simulation.neff_max = neff_max
    simulation.neff_max_offset = neff_max_offset
    simulation.neff_resolution = neff_resolution
    simulation.set_default_initial_field_contour()
    simulation.set_default_Sommerfeld_contour()

    simulation.k_parallel = "default"


def converge_neff_max(simulation,
                      detector="extinction cross section",
                      tolerance=1e-3,
                      tolerance_factor=0.1,
                      tolerance_steps=2,
                      max_iter=30,
                      neff_imag=1e-2,
                      neff_resolution=2e-3,
                      neff_max_increment=0.5,
                      neff_max_offset=0,
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
        tolerance_steps (int):                      Number of consecutive steps at which the tolerance must be met
                                                    during multipole truncation convergence. Default: 2
        max_iter (int):                             Break convergence loops after that number of iterations, even if
                                                    no convergence has been achieved.
        neff_imag (float):                          Extent of the contour into the negative imaginary direction
                                                    (in terms of effective refractive index, n_eff=kappa/omega).
        neff_resolution (float):                    Discretization of the contour (in terms of eff. refractive
                                                    index).
        neff_max_increment (float):                 Increment the neff_max parameter with that step size
        neff_max_offset (float):                    Start neff_max selection from the last estimated singularity
                                                    location plus this value (in terms of effective refractive index)
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

    update_contour(simulation=simulation, neff_imag=neff_imag, neff_max_offset=neff_max_offset, neff_resolution=neff_resolution)

    simulation.k_parallel = flds.reasonable_Sommerfeld_kpar_contour(
        vacuum_wavelength=simulation.initial_field.vacuum_wavelength,
        layer_refractive_indices=simulation.layer_system.refractive_indices,
        neff_imag=neff_imag,
        neff_max_offset=neff_max_offset,
        neff_resolution=neff_resolution)

    neff_max = simulation.k_parallel[-1] / angular_frequency(simulation.initial_field.vacuum_wavelength)
    simulation.neff_max = neff_max
    print("Starting value: neff_max=%f"%neff_max.real)

    if converge_lm:
        with log.LoggerIndented():
            current_value, lmx_list, values = converge_multipole_cutoff(simulation=simulation,
                                                                        detector=detector,
                                                                        tolerance=tolerance*tolerance_factor,
                                                                        tolerance_steps=tolerance_steps,
                                                                        max_iter=max_iter,
                                                                        converge_m=False,
                                                                        ax=ax)
    else:
        current_value = evaluate(simulation, detector)

    x = np.array([neff_max.real])
    y = np.array([current_value.real])
    r = np.array([])
    lmx_lists = [lmx_list]
    lmx_values_lists = [values]
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
                new_value, lmx_list, values = converge_multipole_cutoff(simulation=simulation,
                                                                        detector=detector,
                                                                        tolerance=tolerance*tolerance_factor,
                                                                        tolerance_steps=tolerance_steps,
                                                                        max_iter=max_iter,
                                                                        current_value=None,
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
        lmx_lists.append(lmx_list)
        lmx_values_lists.append(values)
        if ax is not None:
            line1.set_data(x, y)
            line1.set_label('$n_{eff}^{max}$')
            line2.set_data(x[1:], r)
            for axi in ax:
                axi.relim()
                axi.autoscale_view()
            plt.draw()
            plt.pause(0.001)
            for axi in ax[::-2]:
                axi.legend()

        if rel_diff < tolerance:  # then, discard neff_max increment
            neff_max = old_neff_max
            update_contour(simulation=simulation, neff_imag=neff_imag, neff_max=neff_max, neff_resolution=neff_resolution)
            log.write_green(
                "Relative difference smaller than tolerance. Keep neff_max = %g" % neff_max.real)
            if converge_lm:
                for idx, lmx in enumerate(lmx_lists[-2]):
                    update_lmax_mmax(simulation, lmx)
                    rel_diff_1 = abs(
                        lmx_values_lists[-2][idx] - lmx_values_lists[-2][idx + 1]) / abs(lmx_values_lists[-2][idx + 1])
                    rel_diff_2 = abs(
                        lmx_values_lists[-2][idx] - lmx_values_lists[-2][idx + 2]) / abs(lmx_values_lists[-2][idx + 2])
                    current_value = lmx_values_lists[-2][idx]
                    if rel_diff_1 < tolerance and rel_diff_2 < tolerance:
                        break

            if ax is not None:
                titlestr = "relative diff < {:g}, keep $n_{{eff}}^{{max}} = {:g}$".format(tolerance, neff_max.real)
                ax[-1].title.set_text(titlestr)
                ax[-1].title.set_color('g')
                plt.draw()

            return current_value  # most accurate result for old_neff_max
        else:
            titlestr = "no convergence achieved, keep $n_{{eff}}^{{max}} = {:g}$".format(neff_max.real)
            current_value = new_value

    log.write_red("No convergence achieved. Keep neff_max = %g" % neff_max.real)
    if ax is not None:
        ax[-1].title.set_text(titlestr)
        ax[-1].title.set_color('r')
        plt.draw()

    return current_value


def converge_neff_resolution(simulation,
                       detector="extinction cross section",
                       tolerance=1e-3,
                       max_iter=30,
                       neff_imag=1e-2,
                       neff_max=None,
                       neff_resolution=1e-2,
                       ax=None):
    """Find a suitable discretization step size for the multiple scattering Sommerfeld integral contour and update
    the simulation object accordingly.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object
        detector (function or string):              Function that accepts a simulation object and returns a detector
                                                    value the change of which is used to define convergence.
                                                    Alternatively, use "extinction cross section" (default) to have
                                                    the extinction cross section as the detector value.
        tolerance (float):                          Relative tolerance for the detector value change.
        max_iter (int):                             Break convergence loops after that number of iterations, even if
                                                    no convergence has been achieved.
        neff_imag (float):                          Extent of the contour into the negative imaginary direction
                                                    (in terms of effective refractive index, n_eff=kappa/omega).
        neff_max (float):                           Truncation value of contour (in terms of effective refractive
                                                    index).
        neff_resolution (float):                    Discretization of the contour (in terms of eff. refractive
                                                    index) - start value for iteration
        ax (np.array of AxesSubplot):               Array of AxesSubplots where to live-plot convergence output

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
            line1.set_label('$n_{{eff}}^{{max}} = {:g}$, $l_{{max}} = {}$, $m_{{max}} = {}$'.\
                            format(simulation.neff_max.real,
                                   simulation.particle_list[0].l_max,
                                   simulation.particle_list[0].m_max))
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


def converge_angular_resolution(simulation,
                                detector="extinction cross section",
                                tolerance=1e-3,
                                max_iter=30,
                                ax=None):
    """Find a suitable discretization step size for the default angular arrays used for plane wave expansions.

    Args:
        simulation (smuthi.simulation.Simulation):  Simulation object
        detector (function or string):              Function that accepts a simulation object and returns a detector
                                                    value the change of which is used to define convergence.
                                                    Alternatively, use "extinction cross section" (default) to have
                                                    the extinction cross section as the detector value.
        tolerance (float):                          Relative tolerance for the detector value change.
        max_iter (int):                             Break convergence loops after that number of iterations, even if
                                                    no convergence has been achieved.
        ax (np.array of AxesSubplot):               Array of AxesSubplots where to live-plot convergence output

    Returns:
        Detector value for converged settings.
    """

    print("")
    print("-----------------------")
    log.write_blue("Find suitable angular_resolution")

    angular_resolution = simulation.angular_resolution # NOTE: should we have a starting-value flag?
    print("Starting value: angular_resolution=%f rad" % angular_resolution)
    current_value = evaluate(simulation, detector)

    if ax is not None:
        x = np.array([])
        y = np.array([])
        r = np.array([])
        init_x = angular_resolution
        init_y = current_value.real
        line1, = ax[0].plot(np.insert(x, 0, init_x), np.insert(y.real, 0, init_y),'.-')
        line2, = ax[1].plot(x,r,'.-')

    for _ in range(max_iter):
        old_angular_resolution = angular_resolution
        angular_resolution = angular_resolution / 2.0
        simulation.angular_resolution = angular_resolution
        simulation.set_default_angles()

        print("---------------------------------------")
        print("Try angular_resolution = %f rad" % angular_resolution)

        new_value = evaluate(simulation, detector)

        rel_diff = abs(new_value - current_value) / abs(current_value)
        print("Old detector value:", current_value)
        print("New detector value:", new_value)
        print("Relative difference:", rel_diff)
        print("Allowed tolerance:  ", tolerance)

        if ax is not None:
            x = np.append(x, angular_resolution)
            y = np.append(y, new_value)
            r = np.append(r, rel_diff)
            line1.set_data(np.insert(x, 0, init_x), np.insert(y.real, 0, init_y.real))
            line1.set_label('$n_{{eff}}^{{max}} = {:g}$, $l_{{max}} = {}$, $m_{{max}} = {}$, $\delta n_{{eff}} = {}$'.\
                            format(simulation.neff_max.real,
                                   simulation.particle_list[0].l_max,
                                   simulation.particle_list[0].m_max,
                                   simulation.neff_resolution))
            line2.set_data(x, r)
            [( ax.relim(), ax.autoscale_view()) for ax in ax]
            plt.draw()
            plt.pause(0.001)
            ax[0].legend()

        if rel_diff < tolerance:  # in this case: discard halfed neff_resolution
            simulation.angular_resolution = old_angular_resolution
            simulation.set_default_angles()
            log.write_green("Relative difference smaller than tolerance. Keep angular_resolution = %g deg" % old_angular_resolution)
            if ax is not None:
                titlestr = "relative diff < {:g}, keep $\delta \\theta = {:g}$ deg".format(tolerance, old_angular_resolution)
                ax[1].title.set_text(titlestr)
                ax[1].title.set_color('g')
                plt.draw()
            return current_value
        else:
            titlestr = "no convergence achieved, keep $\delta \\theta = {:g}$ deg".format(angular_resolution)
            current_value = new_value

    log.write_red("No convergence achieved. Keep angular_resolution = %g deg" % angular_resolution)
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
                                max_iter=30,
                                neff_imag=1e-2,
                                neff_resolution=1e-2,
                                select_neff_max=True,
                                neff_max_increment=0.5,
                                neff_max_offset=0,
                                neff_max=None,
                                select_neff_resolution=True,
                                select_angular_resolution=False,
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
        tolerance_steps (int):                      Number of consecutive steps at which the tolerance must be met
                                                    during multipole truncation convergence. Default: 2
        max_iter (int):                             Break convergence loops after that number of iterations, even if
                                                    no convergence has been achieved.
        neff_imag (float):                          Extent of the contour into the negative imaginary direction
                                                    (in terms of effective refractive index, n_eff=kappa/omega)
        neff_resolution (float):                    Discretization of the contour (in terms of eff. refractive
                                                    index) - if `select_neff_resolution` is true, this value will be
                                                    eventually overwritten. However, it is required in any case.
                                                    Default: 1e-2
        select_neff_max (logical):                  If set to true (default), the Sommerfeld integral truncation
                                                    parameter `neff_max` is determined automatically with the help
                                                    of a Cauchy convergence criterion.
        neff_max_increment (float):                 Only needed if `select_neff_max` is true.
                                                    Step size with which `neff_max` is incremented.
        neff_max_offset (float):                    Only needed if `select_neff_max` is true.
                                                    Start n_eff selection from the last estimated singularity
                                                    location plus this value (in terms of effective refractive index)
        neff_max (float):                           Only needed if `select_neff_max` is false.
                                                    Truncation value of contour (in terms of effective refractive
                                                    index).
        select_neff_resolution (logical):           If set to true (default), the Sommerfeld integral discretization
                                                    parameter `neff_resolution` is determined automatically with the
                                                    help of a Cauchy convergence criterion.
        select_angular_resolution (logical):        If set to true, the angular resolution step for the default polar
                                                    and azimuthal angles is determined automatically according to a
                                                    Cauchy convergenge criterion.
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

    print("")
    print("----------------------------------------")
    print("Starting automatic paramateter selection")
    print("----------------------------------------")

    if tolerance_factor <= 0 or tolerance_factor > 1:
        raise ValueError('please provide 0 < tolerance_factor <= 1')

    if select_neff_max:
        if relative_convergence:
            if show_plot:
                plt.ion()
                _, ax_array = _init_fig('$l_{max}$', detector,
                                        'multipole cutoff',
                                        tolerance, tolerance_factor * tolerance,
                                        2)
            else:
                ax_array = None

            target_value = converge_neff_max(simulation=simulation,
                                             detector=detector,
                                             tolerance=tolerance,
                                             tolerance_factor=tolerance_factor,
                                             tolerance_steps=tolerance_steps,
                                             max_iter=max_iter,
                                             neff_imag=neff_imag,
                                             neff_resolution=neff_resolution,
                                             neff_max_increment=neff_max_increment,
                                             neff_max_offset=neff_max_offset,
                                             converge_lm=relative_convergence,
                                             ax=ax_array)
            neff_max = simulation.neff_max

            if select_multipole_cutoff:
                if show_plot:
                    plt.ion()
                    _, ax_array = _init_fig('multipole order', detector,
                                            'multipole cutoff selection',
                                            tolerance)
                else:
                    ax_array = None
                converge_multipole_cutoff(simulation=simulation,
                                          detector=detector,
                                          tolerance=tolerance,
                                          tolerance_steps=tolerance_steps,
                                          max_iter=max_iter,
                                          current_value=target_value,
                                          ax=ax_array)

        else:  # select neff_max but no relative convergence
            if select_multipole_cutoff:
                if show_plot:
                    plt.ion()
                    _, ax_array = _init_fig('multipole order', detector,
                                            'multipole cutoff selection',
                                            tolerance)
                else:
                    ax_array = None
                converge_multipole_cutoff(simulation=simulation,
                                          detector=detector,
                                          tolerance=tolerance,
                                          tolerance_steps=tolerance_steps,
                                          max_iter=max_iter,
                                          current_value=None,
                                          ax=ax_array)

            if show_plot:
                plt.ion()
                _, ax_array = _init_fig('$l_{max}$', detector,
                                        'multipole cutoff',
                                        tolerance, tolerance_factor * tolerance,
                                        2)
            else:
                ax_array = None

            target_value = converge_neff_max(simulation=simulation,
                                             detector=detector,
                                             tolerance=tolerance,
                                             tolerance_factor=tolerance_factor,
                                             tolerance_steps=tolerance_steps,
                                             max_iter=max_iter,
                                             neff_imag=neff_imag,
                                             neff_resolution=neff_resolution,
                                             neff_max_increment=neff_max_increment,
                                             neff_max_offset=neff_max_offset,
                                             converge_lm=relative_convergence,
                                             ax=ax_array)
            neff_max = simulation.neff_max

    else:  # no neff_max selection
        ax_array = None
        if select_multipole_cutoff:
            if show_plot:
                plt.ion()
                _, ax_array = _init_fig('multipole order', detector,
                                        'multipole cutoff selection',
                                        tolerance)
            else:
                ax_array = None
            converge_multipole_cutoff(simulation=simulation,
                                      detector=detector,
                                      tolerance=tolerance,
                                      tolerance_steps=tolerance_steps,
                                      max_iter=max_iter,
                                      current_value=None,
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

    if select_angular_resolution:
        if show_plot:
            plt.ion()
            _, ax_array = _init_fig('$\delta \\theta$ [deg]', detector, 'angular resolution selection', tolerance)
        else:
            ax_array = None
        converge_angular_resolution(simulation=simulation,
                                    detector=detector,
                                    tolerance=tolerance,
                                    max_iter=max_iter,
                                    ax=ax_array)

    if show_plot:
        plt.ioff()


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
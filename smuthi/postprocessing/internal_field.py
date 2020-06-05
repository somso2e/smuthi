"""Manage post processing steps to evaluate the electric field inside a sphere"""

import smuthi.fields as flds
import smuthi.fields.expansions as fldex
import smuthi.linearsystem.tmatrix.t_matrix as tmt


def internal_field_piecewise_expansion(vacuum_wavelength, particle_list, layer_system):
    """Compute a piecewise field expansion of the internal field of spheres.

    Args:
        vacuum_wavelength (float):                  vacuum wavelength
        particle_list (list):                       list of smuthi.particles.Particle objects
        layer_system (smuthi.layers.LayerSystem):   stratified medium

    Returns:
        internal field as smuthi.field_expansion.PiecewiseFieldExpansion object

    """
    intfld = fldex.PiecewiseFieldExpansion()

    for particle in particle_list:
        if type(particle).__name__ == 'Sphere':
            i_part = layer_system.layer_number(particle.position[2])
            k_medium = flds.angular_frequency(vacuum_wavelength) * layer_system.refractive_indices[i_part]
            k_particle = flds.angular_frequency(vacuum_wavelength) * particle.refractive_index

            internal_field = fldex.SphericalWaveExpansion(
                k_particle,
                l_max=particle.l_max,
                m_max=particle.m_max,
                kind='regular',
                reference_point=particle.position)

            internal_field.validity_conditions.append(particle.is_inside)

            for m in range(-particle.m_max, particle.m_max+1):
                for l in range(max(1, abs(m)), particle.l_max+1):
                    for tau in range(2):
                        n = flds.multi_to_single_index(tau, l, m, particle.l_max, particle.m_max)
                        b_to_c = (tmt.internal_mie_coefficient(tau, l, k_medium, k_particle, particle.radius)
                                  / tmt.mie_coefficient(tau, l, k_medium, k_particle, particle.radius))
                        internal_field.coefficients[n] = particle.scattered_field.coefficients[n] * b_to_c

            intfld.expansion_list.append(internal_field)
            particle.internal_field = internal_field

    return intfld

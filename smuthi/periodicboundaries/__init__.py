"""This subpackage contains functionality that has to do with light scattering
in periodic particle arrangements.
The __init__ module contains some helper functions (e.g. with respect to
periodic particle arrangements) and is the place to store the default separation
parameter that splits the Ewald lattice sum into real and reciprocal space."""

import numpy as np



"""Default Ewald sum separation parameter eta - needs to be set, e.g. at beginning of a periodic simulation"""
default_Ewald_sum_separation = None


def set_ewald_sum_separation(a1, a2, initial_field=None, particle_list=None, layer_system=None, magM=None):
    """ Set the separation parameter eta that splits the evaluation of the Ewald lattice sum
        into real and reciprocal space.
    Args:
        a1 (numpy.ndarray):                                 lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):                                 lattice vector 2 in Carthesian coordinates
        initial_field (smuthi.initial_field.PlaneWave):     initial plane wave object
        layer_system (smuthi.layers.LayerSystem):           stratified medium
        particle_list (list):                               list of smuthi.particles.Particle objects
        magM (float):                                       maximum tolerated magnitude with regard to
                                                            target accuracy
    Returns:
        Separation parameter eta 
    """
    if magM == None or magM == 1:
        return np.sqrt(np.pi / np.linalg.norm(np.cross(a1, a2)))
    else:
        z_list = [particle.position[2] for particle in particle_list]
        is_list = [layer_system.layer_number(z) for z in z_list]
        assert is_list.count(is_list[0]) == len(is_list)  # all particles in same layer?
        
        k = complex(2 * np.pi * layer_system.refractive_indices[is_list[0]] / initial_field.vacuum_wavelength)
        
        if initial_field.polar_angle < np.pi:
            pwe_exc = initial_field.plane_wave_expansion(layer_system, is_list[0])[0]
        else:
            pwe_exc = initial_field.plane_wave_expansion(layer_system, is_list[0])[1]
        k0t = np.array([pwe_exc.k_parallel[0] * np.cos(pwe_exc.azimuthal_angles)[0],
                        pwe_exc.k_parallel[0] * np.sin(pwe_exc.azimuthal_angles)[0]])     
        
        return (np.sqrt(k ** 2 - np.linalg.norm(k0t) ** 2) / (2 * np.log(magM))).real

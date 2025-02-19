import sys
import smuthi.initial_field as init
import smuthi.particles as part
import smuthi.simulation as simul
import smuthi.layers as lay
import smuthi.postprocessing.far_field as farf


ld = 550
rD1 = [100, -100, 130]
D1 = [1e7, 2e7, 3e7]
rD2 = [-100, 100, 70]
D2 = [3e7, -2e7, 1e7]
rD3 = [-100, 100, -100]
D3 = [-2e7, 3e7, 1e7]

# waypoints = [0, 0.8, 0.8-0.1j, 2.1-0.1j, 2.1, 3]
#neff_max = 3
#neff_discr = 1e-2
#neff_imag = 1e-2

#coord.set_default_k_parallel(vacuum_wavelength=ld, neff_imag=neff_imag, neff_resolution=neff_discr, neff_max=neff_max)

# initialize particle object
sphere1 = part.Sphere(position=[200, 200, 300], refractive_index=2.4 + 0.0j, radius=110, l_max=3, m_max=3)
sphere2 = part.Sphere(position=[-200, -200, 300], refractive_index=2.4 + 0.0j, radius=120, l_max=3, m_max=3)
sphere3 = part.Sphere(position=[-200, 200, 300], refractive_index=2.5 + 0.0j, radius=90, l_max=3, m_max=3)
part_list = [sphere1, sphere2, sphere3]

# initialize layer system object
lay_sys = lay.LayerSystem([0, 400, 0], [2, 1.3, 2])

# initialize dipole object
dipole1 = init.DipoleSource(vacuum_wavelength=ld, dipole_moment=D1, position=rD1)

dipole2 = init.DipoleSource(vacuum_wavelength=ld, dipole_moment=D2, position=rD2)

dipole3 = init.DipoleSource(vacuum_wavelength=ld, dipole_moment=D3, position=rD3)

dipole_collection = init.DipoleCollection(vacuum_wavelength=ld, 
                                          compute_swe_by_pwe=False, 
                                          compute_dissipated_power_by_pwe=False)
dipole_collection.append(dipole1)
dipole_collection.append(dipole2)
dipole_collection.append(dipole3)

# run simulation
simulation = simul.Simulation(layer_system=lay_sys, particle_list=part_list, initial_field=dipole_collection,
                              log_to_terminal=(not sys.argv[0].endswith('nose2')))  # suppress output if called by nose

simulation.run()

# evaluate power
power_list = simulation.initial_field.dissipated_power(particle_list=simulation.particle_list,
                                                       layer_system=simulation.layer_system)

power = sum(power_list)
ff_tup = farf.total_far_field(simulation.initial_field, simulation.particle_list, simulation.layer_system)


def test_energy_conservation():
    ff_power = sum(ff_tup[0].integral())
    err = abs((power - ff_power) / ff_power)
    print('ff power', ff_power)
    print('diss power list', power_list)
    print('diss power', power)
    print('error', err)
    assert err < 5e-4
    print("Test passed.")


dipole_collection2 = init.DipoleCollection(vacuum_wavelength=ld, 
                                          compute_swe_by_pwe=True, 
                                          compute_dissipated_power_by_pwe=True)
dipole_collection2.append(dipole1)
dipole_collection2.append(dipole2)
dipole_collection2.append(dipole3)

# evaluate power
power_list_alt = dipole_collection2.dissipated_power(particle_list=simulation.particle_list,
                                                     layer_system=simulation.layer_system)


def test_alternative_power():
    err = abs((sum(power_list) - sum(power_list_alt)) / sum(power_list))
    print('alt power error', err)
    assert err < 1e-4
    print("Test passed.")
        

if __name__ == '__main__':
    test_energy_conservation()
    test_alternative_power()

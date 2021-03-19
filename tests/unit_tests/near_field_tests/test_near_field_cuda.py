import sys
import smuthi.initial_field as init
import smuthi.particles as part
import smuthi.simulation as simul
import smuthi.layers as lay
import smuthi.postprocessing.scattered_field as sf
import smuthi.utility.cuda as cu
import numpy as np

import unittest


class TestNearFieldCUDA(unittest.TestCase):
    def test_fields_by_gpu_and_cpu(self):
        try:
            import pycuda.autoinit
        except ImportError:
            self.skipTest('PyCUDA is not available')

        ld = 550
        rD = [100, -100, 100]
        D = [1e7, 2e7, 3e7]
        #waypoints = [0, 0.8, 0.8 - 0.1j, 2.1 - 0.1j, 2.1, 4]
        #neff_discr = 2e-2

        #coord.set_default_k_parallel(vacuum_wavelength=ld, neff_waypoints=waypoints, neff_resolution=neff_discr)

        # initialize particle object
        sphere1 = part.Sphere(position=[200, 200, 300], refractive_index=2.4 + 0.0j, radius=110, l_max=3, m_max=3)
        sphere2 = part.Sphere(position=[-200, -200, 300], refractive_index=2.4 + 0.0j, radius=120, l_max=3, m_max=3)
        sphere3 = part.Sphere(position=[-200, 200, 300], refractive_index=2.5 + 0.0j, radius=90, l_max=3, m_max=3)
        part_list = [sphere1, sphere2, sphere3]

        # initialize layer system object
        lay_sys = lay.LayerSystem([0, 400, 0], [1 + 6j, 2.3, 1.5])

        # initialize dipole object
        dipole = init.DipoleSource(vacuum_wavelength=ld, dipole_moment=D, position=rD)

        # run simulation
        simulation = simul.Simulation(layer_system=lay_sys, particle_list=part_list, initial_field=dipole,
                                      log_to_terminal='nose2' not in sys.modules.keys())  # suppress output if called by nose
        simulation.run()

        xarr = np.array([-300, 400, -100, 200])
        yarr = np.array([200, -100, 400, 300])
        zarr = np.array([-50, 200, 600, 700])

        wavelength = 600

        scat_fld_exp = sf.scattered_field_piecewise_expansion(ld, part_list, lay_sys)
        e_x_scat_cpu, e_y_scat_cpu, e_z_scat_cpu = scat_fld_exp.electric_field(xarr, yarr, zarr)
        e_x_init_cpu, e_y_init_cpu, e_z_init_cpu = simulation.initial_field.electric_field(xarr, yarr, zarr, lay_sys)

        h_x_scat_cpu, h_y_scat_cpu, h_z_scat_cpu = scat_fld_exp.magnetic_field(xarr, yarr, zarr, wavelength)
        h_x_init_cpu, h_y_init_cpu, h_z_init_cpu = simulation.initial_field.magnetic_field(xarr, yarr, zarr, lay_sys)

        cu.enable_gpu()
        scat_fld_exp = sf.scattered_field_piecewise_expansion(ld, part_list, lay_sys)
        e_x_scat_gpu, e_y_scat_gpu, e_z_scat_gpu = scat_fld_exp.electric_field(xarr, yarr, zarr)
        e_x_init_gpu, e_y_init_gpu, e_z_init_gpu = simulation.initial_field.electric_field(xarr, yarr, zarr, lay_sys)
        
        h_x_scat_gpu, h_y_scat_gpu, h_z_scat_gpu = scat_fld_exp.magnetic_field(xarr, yarr, zarr, wavelength)
        h_x_init_gpu, h_y_init_gpu, h_z_init_gpu = simulation.initial_field.magnetic_field(xarr, yarr, zarr, lay_sys)

        cpu_electric_fields = (e_x_scat_cpu, e_y_scat_cpu, e_z_scat_cpu,
                                        e_x_init_cpu, e_y_init_cpu, e_z_init_cpu)
        gpu_electric_fields = (e_x_scat_gpu, e_y_scat_gpu, e_z_scat_gpu,
                                        e_x_init_gpu, e_y_init_gpu, e_z_init_gpu)

        cpu_magnetic_fields = (h_x_scat_cpu, h_y_scat_cpu, h_z_scat_cpu,
                                        h_x_init_cpu, h_y_init_cpu, h_z_init_cpu)
        gpu_magnetic_fields = (h_x_scat_gpu, h_y_scat_gpu, h_z_scat_gpu,
                                        h_x_init_gpu, h_y_init_gpu, h_z_init_gpu)

        assert_electric_fields_computed_by_cpu_and_gpu(cpu_electric_fields, gpu_electric_fields)
        assert_magnetic_fields_computed_by_cpu_and_gpu(cpu_magnetic_fields, gpu_magnetic_fields)


def assert_arrays_in_tuples_are_close(first_tuple_with_arrays, second_tuple_with_arrays, rtol=1e-5):
    if len(first_tuple_with_arrays) != len(second_tuple_with_arrays):
        raise Exception("tuples must contain the same number of elements!")

    for i in range(len(first_tuple_with_arrays)):
        np.testing.assert_allclose(np.linalg.norm(first_tuple_with_arrays[i]),
                                    np.linalg.norm(second_tuple_with_arrays[i]),
                                    rtol=rtol)


# we declare these wrappers just to have more clear traceback.
def assert_electric_fields_computed_by_cpu_and_gpu(cpu_electric_fields,
                                                            gpu_electric_fields):
    assert_arrays_in_tuples_are_close(cpu_electric_fields, gpu_electric_fields)


def assert_magnetic_fields_computed_by_cpu_and_gpu(cpu_magnetic_fileds,
                                                            gpu_magnetic_fileds):
    assert_arrays_in_tuples_are_close(cpu_magnetic_fileds, gpu_magnetic_fileds)


if __name__ == '__main__':
    unittest.main()
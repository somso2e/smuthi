import numpy as np
import smuthi.utility.numba_helpers as nh

import unittest

class TestNumba(unittest.TestCase):
    def test_numba_trapz(self):
        test_3d_array = get_complex_random_array((90, 100, 110)).astype(np.complex64)
        test_1d_array = np.random.rand(110)

        default_result = np.trapz(test_3d_array, test_1d_array)
        custom_result = nh.numba_trapz_3dim_array(test_3d_array, test_1d_array)

        np.testing.assert_allclose(default_result, custom_result)


    def test_numba_tensordot(self):
        np.random.seed(1)
        test_2d_x_array = get_complex_random_array((100, 110), seed=1).astype(np.complex64)
        test_2d_y_array = get_complex_random_array((100, 110), seed=2).astype(np.complex64)
        test_2d_z_array = get_complex_random_array((100, 110), seed=3).astype(np.complex64)
        test_1d_x_array = np.random.rand(50).astype(np.float32)
        test_1d_y_array = np.random.rand(50).astype(np.float32)
        test_1d_z_array = np.random.rand(50).astype(np.float32)

        default_result = np.tensordot(test_1d_x_array.copy(), test_2d_x_array.copy(), axes=0)
        default_result += np.tensordot(test_1d_y_array.copy(), test_2d_y_array.copy(), axes=0)
        default_result += np.tensordot(test_1d_z_array.copy(), test_2d_z_array.copy(), axes=0)
        custom_result = nh.numba_3tensordots_1dim_times_2dim( \
                                test_1d_x_array, test_1d_y_array, test_1d_z_array, \
                                test_2d_x_array, test_2d_y_array, test_2d_z_array)

        np.testing.assert_allclose(default_result, custom_result, atol=1e-6, rtol=1e-6)


    def test_evaluate_r_times_eikr(self):
        np.random.seed(1)
        foo_x = get_complex_random_array((110, 120), seed=1).astype(np.complex64)[None, :, :]
        foo_y = get_complex_random_array((110, 120), seed=2).astype(np.complex64)[None, :, :]
        foo_z = get_complex_random_array((110, 120), seed=3).astype(np.complex64)[None, :, :]
        kr = get_complex_random_array((100, 110, 120), seed=4).astype(np.complex64)

        eikr = np.exp(1j * kr)
        default_foo_x_eikr = foo_x * eikr
        default_foo_y_eikr = foo_y * eikr
        default_foo_z_eikr = foo_z * eikr

        numba_foo_x_eikr, numba_foo_y_eikr, numba_foo_z_eikr = \
                nh.evaluate_r_times_eikr(foo_x, foo_y, foo_z, kr)

        np.testing.assert_allclose(default_foo_x_eikr, numba_foo_x_eikr, atol=1e-6, rtol=1e-6)
        np.testing.assert_allclose(default_foo_y_eikr, numba_foo_y_eikr, atol=1e-6, rtol=1e-6)
        np.testing.assert_allclose(default_foo_z_eikr, numba_foo_z_eikr, atol=1e-6, rtol=1e-6)


def get_complex_random_array(dimensions, seed=0):
    np.random.seed(seed)
    real_part = np.random.rand(*dimensions)
    imaginary_part = np.random.rand(*dimensions)
    complex_array = real_part + imaginary_part * 1.0j
    return complex_array

if __name__ == '__main__':
    unittest.main()
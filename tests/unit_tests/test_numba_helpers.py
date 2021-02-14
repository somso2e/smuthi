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


def get_complex_random_array(dimensions, seed=0):
    np.random.seed(seed)
    real_part = np.random.rand(*dimensions)
    imaginary_part = np.random.rand(*dimensions)
    complex_array = real_part + imaginary_part * 1.0j
    return complex_array

if __name__ == '__main__':
    unittest.main()
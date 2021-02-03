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


def get_complex_random_array(dimensions):
    real_part = np.random.rand(*dimensions)
    imaginary_part = np.random.rand(*dimensions)
    complex_array = real_part + imaginary_part * 1.0j
    return complex_array

if __name__ == '__main__':
    unittest.main()
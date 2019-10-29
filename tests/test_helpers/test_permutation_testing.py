"""This module contains classes and functions for testing permutation testing.

Note:
    * Uses pytest style testing.
"""
import numpy as np
import pytest

from lmpy import Matrix

import analyses.helpers.permutation_testing as perm_testing


# .............................................................................
class Test_compare_methods(object):
    """Tests the compare methods for absolute and signed values.
    """
    # .....................................
    def test_compare_absolute_values(self):
        """Tests that compare_absolute_values works how it is expected.
        """
        # Tests for greater than absolute values
        assert perm_testing.compare_absolute_values(1, 2)
        assert perm_testing.compare_absolute_values(1, -2)
        assert perm_testing.compare_absolute_values(-1, 2)
        assert perm_testing.compare_absolute_values(-1, -2)

        # Tests for less than absolute values
        assert not perm_testing.compare_absolute_values(2, 1)
        assert not perm_testing.compare_absolute_values(2, -1)
        assert not perm_testing.compare_absolute_values(-2, 1)
        assert not perm_testing.compare_absolute_values(-2, -1)

        # Tests for equal absolute values
        assert not perm_testing.compare_absolute_values(2, 2)
        assert not perm_testing.compare_absolute_values(2, -2)
        assert not perm_testing.compare_absolute_values(-2, 2)
        assert not perm_testing.compare_absolute_values(-2, -2)

    # .....................................
    def test_compare_signed_values(self):
        """Tests that compare_signed_values works how it is expected.
        """
        # Tests for greater than signed values
        assert perm_testing.compare_signed_values(1, 2)
        assert perm_testing.compare_signed_values(-1, 2)
        assert perm_testing.compare_signed_values(-2, 2)
        assert perm_testing.compare_signed_values(-2, 1)
        assert perm_testing.compare_signed_values(-2, -1)

        # Tests for less than signed values
        assert not perm_testing.compare_signed_values(2, 1)
        assert not perm_testing.compare_signed_values(2, -1)
        assert not perm_testing.compare_signed_values(2, 2)
        assert not perm_testing.compare_signed_values(2, -2)
        assert not perm_testing.compare_signed_values(-2, -2)
        assert not perm_testing.compare_signed_values(1, -2)
        assert not perm_testing.compare_signed_values(-1, -2)


# .............................................................................
class Test_correct_p_values(object):
    """Tests the correct_p_values method.
    """
    # .....................................
    def test_valid(self):
        """Tests that correcting p-values does what is expected.
        """
        uncorrected = Matrix(
            np.array([
                [0.05, 0.1, 0.02], [0.01, 0.05, 0.06], [0.1, 0.01, 0.20]]))
        corrected = perm_testing.correct_p_values(uncorrected)


# .............................................................................
class Test_get_p_values(object):
    """Tests the get_p_values method.
    """
    # .....................................
    def test_with_absolute_value_comparison(self):
        """Tests that getting p-values does what is expected.
        """
        obs_matrix = Matrix(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
        rand_1 = Matrix(np.array([[3, 2, 1], [6, 3, -12], [8, 3, -10]]))
        rand_2 = Matrix(np.array([[9, 23, 1], [4, 2, 9], [-32, -3, 9]]))
        p_vals = perm_testing.get_p_values(
            obs_matrix, [rand_1, rand_2],
            compare_func=perm_testing.compare_absolute_values)
        assert np.all(
            p_vals[:, :, 0] == np.array(
                [[1, 0.5, 0], [0.5, 0, 1], [1, 0, 0.5]]))

    # .....................................
    def test_with_signed_value_comparison(self):
        """Tests that getting p-values does what is expected.
        """
        obs_matrix = Matrix(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
        rand_1 = Matrix(np.array([[3, 2, 1], [6, 3, -12], [8, 3, -10]]))
        rand_2 = Matrix(np.array([[9, 23, 1], [4, 2, 9], [-32, -3, 9]]))
        p_vals = perm_testing.get_p_values(
            obs_matrix, [rand_1, rand_2],
            compare_func=perm_testing.compare_signed_values)
        assert np.all(
            p_vals[:, :, 0] == np.array(
                [[1, 0.5, 0], [0.5, 0, 0.5], [0.5, 0, 0]]))

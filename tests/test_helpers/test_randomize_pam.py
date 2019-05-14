"""This module contains classes and functions for testing pam randomization.

Note:
    * Uses pytest style testing.
"""
import numpy as np
import pytest

from analyses.helpers.randomize_pam import swap_randomize
from analyses.lm_objects.matrix import Matrix


# .............................................................................
def get_random_matrix(*dim_size):
    """Generates a randomized matrix with the shape provided in dim_size.

    Args:
        *dim_size (:obj:`list` of :obj:`int`): Variable length argument list of
            integers representing matrix dimension sizes.
    """
    headers = {}
    i = 0
    for i in range(len(dim_size)):
        headers[str(i)] = [
            'header-{}-{}'.format(i, x) for x in range(dim_size[i])]
    return Matrix(np.random.randint(
        low=0, high=2, size=tuple(dim_size)), headers=headers)


# .............................................................................
class Test_randomize_pam(object):
    """Test class for the swap_randomize method.
    """
    # .....................................
    def test_valid(self):
        """Tests that the function works as expected.
        """
        mtx = get_random_matrix(5, 5)
        rand_mtx = swap_randomize(mtx, 100)

        assert np.sum(mtx.data) == np.sum(rand_mtx.data)

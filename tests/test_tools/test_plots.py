"""This module contains classes and functions for testing plots

Notes:
    * Uses pytest style testing
"""
import os

import numpy as np
import pytest

from analyses.lm_objects.matrix import Matrix
from analyses.lm_objects.tree import TreeWrapper
from analyses.tools.plots import create_distribution_plots


# .............................................................................
class Test_create_distribution_plots(object):
    """Test class for the create_distribution_plots method
    """
    # .....................................
    def test_valid(self, tmpdir):
        """Test the function with valid inputs
        """
        # Create a tree
        tree = TreeWrapper.get(data='(A,(B,((C,D),(E,F))));', schema='newick')
        mtx = Matrix(
            np.random.random((6, 3, 2)),
            headers={'0': ['A', 'B', 'C', 'D', 'E', 'F'],
                     '1': ['label', 'other_val', 'one_more_val']})
        # This should not fail
        output_directory = os.path.join(tmpdir.dirname, 'plots')
        create_distribution_plots(tree, mtx, output_directory)

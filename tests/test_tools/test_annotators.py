"""This module contains classes and functions for testing annotators.

Note:
    * Uses pytest style testing.
"""
import numpy as np
import pytest

from lmpy import Matrix, TreeWrapper

import analyses.tools.annotators as annotators


# .............................................................................
class Test_add_all_annotations(object):
    """Test class for the add_all_annotations method.
    """
    # .....................................
    def test_valid(self):
        """Test the function with valid inputs.
        """
        # Create a tree
        tree = TreeWrapper.get(data='(A,(B,((C,D),(E,F))));', schema='newick')
        mtx = Matrix(
            np.random.random((6, 3, 1)),
            headers={'0': ['A', 'B', 'C', 'D', 'E', 'F'],
                     '1': ['label', 'other_val', 'one_more_val']})
        # This should not fail
        annotators.add_all_annotations(tree, mtx, update=True)


# .............................................................................
class Test_annotate_tree_with_label(object):
    """Test class for the annotate_tree_with_label method.
    """
    # .....................................
    def test_valid(self):
        """Test the function with valid inputs.
        """
        # Create a tree
        tree = TreeWrapper.get(data='(A,(B,((C,D),(E,F))));', schema='newick')
        mtx = Matrix(
            np.random.random((6, 2, 1)),
            headers={'0': ['A', 'B', 'C', 'D', 'E', 'F'],
                     '1': ['label', 'other_val']})
        # This should not fail
        annotators.annotate_tree_with_label(tree, mtx, label_column=0)

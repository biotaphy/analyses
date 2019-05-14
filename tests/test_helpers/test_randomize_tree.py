"""This module contains classes and functions for testing tree randomization.

Note:
    * Uses pytest style testing.
"""
import pytest

from analyses.helpers.randomize_tree import randomize_tree
import analyses.lm_objects.tree as tree


# .............................................................................
class Test_randomize_tree(object):
    """test class for randomize_tree method.
    """
    # .....................................
    def test_valid_newick_tree(self, valid_newick_tree):
        """Attempt to randomize a newick tree from file.

        Args:
            valid_newick_tree (pytest.fixture): A parameterized pytest fixture
                that provides valid newick trees, one at a time, to this test
                function.
        """
        loaded_tree = tree.TreeWrapper.from_filename(valid_newick_tree)
        assert isinstance(loaded_tree, tree.TreeWrapper)
        randomized_tree = randomize_tree(loaded_tree)

        for tax in randomized_tree.taxon_namespace:
            assert tax in loaded_tree.taxon_namespace

    # .....................................
    def test_valid_nexml_tree(self, valid_nexml_tree):
        """Attempt to randomize a nexml tree from file.

        Args:
            valid_nexml_tree (pytest.fixture): A parameterized pytest fixture
                that provides valid nexml trees, one at a time, to this test
                function.
        """
        loaded_tree = tree.TreeWrapper.from_filename(valid_nexml_tree)
        assert isinstance(loaded_tree, tree.TreeWrapper)
        randomized_tree = randomize_tree(loaded_tree)

        for tax in randomized_tree.taxon_namespace:
            assert tax in loaded_tree.taxon_namespace

    # .....................................
    def test_valid_nexus_tree(self, valid_nexus_tree):
        """Attempt to randomize a nexus tree from file.

        Args:
            valid_nexus_tree (pytest.fixture): A parameterized pytest fixture
                that provides valid nexus trees, one at a time, to this test
                function.
        """
        loaded_tree = tree.TreeWrapper.from_filename(valid_nexus_tree)
        assert isinstance(loaded_tree, tree.TreeWrapper)
        randomized_tree = randomize_tree(loaded_tree)

        for tax in randomized_tree.taxon_namespace:
            assert tax in loaded_tree.taxon_namespace

"""This module tests the ancestral_reconstruction/lm_objects/tree.py module

The functions in this module are pytest style tests for the tree.py module
"""
import dendropy
import pytest

from ancestral_reconstruction.lm_objects import tree


# .............................................................................
class Test_PhyloTreeKeys(object):
    """Test the PhyloTreeKeys class

    This is a simple class that really just contains constants
    """
    # .....................................
    def test_get_constants(self):
        """Test that the constants can be retrieved
        """
        assert tree.PhyloTreeKeys.MTX_IDX is not None
        assert tree.PhyloTreeKeys.SQUID is not None


# .............................................................................
class Test_LmTreeException(object):
    """Test the LmException class
    """
    # .....................................
    def test_exception(self):
        """Attempt to throw the exception
        """
        with pytest.raises(tree.LmTreeException):
            raise tree.LmTreeException('Test exception')


# .............................................................................
class Test_TreeWrapper(object):
    """Test the TreeWrapper class

    The TreeWrapper class is an extension of dendropy.Tree.  Test that the
    functions do not break the existing functionality of the class and that
    they work properly
    """
    # .....................................
    def test_from_base_tree(self, data_files):
        """Attempt to get a tree using Dendropy and then wrap it

        Try to retrieve a tree using Dendropy and then send it through the
        wrapper function to determine if it produces a correctly wrapped tree
        """
        schemas = ['newick', 'nexus']
        for schema in schemas:
            for tree_filename in data_files.get_trees(schema, True):
                dendropy_tree = dendropy.Tree.get(path=tree_filename,
                                                  schema=schema)
                wrapped_tree = tree.TreeWrapper.from_base_tree(dendropy_tree)
                assert isinstance(wrapped_tree, tree.TreeWrapper)

    # .....................................
    def test_add_node_labels_no_prefix_no_overwrite(self):
        """Test that node labels are added correctly to a tree

        Attempt to add node labels without a prefix.  Any existing labels
        should be retained
        """
        newick_string = '(A,((B,C)testnode,(G,(D,(E,F)))));'
        existing_labels = ['A', 'B', 'C', 'testnode', 'G', 'D', 'E', 'F']
        my_tree = tree.TreeWrapper.get(data=newick_string, schema='newick')

        # Add node labels - Should all be integers
        my_tree.add_node_labels()
        print(my_tree.as_string(schema='newick'))

        for node in my_tree.nodes():
            if node.label is not None:
                label = node.label
            else:
                label = node.taxon.label
            # Check if the node label is in the existing labels
            if label in existing_labels:
                # Remove the label so we can verify that they were maintained
                existing_labels.pop(existing_labels.index(label))
            else:
                # Check that the label is the correct format
                assert str(int(label)) == label

        # Now verify that all of the previous labels were maintained
        assert len(existing_labels) == 0

    # .....................................
    def test_add_node_labels_with_prefix_no_overwrite(self):
        """Test that node labels are added correctly to a tree

        Attempt to add node labels with a prefix.  Any existing labels should
        be retained
        """
        node_prefix = 'nd_'
        newick_string = '(A,((B,C)testnode,(G,(D,(E,F)))));'
        existing_labels = ['A', 'B', 'C', 'testnode', 'G', 'D', 'E', 'F']
        my_tree = tree.TreeWrapper.get(data=newick_string, schema='newick')

        # Add node labels - Should all be integers
        my_tree.add_node_labels(prefix=node_prefix)
        print(my_tree.as_string(schema='newick'))

        for node in my_tree.nodes():
            if node.label is not None:
                label = node.label
            else:
                label = node.taxon.label
            # Check if the node label is in the existing labels
            if label in existing_labels:
                # Remove the label so we can verify that they were maintained
                existing_labels.pop(existing_labels.index(label))
            else:
                # Check that the label is the correct format
                assert label.startswith(node_prefix)
                # Check that we can create an integer from the suffix
                _ = int(label.strip(node_prefix))

        # Now verify that all of the previous labels were maintained
        assert len(existing_labels) == 0

    # .....................................
    def test_add_node_labels_no_prefix_yes_overwrite(self):
        """Test that node labels are added correctly to a tree

        Attempt to add node labels without a prefix.  Any existing node labels
        should not be retained
        """
        newick_string = '(A,((B,C)testnode,(G,(D,(E,F)))));'
        existing_labels = ['A', 'B', 'C', 'testnode', 'G', 'D', 'E', 'F']
        node_labels = ['testnode']  # Should not exist in modified version
        my_tree = tree.TreeWrapper.get(data=newick_string, schema='newick')

        # Add node labels - Should all be integers
        my_tree.add_node_labels(overwrite=True)
        print(my_tree.as_string(schema='newick'))

        for node in my_tree.nodes():
            if node.label is not None:
                label = node.label
            else:
                label = node.taxon.label
            # Check if the node label is in the existing labels
            if label in existing_labels:
                # Remove the label so we can verify that they were maintained
                existing_labels.pop(existing_labels.index(label))
            else:
                # Check that the label is the correct format
                assert str(int(label)) == label

        # Now verify that existing node labels were changed
        assert len(existing_labels) == len(node_labels)
        for node_label in existing_labels:
            assert node_label in node_labels

        for node_label in node_labels:
            assert node_label in existing_labels

    # .....................................
    def test_add_node_labels_with_prefix_yes_overwrite(self):
        """Test that node labels are added correctly to a tree

        Attempt to add node labels with a prefix.  Any existing labels should
        not be retained
        """
        node_prefix = 'nd_'
        newick_string = '(A,((B,C)testnode,(G,(D,(E,F)))));'
        existing_labels = ['A', 'B', 'C', 'testnode', 'G', 'D', 'E', 'F']
        node_labels = ['testnode']  # Should not exist in modified version
        my_tree = tree.TreeWrapper.get(data=newick_string, schema='newick')

        # Add node labels - Should all be integers
        my_tree.add_node_labels(prefix=node_prefix, overwrite=True)
        print(my_tree.as_string(schema='newick'))

        for node in my_tree.nodes():
            if node.label is not None:
                label = node.label
            else:
                label = node.taxon.label
            # Check if the node label is in the existing labels
            if label in existing_labels:
                # Remove the label so we can verify that they were maintained
                existing_labels.pop(existing_labels.index(label))
            else:
                # Check that the label is the correct format
                assert label.startswith(node_prefix)
                # Check that we can create an integer from the suffix
                _ = int(label.strip(node_prefix))

        # Now verify that existing node labels were changed
        assert len(existing_labels) == len(node_labels)
        for node_label in existing_labels:
            assert node_label in node_labels

        for node_label in node_labels:
            assert node_label in existing_labels

"""
@summary: This module contains classes and functions for testing ancestral
            state reconstruction
@note: Uses pytest style testing
"""
import dendropy
import os
import pytest

from ancestral_reconstruction.analysis import anc_dp
import ancestral_reconstruction.helpers.data_readers as data_readers


# .............................................................................
def test_calculate_continuous_ancestral_states(data_files):
    """
    @summary: Tests the calculate_continusous_ancestral_states method
    @note: This test will need to evolve as the output format changes.  It
            will probably be better to return a data structure with various
            values for each node rather than assigning the value to the node
            label
    """
    # Get the data files
    packages = data_files.get_packages()
    assert len(packages) > 0
    for tree_filename, alignment_filename, results_filename in packages:
        # Process the tree file
        _, tree_ext = os.path.splitext(tree_filename)
        if tree_ext == '.nex':
            tree_schema = 'nexus'
        elif tree_ext == '.xml':
            tree_schema = 'nexml'
        elif tree_ext == '.tre':
            tree_schema = 'newick'
        else:
            raise Exception(
                'Cannot handle tree with extension: {}'.format(tree_ext))
        tree = dendropy.Tree.get(path=tree_filename, schema=tree_schema)

        # Process the alignment file
        _, align_ext = os.path.splitext(alignment_filename)
        if align_ext == '.csv':
            with open(alignment_filename) as align_file:
                sequences, headers = data_readers.read_csv_alignment_flo(
                                                                align_file)
        elif align_ext == '.json':
            with open(alignment_filename) as align_file:
                sequences, headers = data_readers.read_json_alignment_flo(
                                                                align_file)
        elif align_ext == '.phylip':
            with open(alignment_filename) as align_file:
                sequences = data_readers.read_phylip_alignment_flo(
                                                                align_file)
        elif align_ext == '.tbl':
            with open(alignment_filename) as align_file:
                sequences = data_readers.read_table_alignment_flo(
                                                                align_file)
        else:
            raise Exception(
                'Cannot handle alignments with extension: {}'.format(
                    align_ext))

        # Run analysis
        anc_dp.calculate_continuous_ancestral_states(tree, sequences)

        # Check the output
        _, out_tree_ext = os.path.splitext(results_filename)
        if out_tree_ext == '.nex':
            out_tree_schema = 'nexus'
        elif out_tree_ext == '.xml':
            out_tree_schema = 'nexml'
        elif out_tree_ext == '.tre':
            out_tree_schema = 'newick'
        else:
            raise Exception(
                'Cannot handle tree with extension: {}'.format(out_tree_ext))
        out_tree_string = tree.as_string(schema=out_tree_schema)

        # Compare with results
        with open(results_filename) as results_file:
            results_string = results_file.read()

        assert out_tree_string.strip() == results_string.strip()

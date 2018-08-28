"""
@summary: This module contains classes and functions for testing ancestral
            state reconstruction
@note: Uses pytest style testing
"""
import dendropy
import numpy as np
import os
import pytest

from ancestral_reconstruction.analysis import anc_dp
import ancestral_reconstruction.helpers.data_readers as data_readers
from ancestral_reconstruction.lm_objects.tree import TreeWrapper


# .............................................................................
def test_calculate_continuous_ancestral_states_invalid(data_files):
    """
    @summary: Tests the calculate_continusous_ancestral_states method with
                invalid data
    @note: This test will need to evolve as the output format changes.  It
            will probably be better to return a data structure with various
            values for each node rather than assigning the value to the node
            label
    """
    # Get the data files
    packages = data_files.get_packages(False)
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

        char_mtx = data_readers.get_character_matrix_from_sequences_list(
                    sequences)
        # Run analysis
        with pytest.raises(Exception):
            anc_dp.calculate_continuous_ancestral_states(tree, char_mtx)


# .............................................................................
def test_calculate_continuous_ancestral_states_valid(data_files):
    """
    @summary: Tests the calculate_continusous_ancestral_states method
    @note: This test will need to evolve as the output format changes.  It
            will probably be better to return a data structure with various
            values for each node rather than assigning the value to the node
            label
    """
    # Get the data files
    packages = data_files.get_packages(True)
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
        # tree = dendropy.Tree.get(path=tree_filename, schema=tree_schema)
        tree = TreeWrapper.get(path=tree_filename, schema=tree_schema)

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

        char_mtx = data_readers.get_character_matrix_from_sequences_list(
                    sequences)
        # Run analysis
        # tree.add_node_labels()
        _, anc_mtx = anc_dp.calculate_continuous_ancestral_states(tree,
                                                                  char_mtx)

        # New testing method
        # (For now) assume that results file is csv with row headers for node
        #    labels and column headers for variables
        results = []
        h = None
        with open(results_filename) as results_file:
            for line in results_file:
                if h is None:
                    # Get headers
                    h = line.strip().split(',')[1:]
                else:
                    # Add result (without label) to list
                    node_result = [
                        float(i) for i in line.strip().split(',')[1:]]
                    results.append(np.array(node_result, dtype=float))

        # Look for all results (only maximum likelihood)
        for row in anc_mtx.data[:, :, 0]:
            found = False
            for i in range(len(results)):
                # Allow for some wiggle room with decimal precision
                if np.all(np.isclose(row, results[i])):
                    found = True
                    results.pop(i)
                    break
            if not found:
                raise Exception(
                    'Could not find expected result: {} in results'.format(
                        row))

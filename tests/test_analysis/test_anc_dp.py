"""This module is used for testing ancestral state reconstruction

Notes:
    * Uses pytest style testing
"""
import os

import dendropy
import numpy as np
import pytest

from analyses.ancestral_state import anc_dp
import analyses.helpers.data_readers as data_readers
from analyses.lm_objects.tree import TreeWrapper


# .............................................................................
class Test_calculate_continuous_ancestral_states(object):
    """Tests ancestral state reconstruction
    """
    # .....................................
    def test_package_invalid(self, data_files):
        """Test calculate_continusous_ancestral_states method with invalid data

        Args:
            data_files : A pytest fixture defined in conftest.py

        Note:
            * This test will need to evolve as the output format changes.  It
                will probably be better to return a data structure with various
                values for each node rather than assigning the value to the
                node label
        """
        # Get the data files
        packages = data_files.get_ancestral_state_packages(False)
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

    # .....................................
    def test_package_valid(self, data_files):
        """Tests the calculate_continusous_ancestral_states method

        Args:
            data_files : A pytest fixture defined in conftest.py

        Note:
            * This test will need to evolve as the output format changes.  It
                will probably be better to return a data structure with various
                values for each node rather than assigning the value to the
                node label
        """
        # Get the data files
        packages = data_files.get_ancestral_state_packages(True)
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
            _, anc_mtx = anc_dp.calculate_continuous_ancestral_states(
                tree, char_mtx, calc_std_err=True, sum_to_one=False)

            # New testing method
            # (For now) assume that results file is csv with row headers for
            #    node labels and column headers for variables
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


# .............................................................................
class Test_ancestal_distribution(object):
    """Tests ancestral distribution reconstruction
    """
    # .....................................
    def test_package_invalid(self, data_files):
        """Test calculate_ancestral_distributions method with invalid data

        Args:
            data_files : A pytest fixture defined in conftest.py

        Note:
            * This test will need to evolve as the output format changes.  It
                will probably be better to return a data structure with various
                values for each node rather than assigning the value to the
                node label
        """
        # Get the data files
        packages = data_files.get_ancestral_distribution_packages(False)
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
                anc_dp.calculate_ancestral_distributions(tree, char_mtx)

    # .....................................
    def test_package_valid(self, data_files):
        """Tests the calculate_ancestral_distributions method

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        # Get the data files
        packages = data_files.get_ancestral_distribution_packages(True)
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
            _, anc_mtx = anc_dp.calculate_ancestral_distributions(tree,
                                                                  char_mtx)

            # Testing method
            # Assume that the results file is a csv with row headers for node
            #    labels and output layer (maximum_likeliehood / standard_error)
            #    and column headers for variables
            ml_results = []
            std_err_results = []
            h = None
            with open(results_filename) as results_file:
                for line in results_file:
                    if h is None:
                        # Get headers
                        h = line.strip().split(',')[1:]
                    else:
                        # Add result (without label) to appropriate list
                        parts = line.strip().split(',')
                        layer = parts[1].lower()
                        values = np.array([float(i) for i in parts[2:]],
                                          dtype=float)
                        if layer == 'maximum_likelihood':
                            ml_results.append(values)
                        else:
                            std_err_results.append(values)
            assert(len(ml_results) == len(std_err_results))
            print('ml results')
            print(ml_results)
            print('std err results')
            print(std_err_results)

            # Look for all results (ml and std err results should match rows)
            for row_idx in range(anc_mtx.data.shape[0]):
                found = False
                # Get rows from data
                ml_row = anc_mtx.data[row_idx, :, 0]
                std_err_row = anc_mtx.data[row_idx, :, 1]

                for i in range(len(ml_results)):
                    print(ml_results[i])
                    print(std_err_results[i])
                    if np.all(np.isclose(ml_row, ml_results[i])) and \
                            np.all(np.isclose(std_err_row,
                                              std_err_results[i])):
                        found = True
                        ml_results.pop(i)
                        std_err_results.pop(i)
                        break
                if not found:
                    raise Exception(
                        'Could not find {}, {} in results'.format(ml_row,
                                                                  std_err_row))

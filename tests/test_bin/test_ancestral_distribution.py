"""This module tests the ancestral_distribution executable script
"""
import numpy as np
import os
import pytest
import subprocess
import sys


# .............................................................................
class Test_ancestral_distribution(object):
    """Test ancestral distribution reconstruction
    """
    base_dir = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), '../..')
    script_path = os.path.join(base_dir, 'bin/ancestral_distribution.py')

    # .....................................
    def test_package_valid(self, valid_ancestral_distribution_package, tmpdir):
        """Tests the calculate_ancestral_distributions method.

        Args:
            valid_ancestral_distribution_package (pytest.fixture): A pytest
                fixture that is parameterized to provide a valid ancestral
                distribution package, one at a time, and create new test
                functions for each.
            tmpdir (pytest.fixture): A built-in pytest fixture providing a
                temporary directory for this test function.

        Raises:
            IOError: When the tree or alignment cannot be loaded for the
                specified file extension.
            Exception: When a specified successful result value cannot be
                found.
        """
        # Get the data files
        (tree_filename, alignment_filename, results_filename
         ) = valid_ancestral_distribution_package
        # Process the tree file
        _, tree_ext = os.path.splitext(tree_filename)
        if tree_ext == '.nex':
            tree_schema = 'nexus'
        elif tree_ext == '.xml':
            tree_schema = 'nexml'
        elif tree_ext == '.tre':
            tree_schema = 'newick'
        else:
            raise IOError(
                'Cannot handle tree with extension: {}'.format(tree_ext))

        # Process alignment file
        _, align_ext = os.path.splitext(alignment_filename)
        if align_ext == '.csv':
            alignment_format = 'csv'
        elif align_ext == '.json':
            alignment_format = 'json'
        elif align_ext == '.phylip':
            alignment_format = 'phylip'
        elif align_ext == '.tbl':
            alignment_format = 'table'
        else:
            raise IOError(
                'Cannot handle alignments with extension: {}'.format(
                    align_ext))

        csv_filename = os.path.join(tmpdir.dirname, 'test_out.csv')
        out_tree_filename = os.path.join(tmpdir.dirname, 'test_out.nex')
        print(sys.path)
        cmd = '{} -c {} {} {} {} {} {} nexus'.format(
            self.script_path, csv_filename, tree_filename,
            tree_schema, alignment_filename, alignment_format,
            out_tree_filename)

        # Call process
        cmd2 = 'export PYTHONPATH={}; {}'.format(self.base_dir, cmd)
        res = subprocess.check_call(cmd2, shell=True)

        assert res == 0

        # Load output matrix
        out_ml_results = []
        out_std_err_results = []
        h = None
        with open(csv_filename) as output_file:
            for line in output_file:
                if h is None:
                    # Get headers
                    h = line.strip().split(',')[1:]
                else:
                    # Add result (without label) to appropriate list
                    parts = line.strip().split(',')
                    layer = parts[1].lower()
                    values = np.array(
                        [float(i) for i in parts[2:]], dtype=np.float)
                    if layer == 'maximum_likelihood':
                        out_ml_results.append(values)
                    else:
                        out_std_err_results.append(values)
        assert(len(out_ml_results) == len(out_std_err_results))

        # Load desired result matrix
        test_ml_results = []
        test_std_err_results = []
        h = None
        with open(csv_filename) as output_file:
            for line in output_file:
                if h is None:
                    # Get headers
                    h = line.strip().split(',')[1:]
                else:
                    # Add result (without label) to appropriate list
                    parts = line.strip().split(',')
                    layer = parts[1].lower()
                    values = np.array(
                        [float(i) for i in parts[2:]], dtype=np.float)
                    if layer == 'maximum_likelihood':
                        test_ml_results.append(values)
                    else:
                        test_std_err_results.append(values)
        assert(len(test_ml_results) == len(test_std_err_results))

        # Compare results
        for i in range(len(out_ml_results)):
            ml_row = out_ml_results[i]
            std_err_row = out_std_err_results[i]

            for i in range(len(test_ml_results)):
                if np.all(np.isclose(ml_row, test_ml_results[i])) and \
                        np.all(np.isclose(
                            std_err_row, test_std_err_results[i])):
                    found = True
                    test_ml_results.pop(i)
                    test_std_err_results.pop(i)
                    break
            if not found:
                raise Exception(
                    'Could not find {}, {} in results'.format(
                        ml_row, std_err_row))

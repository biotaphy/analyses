"""Test configuration fixtures.
"""
import glob
import os


import pytest

# .............................................................................
# .                                 Constants                                 .
# .............................................................................
ALIGNMENTS_DIR = 'alignments'
ANC_STATE_PACKAGES_DIR = 'ancestral_state_packages'
ANC_DIST_PACKAGES_DIR = 'ancestral_distribution_packages'
PHYLO_BETA_DIV_PACKAGES_DIR = 'phylo_beta_diversity'
TREES_DIR = 'trees'

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_DATA_PATH = os.path.join(THIS_DIR, 'data_dir')


# .............................................................................
class SampleDataFiles(object):
    """This class is used to retrieve sample data for the tests.

    Note:
        * For test files, the format should be something like:
            "(in)valid_{name}.{extension}".
    """
    # .....................................
    def get_alignments(self, fmt, is_valid):
        """Gets an alignment file from the sample data.

        Args:
            fmt (str):  The format of the file you want (csv, json, phylip,
                table).
            is_valid (bool): Return valid data files if true, invalid files if
                not.

        Returns:
            A list of alignments matching the arguments.
        """
        ALIGNMENTS_PATH = os.path.join(SAMPLE_DATA_PATH, ALIGNMENTS_DIR)
        return glob.iglob(
            self._get_glob_string(ALIGNMENTS_PATH, is_valid,
                                  self._get_format_extension(fmt)))

    # .....................................
    def get_ancestral_distribution_packages(self, is_valid):
        """Gets a list of the available packages for testing.

        Args:
            is_valid (bool): Return valid data files if true, invalid files if
                not.

        Returns:
            A list of (tree, alignment, results) tuples.
        """
        packages_path = os.path.join(SAMPLE_DATA_PATH, ANC_DIST_PACKAGES_DIR)
        return self._get_packages(packages_path, is_valid)

    # .....................................
    def get_ancestral_state_packages(self, is_valid):
        """Gets a list of the available packages for testing.

        Args:
            is_valid (bool): Return valid data files if true, invalid files if
                not.

        Returns:
            A list of (tree, alignment, results) tuples.
        """
        packages_path = os.path.join(SAMPLE_DATA_PATH, ANC_STATE_PACKAGES_DIR)
        return self._get_packages(packages_path, is_valid)

    # .....................................
    def get_phylo_beta_diversity_packages(self, is_valid):
        """Gets a list of the available phylo beta diversity packages.

        Args:
            is_valid (bool): Return valid data files if true, invalid files if
                not.

        Returns:
            A list of (pam, tree, beta_jac_csv, beta_jne_csv, beta_jtu_csv,
                phylo_beta_jac_csv, phylo_beta_jne_csv, phylo_beta_jtu_csv,
                ) tuples.
        """
        packages_path = os.path.join(
            SAMPLE_DATA_PATH, PHYLO_BETA_DIV_PACKAGES_DIR)
        if is_valid:
            valid_str = 'valid'
        else:
            valid_str = 'invalid'
        package_dirs = glob.glob(os.path.join(packages_path, valid_str, '*'))
        packages = []
        for pkg_dir in package_dirs:
            pam_fn = None
            tree_fn = None
            beta_jac_fn = None
            beta_jne_fn = None
            beta_jtu_fn = None
            phylo_beta_jac_fn = None
            phylo_beta_jne_fn = None
            phylo_beta_jtu_fn = None
            beta_sim_fn = None
            beta_sne_fn = None
            beta_sor_fn = None
            phylo_beta_sim_fn = None
            phylo_beta_sne_fn = None
            phylo_beta_sor_fn = None

            for fn in glob.glob(os.path.join(pkg_dir, '*')):
                basename = os.path.basename(fn)
                if basename.lower().startswith('beta_jac'):
                    beta_jac_fn = fn
                elif basename.lower().startswith('beta_jne'):
                    beta_jne_fn = fn
                elif basename.lower().startswith('beta_jtu'):
                    beta_jtu_fn = fn
                elif basename.lower().startswith('beta_sim'):
                    beta_sim_fn = fn
                elif basename.lower().startswith('beta_sne'):
                    beta_sne_fn = fn
                elif basename.lower().startswith('beta_sor'):
                    beta_sor_fn = fn
                elif basename.lower().startswith('pam'):
                    pam_fn = fn
                elif basename.lower().startswith('phylo_beta_jac'):
                    phylo_beta_jac_fn = fn
                elif basename.lower().startswith('phylo_beta_jne'):
                    phylo_beta_jne_fn = fn
                elif basename.lower().startswith('phylo_beta_jtu'):
                    phylo_beta_jtu_fn = fn
                elif basename.lower().startswith('phylo_beta_sim'):
                    phylo_beta_sim_fn = fn
                elif basename.lower().startswith('phylo_beta_sne'):
                    phylo_beta_sne_fn = fn
                elif basename.lower().startswith('phylo_beta_sor'):
                    phylo_beta_sor_fn = fn
                elif basename.lower().startswith('tree'):
                    tree_fn = fn

            if all([pam_fn, tree_fn, beta_jac_fn, beta_jne_fn, beta_jtu_fn,
                   beta_sim_fn, beta_sne_fn, beta_sor_fn, phylo_beta_jac_fn,
                   phylo_beta_jne_fn, phylo_beta_jtu_fn, phylo_beta_sim_fn,
                   phylo_beta_sne_fn, phylo_beta_sor_fn]):
                packages.append(
                    (pam_fn, tree_fn, beta_jac_fn, beta_jne_fn, beta_jtu_fn,
                     beta_sim_fn, beta_sne_fn, beta_sor_fn, phylo_beta_jac_fn,
                     phylo_beta_jne_fn, phylo_beta_jtu_fn, phylo_beta_sim_fn,
                     phylo_beta_sne_fn, phylo_beta_sor_fn))

        return packages

    # .....................................
    def get_trees(self, fmt, is_valid):
        """Gets an alignment file from the sample data.

        Args:
            fmt (str): The format of the file you want (newick, nexus).
            is_valid (bool): Return valid data files if true, invalid files if
                not.

        Returns:
            A list of tree filenames.
        """
        TREE_PATH = os.path.join(SAMPLE_DATA_PATH, TREES_DIR)
        return glob.iglob(
            self._get_glob_string(TREE_PATH, is_valid,
                                  self._get_format_extension(fmt)))

    # .....................................
    def _get_format_extension(self, fmt):
        """Get the file extension for a format.

        Args:
            fmt (str): The name of a format to get the file extension for.

        Returns:
            str: A file extension for the specified format.

        Raises:
            Exception: If the specified format was not found.
        """
        if fmt.lower() == 'csv':
            return '.csv'
        elif fmt.lower() == 'json':
            return '.json'
        elif fmt.lower() == 'newick':
            return '.tre'
        elif fmt.lower() == 'nexus':
            return '.nex'
        elif fmt.lower() == 'nexml':
            return '.xml'
        elif fmt.lower() == 'phylip':
            return '.phylip'
        elif fmt.lower() == 'table':
            return '.tbl'
        else:
            raise Exception('Cannot handle format: {}'.format(fmt))

    # .....................................
    def _get_glob_string(self, search_dir, is_valid, fmt_ext):
        """Get a glob string for returning files.

        Args:
            search_dir (str): A directory to search for files.
            is_valid (bool): If True, look for files with 'valid' in filename
                otherwise look for 'invalid'.
            fmt_ext (str): Look for files with this file extension.

        Returns:
            A string that can be sent to glob to retrieve files matching the
                provided parameters.
        """
        if is_valid:
            valid_str = 'valid'
        else:
            valid_str = 'invalid'
        return os.path.join(search_dir, '{}_*{}'.format(valid_str, fmt_ext))

    # .....................................
    def _get_packages(self, packages_path, is_valid):
        """Gets a list of the available packages for testing.

        Args:
            packages_path (str): A directory to look for packages within.
            is_valid (bool): Return valid data files if true, invalid files if
                not.


        Returns:
            A list of (tree, alignment, results) tuples.
        """
        packages = []
        package_dirs = glob.glob(
            self._get_glob_string(packages_path, is_valid, ''))
        for pkg_dir in package_dirs:
            tree_fn = None
            align_fn = None
            results_fn = None
            for fn in glob.glob(os.path.join(pkg_dir, '*')):
                basename = os.path.basename(fn)
                if basename.lower().startswith('tree'):
                    tree_fn = fn
                elif basename.lower().startswith('align'):
                    align_fn = fn
                elif basename.lower().startswith('result'):
                    results_fn = fn
            if tree_fn is not None and align_fn is not None and \
                    results_fn is not None:
                packages.append((tree_fn, align_fn, results_fn))
        return packages


# .............................................................................
@pytest.fixture(scope="session")
def data_files():
    """Gets test fixture used to retrieve sample data files.

    Returns:
        A `SampleDataFiles` object.
    """
    return SampleDataFiles()


# .............................................................................
def pytest_generate_tests(metafunc):
    """Pytest function for generating tests

    Args:
        metafunc (:obj:`pytest.Metafunc`): Pytest metafunc object passed to
            test hook.

    Note:
        * We are catching this function to parameterize tests here in a central
            location rather than for each test instance
    """
    df = SampleDataFiles()
    # Tuples of (fixture name, parameterization lists)
    fixture_tuples = [
        ('invalid_ancestral_distribution_package',
            df.get_ancestral_distribution_packages(False)),
        ('invalid_ancestral_state_package',
            df.get_ancestral_state_packages(False)),
        ('invalid_csv_alignment', df.get_alignments('csv', False)),
        ('invalid_json_alignment', df.get_alignments('json', False)),
        ('invalid_phylip_alignment', df.get_alignments('phylip', False)),
        ('invalid_table_alignment', df.get_alignments('table', False)),
        ('valid_ancestral_distribution_package',
            df.get_ancestral_distribution_packages(True)),
        ('valid_ancestral_state_package',
            df.get_ancestral_state_packages(True)),
        ('valid_csv_alignment', df.get_alignments('csv', True)),
        ('valid_json_alignment', df.get_alignments('json', True)),
        ('valid_newick_tree', df.get_trees('newick', True)),
        ('valid_nexml_tree', df.get_trees('nexml', True)),
        ('valid_nexus_tree', df.get_trees('nexus', True)),
        ('valid_phylip_alignment', df.get_alignments('phylip', True)),
        ('valid_phylo_beta_diversity_package',
         df.get_phylo_beta_diversity_packages(True)),
        ('valid_table_alignment', df.get_alignments('table', True))
    ]
    for fixture_name, fixture_values in fixture_tuples:
        if fixture_name in metafunc.fixturenames:
            metafunc.parametrize(fixture_name, fixture_values)

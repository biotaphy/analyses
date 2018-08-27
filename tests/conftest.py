"""
@summary: Test configuration fixtures
"""
import glob
import os
import pytest

# .............................................................................
# .                                 Constants                                 .
# .............................................................................
ALIGNMENTS_DIR = 'alignments'
PACKAGES_DIR = 'packages'
TREES_DIR = 'trees'

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_DATA_PATH = os.path.join(THIS_DIR, 'data_dir')


# .............................................................................
class SampleDataFiles(object):
    """
    @summary: This class is used to retrieve sample data for the tests
    @note: For test files, the format should be something like:
                "(in)valid_{name}.{extension}"
    """
    # .....................................
    def get_alignments(self, fmt, is_valid):
        """
        @summary: Get an alignment file from the sample data
        @param fmt: The format of the file you want (csv, json, phylip, table)
        @param is_valid: Return valid data files if true, invalid files if not
        @note: Will return a list
        """
        ALIGNMENTS_PATH = os.path.join(SAMPLE_DATA_PATH, ALIGNMENTS_DIR)
        return glob.iglob(
            self._get_glob_string(ALIGNMENTS_PATH, is_valid,
                                  self._get_format_extension(fmt)))

    # .....................................
    def get_packages(self, is_valid):
        """
        @summary: Get a list of the available packages for testing
        @note: All packages are assumed to be correct as of now
        @note: Will return list of (tree, alignment, results) tuples
        """
        PACKAGES_PATH = os.path.join(SAMPLE_DATA_PATH, PACKAGES_DIR)
        packages = []
        package_dirs = glob.glob(self._get_glob_string(PACKAGES_PATH,
                                                       is_valid, ''))
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

    # .....................................
    def get_trees(self, fmt, is_valid):
        """
        @summary: Get an alignment file from the sample data
        @param fmt: The format of the file you want (newick, nexus)
        @param is_valid: Return valid data files if true, invalid files if not
        @note: Will return a list
        """
        TREE_PATH = os.path.join(SAMPLE_DATA_PATH, TREES_DIR)
        return glob.iglob(
            self._get_glob_string(TREE_PATH, is_valid,
                                  self._get_format_extension(fmt)))

    # .....................................
    def _get_format_extension(self, fmt):
        if fmt.lower() == 'csv':
            return '.csv'
        elif fmt.lower() == 'json':
            return '.json'
        elif fmt.lower() == 'newick':
            return '.tre'
        elif fmt.lower() == 'nexus':
            return '.nex'
        elif fmt.lower() == 'phylip':
            return '.phylip'
        elif fmt.lower() == 'table':
            return '.tbl'
        else:
            raise Exception('Cannot handle format: {}'.format(fmt))

    # .....................................
    def _get_glob_string(self, search_dir, is_valid, fmt_ext):
        if is_valid:
            valid_str = 'valid'
        else:
            valid_str = 'invalid'
        return os.path.join(search_dir, '{}_*{}'.format(valid_str, fmt_ext))


# .............................................................................
@pytest.fixture(scope="session")
def data_files():
    return SampleDataFiles()

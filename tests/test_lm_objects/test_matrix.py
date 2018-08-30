"""This module tests the ancestral_reconstruction/lm_objects/matrix.py module

The functions in this module are pytest style tests for the matrix.py module
"""
import io
import numpy as np
import pytest

from ancestral_reconstruction.lm_objects import matrix


# .............................................................................
def get_random_matrix(*dim_size):
    """Generates a randomized matrix with the shape provided in dim_size
    """
    headers = {}
    i = 0
    for i in range(len(dim_size)):
        headers[str(i)] = [
            'header-{}-{}'.format(i, x) for x in range(dim_size[i])]
    return matrix.Matrix(np.random.rand(*dim_size), headers=headers)


# .............................................................................
class Test_Matrix(object):
    """Test the Matrix class
    """
    # .....................................
    def test_load(self):
        """Test the load class method
        """
        orig_mtx = get_random_matrix(5, 5)

        # Create a file like object and save original matrix
        mtx_stringio = io.BytesIO()
        orig_mtx.save(mtx_stringio)
        mtx_stringio.seek(0)

        # Attempt to load matrix
        loaded_mtx = matrix.Matrix.load(mtx_stringio)

        # Verify data and headers are the same
        assert np.isclose(loaded_mtx.data.sum(), orig_mtx.data.sum())
        assert loaded_mtx.get_headers() == orig_mtx.get_headers()

    # .....................................
    def test_load_new(self):
        pass

    # .....................................
    def test_load_from_csv(self):
        pass

    # .....................................
    def test_concatenate(self):
        """Test the concatenate function
        """
        mtx1 = get_random_matrix(3, 2)
        mtx2 = get_random_matrix(3, 3)
        mtx3 = get_random_matrix(6, 5)

        # Concatenate matrices
        mtx4 = matrix.Matrix.concatenate([mtx1, mtx2], axis=1)
        mtx5 = matrix.Matrix.concatenate([mtx4, mtx3], axis=0)

        # Check that shapes are what we expect
        assert mtx4.data.shape == (mtx1.data.shape[0],
                                   mtx1.data.shape[1] + mtx2.data.shape[1])
        assert mtx5.data.shape == (mtx4.data.shape[0] + mtx3.data.shape[0],
                                   mtx4.data.shape[1])

    # .....................................
    def test_append(self):
        """Test append
        """
        x1, y1 = (3, 4)
        x2, y2 = (6, 4)
        mtx1 = get_random_matrix(x1, y1)
        mtx2 = get_random_matrix(x2, y2)
        mtx1.append(mtx2, axis=0)

        # Verify the new shape
        assert mtx1.data.shape == (x1 + x2, y1)

    # .....................................
    def test_flatten_2D(self):
        """Test flatten_2D method
        """
        x, y, z = 5, 5, 3
        mtx = get_random_matrix(x, y, z)
        flat_mtx = mtx.flatten_2D()

        # Test that there are two dimensions of headers and data
        assert len(flat_mtx.data.shape) == 2
        assert len(flat_mtx.get_headers().keys()) == 2

        # Test that the length of the row headers is now 2 (accounting for
        #    the flattened dimension
        for h in flat_mtx.get_row_headers():
            assert len(h) == 2

        # Check that the data shape in the row dimension is z * old shape
        assert flat_mtx.data.shape[0] == x * z

    # .....................................
    def test_get_column_headers(self):
        """Test get_column_headers
        """
        mtx = get_random_matrix(3, 8)
        col_headers = mtx.get_column_headers()
        assert isinstance(col_headers, list)
        assert len(col_headers) == mtx.data.shape[1]

    # .....................................
    def test_get_headers(self):
        """Test get_headers
        """
        mtx = get_random_matrix(2, 2, 4)
        headers = mtx.get_headers()
        assert isinstance(headers, dict)
        for i in range(len(mtx.data.shape)):
            assert len(headers[str(i)]) == mtx.data.shape[i]
        assert mtx.get_headers(axis=9999) is None

    # .....................................
    def test_get_row_headers(self):
        """Test get_row_headers
        """
        mtx = get_random_matrix(3, 8)
        row_headers = mtx.get_row_headers()
        assert isinstance(row_headers, list)
        assert len(row_headers) == mtx.data.shape[0]

    # .....................................
    def test_save(self):
        pass

    # .....................................
    def test_set_column_headers(self):
        pass

    # .....................................
    def test_set_headers(self):
        pass

    # .....................................
    def test_set_row_headers(self):
        pass

    # .....................................
    def test_slice(self):
        pass

    # .....................................
    def test_slice_by_header(self):
        pass

    # .....................................
    def test_write_csv(self):
        pass


# .............................................................................
class Test_ArrayStream(object):
    """Test the ArrayStream class
    """
    # .....................................
    def test_gen_1d(self):
        """Test the gen function with a 1d array

        Note:
            gen() is an iterator and can be tested as such
        """
        arr = np.array([1, 2, 3])
        stream = matrix.ArrayStream(arr)
        assert len(stream) == arr.shape[0]
        for i in stream:
            assert i is not None

    # .....................................
    def test_gen_2d(self):
        """Test the gen function with a 2d array

        Note:
            gen() is an iterator and can be tested as such
        """
        arr = np.array([[1, 2, 3], [4, 5, 6]])
        stream = matrix.ArrayStream(arr)
        assert len(stream) == arr.shape[0]
        for i in stream:
            assert i is not None

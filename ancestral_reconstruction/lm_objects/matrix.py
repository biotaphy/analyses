"""
@summary Module that contains a Matrix class that has header information
@author CJ Grady
@status: alpha
@version: 2.0
@license: gpl2
@copyright: Copyright (C) 2018, University of Kansas Center for Research

          Lifemapper Project, lifemapper [at] ku [dot] edu,
          Biodiversity Institute,
          1345 Jayhawk Boulevard, Lawrence, Kansas, 66045, USA

          This program is free software; you can redistribute it and/or modify
          it under the terms of the GNU General Public License as published by
          the Free Software Foundation; either version 2 of the License, or (at
          your option) any later version.

          This program is distributed in the hope that it will be useful, but
          WITHOUT ANY WARRANTY; without even the implied warranty of
          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
          General Public License for more details.

          You should have received a copy of the GNU General Public License
          along with this program; if not, write to the Free Software
          Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
          02110-1301, USA.
@todo: Use https://docs.scipy.org/doc/numpy/user/basics.subclassing.html when
          changing this to subclass numpy.ndarray
@todo: Handle multiple rows / columns / etc of headers like:
            (PAM x, y, site ids)
@todo: Load should handle compressed and not compressed
"""
from copy import deepcopy
import io
import json
import numpy as np


HEADERS_KEY = 'headers'
DATA_KEY = 'data'
VERSION_KEY = 'version'
VERSION = '2.0.0'


# .............................................................................
class Matrix(object):
    """
    @summary: Lifemapper wrapper for Numpy ndarrays that adds headers
    """
    # ...........................
    def __init__(self, mtx, headers=None):
        """
        @summary: Constructor
        @param mtx: A matrix (like) object to use as the base data for the
                    Matrix.  This can be None if the data has not been
                    initialized
        @param headers: Optional headers for this matrix.  This may be either a
                          list of lists, where the index of a list in the lists
                          will be treated as the axis
                          (ex. [['Row 1', 'Row 2', 'Row 3'],
                                ['Column 1', 'Column 2']])
                          Or this could be a dictionary where the key is used
                          for the axis.  (Ex:
                          {
                             '1' : ['Column 1', 'Column 2'],
                             '0' : ['Row 1', 'Row 2', 'Row 3']
                          }
        @note: If the headers for an axis are a string and not a list, it will
                    be treated as a file name
        """
        self.data = mtx
        self.headers = {}
        if headers is not None:
            self.set_headers(headers)

    # ...........................
    @classmethod
    def load(cls, flo):
        """
        @summary: Attempt to load a Matrix object from a file
        @param flo: A file-like object for the stored Matrix, Numpy array, or ?
        """
        # Try loading the Matrix object
        try:
            return cls.load_new(flo)
        except:
            # Try loading a numpy array
            try:
                # Seek back to start of file
                flo.seek(0)
                data = np.load(flo)
                return cls(data)
            except Exception as e:
                raise IOError('{} : {}'.format(
                    'Cannot load matrix data from file-like object provided',
                    str(e)))

    # ...........................
    @classmethod
    def load_new(cls, flo):
        """
        @summary: Attempt to load a Matrix object from a file
        @param flo: A file-like object with matrix data
        """
        header_lines = []
        data_lines = []
        do_headers = True
        data_stream = io.BytesIO()
        for line in flo:
            try:
                str_line = line.decode('utf-8')
                if do_headers:
                    if str_line.startswith(DATA_KEY):
                        do_headers = False
                    else:
                        header_lines.append(str_line)
            except:
                data_stream.write(line)
                # data_lines.append(line)

        my_obj = json.loads(''.join(header_lines))

        headers = my_obj[HEADERS_KEY]
        # Load returns a tuple if compressed
        data_stream.seek(0)
        tmp = np.load(data_stream)

        # Get data from temp object, we will return data item for first key
        data = tmp[list(tmp.keys())[0]]
        return cls(data, headers=headers)

    # ...........................
    @classmethod
    def load_from_csv(cls, flo):
        """
        @summary: Loads a Matrix from a CSV file
        @param flo: A string (filename) or file-like object containing a CSV
        """
        raise Exception('Not implemented')

    # ...........................
    @classmethod
    def concatenate(cls, mtx_list, axis=0):
        """
        @summary: Concatenates multiple Matrix objects together to form a new
                        Matrix object
        @param mtx_list: A List of Matrix objects to concatenate together
        @param axis: The axis to concatenate these Matrix objects on
        @note: Assumes that headers for other axes are the same
        """
        mtx_objs = []
        axis_headers = []
        for mtx in mtx_list:
            if not isinstance(mtx, Matrix):
                mtx = Matrix(mtx)
            if mtx.data is not None:
                # Make sure we reshape if necessary if adding new axis
                #    (stacking)
                if mtx.data.ndim < axis + 1:  # Add 1 since zero-based
                    new_shape = list(mtx.data.shape) + [1]
                    mtx.data = mtx.data.reshape(new_shape)
                    mtx.set_headers([''], axis=str(axis))

                h = mtx.get_headers(axis=str(axis))
                if h is None:
                    h = ['']
                axis_headers.extend(h)
                mtx_objs.append(mtx.data)

        # Create a new data matrix
        new_data = np.concatenate(mtx_objs, axis=axis)
        # Use the first Matrix's headers as the base
        new_headers = mtx_list[0].get_headers()
        # Replace the axis of headers with the concatenated version
        new_headers[str(axis)] = axis_headers
        return cls(new_data, headers=new_headers)

    # ...........................
    def append(self, mtx, axis=0):
        """
        @summary: Appends the provided Matrix object to this one
        @param mtx: The Matrix object to append to this one
        @param axis: The axis to append this matrix on
        @note: Only keeps the headers for the append axis, assumes the other
                   axes are the same
        """
        self.data = np.append(self.data, mtx.data, axis=axis)
        self.headers[str(axis)].append(mtx.get_headers(axis=axis))

    # ...........................
    def flatten_2D(self):
        """
        @summary: Flattens a higher dimension Matrix object into a 2D matrix
        """
        flat_mtx = self
        while flat_mtx.data.ndim > 2:
            # More than two dimensions so we must flatten
            old_shape = flat_mtx.data.shape
            old_num_rows = old_shape[0]
            new_shape = tuple([old_shape[0]*old_shape[2],
                               old_shape[1]] + list(old_shape[3:]))
            new_mtx = Matrix(np.zeros(new_shape))

            old_rh = flat_mtx.get_row_headers()
            new_rh = []

            # Get old headers
            try:
                old_headers = flat_mtx.get_headers(axis=2)
            except KeyError:
                old_headers = [''] * old_shape[2]

            # Set data and headers
            for i in range(old_shape[2]):
                oh = old_headers[i]
                # Set data
                start_row = i * old_num_rows
                end_row = (i+1) * old_num_rows
                new_mtx.data[start_row:end_row, :] = flat_mtx.data[:, :, i]

                # Set row headers
                for rh in old_rh:
                    if not isinstance(rh, list):
                        rh = [rh]
                    new_rh.append(rh+[oh])

            # Set the headers on the new matrix
            new_mtx.set_row_headers(new_rh)
            new_mtx.set_column_headers(flat_mtx.get_column_headers())

            # Higher order headers
            for axis in flat_mtx.headers.keys():
                if int(axis) > 2:
                    # Reduce the key of the axis by one and set headers on
                    #   new matrix
                    new_mtx.set_headers(flat_mtx.get_headers(axis=axis),
                                        axis=str(int(axis) - 1))

            flat_mtx = new_mtx

        return flat_mtx

    # ...........................
    def get_column_headers(self):
        """
        @summary: Shortcut to get column headers
        @todo: Throw a different exception if no column header?
        """
        return self.get_headers(axis=1)

    # ...........................
    def get_headers(self, axis=None):
        """
        @summary: Get the headers associated with this Matrix, optionally
                     limited to a specific axis
        @param axis: If provided, return headers for this axis, else, return
                        all
        """
        if axis is None:
            return self.headers
        else:
            if str(axis) in self.headers.keys():
                return self.headers[str(axis)]
            else:
                return None

    # ...........................
    def get_row_headers(self):
        """
        @summary: Shortcut to get row headers
        @todo: Throw a different exception if no row headers?
        """
        return self.get_headers(axis=0)

    # ...........................
    def save(self, flo):
        """
        @summary: Saves the Matrix object in a JSON / Numpy hybrid format to
                    the file-like object
        @param flo: The file-like object to write to
        """
        my_obj = {}
        my_obj[HEADERS_KEY] = self.headers
        my_obj[VERSION_KEY] = VERSION
        try:
            flo.write('{}\n'.format(json.dumps(my_obj, indent=3,
                                               default=float)))
            flo.write('{}\n'.format(DATA_KEY))
            np.savez_compressed(flo, self.data)
        except:
            my_obj_str = '{}\n{}\n'.format(
                json.dumps(my_obj, indent=3, default=float), DATA_KEY)
            flo.write(bytes(my_obj_str, 'utf8'))
            np.savez_compressed(flo, self.data)

    # ...........................
    def set_column_headers(self, headers):
        """
        @summary: Shortcut to set column headers
        """
        self.set_headers(headers, axis=1)

    # ...........................
    def set_headers(self, headers, axis=None):
        """
        @summary: Set the headers for this Matrix, optionally for a specific
                    axis
        @param headers: Matrix headers.  Can be a list of lists, a dictionary
                          of lists, or if axis is provided, a single list
        @param axis: If provided, set the headers for a specific axis, else,
                      process as if it is for the entire Matrix
        @todo: Validate input for single axis operation?
        @note: Resets headers dictionary when setting values for all headers
        @note: Duck types to use list of lists or dictionary to set values for
                different axes
        """
        if axis is not None:
            self.headers[str(axis)] = headers
        else:
            self.headers = {}
            try:
                headers_keys = headers.keys()
            except:  # Not a dictionary
                # Check if first item is a list
                if isinstance(headers[str(0)], list):
                    # Assume list of lists
                    headers_keys = range(len(headers))
                else:
                    # Convert to a list
                    headers = [headers]
                    headers_keys = [0]

            # We should have a list of keys, which could be either dictionary
            #    keys or list indices
            for k in headers_keys:
                self.headers[str(k)] = headers[str(k)]

    # ...........................
    def set_row_headers(self, headers):
        """
        @summary: Shortcut to set row headers
        """
        self.set_headers(headers, axis=0)

    # ...........................
    def slice(self, *args):
        """
        @summary: Subsets the matrix and returns a new instance
        @param *args: These are iterables for the indices to retrieve
        @note: The first parameter will be for axis 0, second for axis 1, etc
        """
        new_data = np.copy(self.data)
        new_headers = deepcopy(self.headers)
        # For each arg in the list
        for i in range(len(args)):
            # Subset the data matrix
            new_data = new_data.take(args[i], axis=i)
            # Subset the headers
            tmp = []
            for j in args[i]:
                tmp.append(new_headers[str(i)][j])
            new_headers[str(i)] = tmp
        return Matrix(new_data, headers=new_headers)

    # ...........................
    def slice_by_header(self, header, axis):
        """
        @summary: Gets a slice of the Matrix matching the header provided
        @param header: The name of a header to use for slicing
        @param axis: The axis to find this header
        @raise ValueError: If the header is not found for the specified axis
        @todo: Add capability to slice over multiple axes and multiple headers
                Maybe combine with other slice method and provide method to
                search for header indices
        """
        idx = self.headers[str(axis)].index(header)

        new_data = np.copy(np.take(self.data, idx, axis=axis))

        # Need to reshape the result.  Take the existing shape and change the
        #    query axis to 1
        new_shape = list(self.data.shape)
        new_shape[axis] = 1
        new_data = new_data.reshape(new_shape)

        # Copy the headers and set the header for the axis to just be the
        #    search header
        new_headers = deepcopy(self.headers)
        new_headers[str(axis)] = [header]

        # Return a new Matrix
        return Matrix(new_data, headers=new_headers)

    # ...........................
    def write_csv(self, flo, *slice_args):
        """
        @summary: Write the Matrix object to a CSV file-like object
        @param flo: The file-like object to write to
        @param slice_args: If provided, perform a slice operation and use the
                         resulting matrix for writing
        @todo: Handle header overlap (where the header for one axis is for
                  another axis header
        @note: Currently only works for 2-D tables
        """
        if list(slice_args):
            mtx = self.slice(*slice_args)
        else:
            mtx = self

        if mtx.data.ndim > 2:
            mtx = mtx.flatten_2D()

        # .....................
        # Inner functions

        # .....................
        def already_lists(x):
            """
            @summary: Use this function for processing headers when they are
                        already lists
            """
            return x

        # .....................
        def make_lists(x):
            """
            @summary: Use this function for processing headers when they are
                        not lists to start
            """
            return [x]

        # .....................
        def csv_generator():
            """
            @summary: This function is a generator that yields rows of values
                       to be output as CSV
            """
            try:
                row_headers = mtx.headers['0']
            except:
                # No row headers
                row_headers = [[] for _ in range(mtx.data.shape[0])]

            if isinstance(row_headers[0], list):
                listify = already_lists
            else:
                listify = make_lists

            # Start with the header row, if we have one
            if '1' in mtx.headers and mtx.headers['1']:
                # Add a blank entry if we have row headers
                header_row = ['']*len(listify(
                                    row_headers[0]) if row_headers else [])
                header_row.extend(mtx.headers['1'])
                yield header_row
            # For each row in the data set
            for i in range(mtx.data.shape[0]):
                # Add the row headers if exists
                row = []
                row.extend(listify(row_headers[i]))
                # Get the data from the data array
                row.extend(mtx.data[i].tolist())
                yield row

        # .....................
        # Main write_csv function
        for row in csv_generator():
            flo.write(u"{0}\n".format(','.join([str(v) for v in row])))


# .............................................................................
class ArrayStream(list):
    """
    @summary: Generator class for a numpy array for JSON serialization
    @note: This is done to save memory rather than creating a list of the
             entire array / matrix.  It is used by the JSON encoder to write
             the data to file
    """
    # ...........................
    def __init__(self, x):
        """
        @summary: Constructor
        @param x: The numpy array to stream
        """
        self.x = x
        self.my_len = self.x.shape[0]

    # ...........................
    def __iter__(self):
        """
        @summary: Iterator for array
        """
        return self.gen()

    # ...........................
    def __len__(self):
        """
        @summary: Length function
        """
        return self.my_len

    # ...........................
    def gen(self):
        """
        @summary: Generator function.  Loop over array and create ArrayStreams
                   for sub arrays
        """
        n = 0

        while n < self.my_len:
            if isinstance(self.x[n], np.ndarray):
                yield ArrayStream(self.x[n])
            else:
                yield self.x[n]
            n += 1

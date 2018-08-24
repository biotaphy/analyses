"""
@summary: This module contains classes and functions for testing the provided
             data readers
@note: Uses pytest style testing
"""
import pytest
from StringIO import StringIO

import ancestral_reconstruction.helpers.data_readers as dr
from ancestral_reconstruction.helpers.sequence import Sequence


# .............................................................................
class Test_create_sequence_list_from_dict(object):
    """
    @summary: Tests the create_sequence_list_from_dict method
    """
    # .....................................
    def test_invalid_dict(self):
        pass

    # .....................................
    def test_valid_dict(self):
        pass


# .............................................................................
class Test_read_csv_alignment_flo(object):
    """
    @summary: Tests the dr.read_csv_alignment_flo method
    """
    # .....................................
    def test_file_invalid(self, data_files):
        """
        @summary: Tests that the invalid alignment files fail properly
        @param data_files: A pytest fixture defined in conftest.py
        """
        invalid_csv_alignments = data_files.get_alignments('csv', False)
        for fn in invalid_csv_alignments:
            with open(fn) as in_csv:
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_csv_alignment_flo(in_csv)

    # .....................................
    def test_file_valid(self, data_files):
        """
        @summary: Tests that the valid alignment files do not fail
        @param data_files: A pytest fixture defined in conftest.py
        """
        valid_csv_alignments = data_files.get_alignments('csv', True)
        for fn in valid_csv_alignments:
            with open(fn) as in_csv:
                try:
                    sequence_list, headers = dr.read_csv_alignment_flo(in_csv)
                    assert len(headers) > 0
                    assert len(sequence_list) > 0
                    for i in sequence_list:
                        assert isinstance(i, Sequence)
                        assert isinstance(i.name, basestring)
                        assert len(i.cont_values) > 0
                except Exception as e:
                    print 'Raised exception:', str(e)
                    assert False

    # .....................................
    def test_stringio_invalid(self, data_files):
        """
        @summary: Tests that the invalid alignment files fail properly when
                   loaded into StringIO objects
        @param data_files: A pytest fixture defined in conftest.py
        """
        invalid_csv_alignments = data_files.get_alignments('csv', False)
        for fn in invalid_csv_alignments:
            with open(fn) as in_csv:
                csv_stringio = StringIO()
                csv_stringio.write(in_csv.read())
                csv_stringio.seek(0)
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_csv_alignment_flo(csv_stringio)

    # .....................................
    def test_stringio_valid(self, data_files):
        """
        @summary: Tests that the valid alignment files do not fail when loaded
                     into StringIO objects
        @param data_files: A pytest fixture defined in conftest.py
        """
        valid_csv_alignments = data_files.get_alignments('csv', True)
        for fn in valid_csv_alignments:
            with open(fn) as in_csv:
                csv_stringio = StringIO()
                csv_stringio.write(in_csv.read())
                csv_stringio.seek(0)
                try:
                    sequence_list, headers = dr.read_csv_alignment_flo(
                                                                csv_stringio)
                    assert len(headers) > 0
                    assert len(sequence_list) > 0
                    for i in sequence_list:
                        assert isinstance(i, Sequence)
                        assert isinstance(i.name, basestring)
                        assert len(i.cont_values) > 0
                except Exception as e:
                    print 'Raised exception:', str(e)
                    assert False


# .............................................................................
class Test_read_json_alignment_flo(object):
    """
    @summary: Tests the dr.read_json_alignment_flo method
    """
    # .....................................
    def test_file_invalid(self, data_files):
        """
        @summary: Tests that the invalid alignment files fail properly
        @param data_files: A pytest fixture defined in conftest.py
        """
        invalid_json_alignments = data_files.get_alignments('json', False)
        for fn in invalid_json_alignments:
            with open(fn) as in_json:
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_json_alignment_flo(in_json)

    # .....................................
    def test_file_valid(self, data_files):
        pass

    # .....................................
    def test_stringio_invalid(self, data_files):
        pass

    # .....................................
    def test_stringio_valid(self, data_files):
        pass


# .............................................................................
class Test_read_phylip_alignment_flo(object):
    """
    @summary: Tests the dr.read_phylip_alignment_flo method
    """
    # .....................................
    def test_file_invalid(self, data_files):
        """
        @summary: Tests that the invalid alignment files fail properly
        @param data_files: A pytest fixture defined in conftest.py
        """
        invalid_phylip_alignments = data_files.get_alignments('phylip', False)
        for fn in invalid_phylip_alignments:
            with open(fn) as in_phylip:
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_phylip_alignment_flo(in_phylip)

    # .....................................
    def test_file_valid(self, data_files):
        pass

    # .....................................
    def test_stringio_invalid(self, data_files):
        pass

    # .....................................
    def test_stringio_valid(self, data_files):
        pass


# .............................................................................
class Test_read_phylip_continuous_alignment_flo(object):
    """
    @summary: Tests the read_phylip_continuous_alignment_flo method
    """
    # .....................................
    def test_file_invalid(self, data_files):
        pass

    # .....................................
    def test_file_valid(self, data_files):
        pass

    # .....................................
    def test_stringio_invalid(self, data_files):
        pass

    # .....................................
    def test_stringio_valid(self, data_files):
        pass


# .............................................................................
class Test_read_table_continuous_alignment_flo(object):
    """
    @summary: Tests the read_table_continuous_alignment_flo method
    """
    # .....................................
    def test_file_invalid(self, data_files):
        """
        @summary: Tests that the invalid alignment files fail properly
        @param data_files: A pytest fixture defined in conftest.py
        """
        invalid_table_alignments = data_files.get_alignments('table', False)
        for fn in invalid_table_alignments:
            with open(fn) as in_table:
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_table_alignment_flo(in_table)

    # .....................................
    def test_file_valid(self, data_files):
        pass

    # .....................................
    def test_stringio_invalid(self, data_files):
        pass

    # .....................................
    def test_stringio_valid(self, data_files):
        pass

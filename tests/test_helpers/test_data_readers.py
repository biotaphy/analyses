"""This module contains classes and functions for testing data readers.

Note:
    * Uses pytest style testing.
"""
try:
    from StringIO import StringIO
    string_formats = (basestring)
except:
    # Python 3
    from io import StringIO
    string_formats = (str)

import pytest

from lmpy import Matrix

import analyses.helpers.data_readers as dr
from analyses.helpers.sequence import Sequence


# .............................................................................
class Test_create_sequence_list_from_dict(object):
    """Test class for the create_sequence_list_from_dict method.
    """
    # .....................................
    def test_invalid_dict(self):
        """Tests that function fails with a dict that has the wrong structure.
        """
        test_dict = {
            "name1": [1, 2, 3],
            "name2": "bad value"
        }
        with pytest.raises(dr.AlignmentIOError):
            dr.create_sequence_list_from_dict(test_dict)

    # .....................................
    def test_valid_dict(self):
        """Tests that the function operates correctly with valid input data.
        """
        test_dict = {
            "A": [0.9, 0.2, 0.2, 0.3, 0.4, 0.4],
            "B": [0.01, 0.1, 0.2, 0.3, 0.4, 0.4],
            "C": [0.8, 0.1, 0.2, 0.3, 0.4, 0.4],
            "D": [0.3, 0.1, 0.2, 0.3, 0.4, 0.4],
            "E": [0.001, 0.1, 0.2, 0.3, 0.4, 0.4],
            "F": [0.11, 0.1, 0.2, 0.3, 0.4, 0.4],
            "G": [0.99, 0.2, 0.2, 0.3, 0.4, 0.4]
        }
        sequence_list, headers = dr.create_sequence_list_from_dict(test_dict)
        for i in sequence_list:
            assert isinstance(i, Sequence)
            assert isinstance(i.name, string_formats)
            assert len(i.cont_values) > 0


# .............................................................................
class Test_get_character_matrix_from_sequences_list(object):
    """Test class for get_character_matrix_from_sequences_list method.
    """
    # .....................................
    def test_valid_no_headers(self):
        """Test the function with valid Sequences and no variable headers.
        """
        seq1 = Sequence(name='seq1')
        seq1.set_cont_values([1.0, 2.0, 3.0])
        seq2 = Sequence(name='seq2')
        seq2.set_cont_values([3.0, 3.2, 4.1])
        test_sequence_list = [seq1, seq2]
        mtx = dr.get_character_matrix_from_sequences_list(test_sequence_list)
        assert isinstance(mtx, Matrix)

    # .....................................
    def test_valid_with_headers(self):
        """Test the function with valid Sequences and include headers.
        """
        seq1 = Sequence(name='seq1')
        seq1.set_cont_values([1.0, 2.0, 3.0])
        seq2 = Sequence(name='seq2')
        seq2.set_cont_values([3.0, 3.2, 4.1])
        test_sequence_list = [seq1, seq2]
        col_headers = ['Val 1', 'Val 2', 'Val 3']
        mtx = dr.get_character_matrix_from_sequences_list(
            test_sequence_list, var_headers=col_headers)
        assert isinstance(mtx, Matrix)
        assert mtx.get_column_headers() == col_headers


# .............................................................................
class Test_load_alignment_from_filename(object):
    """Tests the dr.load_alignment_from_filename method.
    """
    # .....................................
    def test_csv_file(self, valid_csv_alignment):
        """Tests that a valid CSV file can be loaded into an alignment.

        Args:
            valid_csv_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid csv alignment filenames.
        """
        sequences, headers = dr.load_alignment_from_filename(
            valid_csv_alignment)
        assert len(headers) > 0
        assert len(sequences) > 0
        for i in sequences:
            assert isinstance(i, Sequence)
            assert isinstance(i.name, string_formats)
            assert len(i.cont_values) > 0

    # .....................................
    def test_invalid_file(self):
        """Tests that an invalid file fails properly on load.
        """
        with pytest.raises(RuntimeError):
            dr.load_alignment_from_filename('./file_with_bad_ext.bad')

    # .....................................
    def test_json_file(self, valid_json_alignment):
        """Tests that a valid JSON file can be loaded into an alignment.

        Args:
            valid_json_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid json alignment filenames.
        """
        sequences, headers = dr.load_alignment_from_filename(
            valid_json_alignment)
        if headers is not None:
            assert len(headers) > 0
        assert len(sequences) > 0
        for i in sequences:
            assert isinstance(i, Sequence)
            assert isinstance(i.name, string_formats)
            assert len(i.cont_values) > 0

    # .....................................
    def test_phylip_file(self, valid_phylip_alignment):
        """Tests that a valid phylip file can be loaded into an alignment.

        Args:
            valid_phylip_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid phylip alignment filenames.
        """
        sequences, headers = dr.load_alignment_from_filename(
            valid_phylip_alignment)
        assert headers is None
        assert len(sequences) > 0
        for i in sequences:
            assert isinstance(i, Sequence)
            assert isinstance(i.name, string_formats)
            assert i.seq is not None

    # .....................................
    def test_table_file(self, valid_table_alignment):
        """Tests that a valid table file can be loaded into an alignment.

        Args:
            valid_table_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid table alignment filenames.
        """
        sequences, headers = dr.load_alignment_from_filename(
            valid_table_alignment)
        assert headers is None
        assert len(sequences) > 0
        for i in sequences:
            assert isinstance(i, Sequence)
            assert isinstance(i.name, string_formats)
            assert len(i.cont_values) > 0


# .............................................................................
class Test_read_csv_alignment_flo(object):
    """Tests the dr.read_csv_alignment_flo method.
    """
    # .....................................
    def test_file_invalid(self, invalid_csv_alignment):
        """Tests that the invalid alignment files fail properly.

        Args:
            invalid_csv_alignment (pytest.fixture): A parameterized pytest
                fixture providing invalid csv alignment filenames.
        """
        with open(invalid_csv_alignment) as in_csv:
            with pytest.raises(dr.AlignmentIOError):
                dr.read_csv_alignment_flo(in_csv)

    # .....................................
    def test_file_valid(self, valid_csv_alignment):
        """Tests that the valid alignment files do not fail.

        Args:
            valid_csv_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid csv alignment filenames.
        """
        with open(valid_csv_alignment) as in_csv:
            try:
                sequence_list, headers = dr.read_csv_alignment_flo(in_csv)
                assert len(headers) > 0
                assert len(sequence_list) > 0
                for i in sequence_list:
                    assert isinstance(i, Sequence)
                    assert isinstance(i.name, string_formats)
                    assert len(i.cont_values) > 0
            except Exception as e:
                print('Raised exception: {}'.format(str(e)))
                assert False

    # .....................................
    def test_stringio_invalid(self, invalid_csv_alignment):
        """Tests that invalid alignment files cause the proper failure.

        Tests that the invalid alignment files fail properly when loaded into
        StringIO objects.

        Args:
            invalid_csv_alignment (pytest.fixture): A parameterized pytest
                fixture providing invalid csv alignment filenames.
        """
        with open(invalid_csv_alignment) as in_csv:
            csv_stringio = StringIO()
            csv_stringio.write(in_csv.read())
            csv_stringio.seek(0)
            with pytest.raises(dr.AlignmentIOError):
                dr.read_csv_alignment_flo(csv_stringio)

    # .....................................
    def test_stringio_valid(self, valid_csv_alignment):
        """Tests that alignments are read properly with valid input data.

        Tests that the valid alignment files do not fail when loaded into
        StringIO objects.

        Args:
            valid_csv_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid csv alignment filenames.
        """
        with open(valid_csv_alignment) as in_csv:
            csv_stringio = StringIO()
            csv_stringio.write(in_csv.read())
            csv_stringio.seek(0)
            try:
                (sequence_list, headers) = dr.read_csv_alignment_flo(
                    csv_stringio)
                assert len(headers) > 0
                assert len(sequence_list) > 0
                for i in sequence_list:
                    assert isinstance(i, Sequence)
                    assert isinstance(i.name, string_formats)
                    assert len(i.cont_values) > 0
            except Exception as e:
                print('Raised exception: {}'.format(str(e)))
                assert False


# .............................................................................
class Test_read_json_alignment_flo(object):
    """Test class for the dr.read_json_alignment_flo method.
    """
    # .....................................
    def test_file_invalid(self, invalid_json_alignment):
        """Tests that the invalid alignment files fail properly.

        Args:
            invalid_json_alignment (pytest.fixture): A parameterized pytest
                fixture providing invalid json alignment filenames.
        """
        with open(invalid_json_alignment) as in_json:
            with pytest.raises(dr.AlignmentIOError):
                dr.read_json_alignment_flo(in_json)

    # .....................................
    def test_file_valid(self, valid_json_alignment):
        """Tests that the valid alignment files do not fail.

        Args:
            valid_json_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid json alignment filenames.
        """
        with open(valid_json_alignment) as in_json:
            try:
                (sequence_list, headers
                 ) = dr.read_json_alignment_flo(in_json)
                if headers is not None:
                    assert len(headers) > 0
                assert len(sequence_list) > 0
                for i in sequence_list:
                    assert isinstance(i, Sequence)
                    assert isinstance(i.name, string_formats)
                    assert len(i.cont_values) > 0
            except Exception as e:
                print('Raised exception: {}'.format(str(e)))
                assert False

    # .....................................
    def test_stringio_invalid(self, invalid_json_alignment):
        """Tests that invalid alignment files fail properly.

        Tests that the invalid alignment files fail properly when loaded into
        StringIO objects.

        Args:
            invalid_json_alignment (pytest.fixture): A parameterized pytest
                fixture providing invalid json alignment filenames.
        """
        with open(invalid_json_alignment) as in_json:
            json_stringio = StringIO()
            json_stringio.write(in_json.read())
            json_stringio.seek(0)
            with pytest.raises(dr.AlignmentIOError):
                dr.read_json_alignment_flo(json_stringio)

    # .....................................
    def test_stringio_valid(self, valid_json_alignment):
        """Tests that valid JSON alignment files do not fail.

        Tests that the valid alignment files do not fail when loaded into
        StringIO objects.

        Args:
            valid_json_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid json alignment filenames.
        """
        with open(valid_json_alignment) as in_json:
            json_stringio = StringIO()
            json_stringio.write(in_json.read())
            json_stringio.seek(0)
            try:
                (sequence_list, headers
                 ) = dr.read_json_alignment_flo(json_stringio)
                if headers is not None:
                    assert len(headers) > 0
                assert len(sequence_list) > 0
                for i in sequence_list:
                    assert isinstance(i, Sequence)
                    assert isinstance(i.name, string_formats)
                    assert len(i.cont_values) > 0
            except Exception as e:
                print('Raised exception: {}'.format(str(e)))
                assert False


# .............................................................................
class Test_read_phylip_alignment_flo(object):
    """Test class for the dr.read_phylip_alignment_flo method.
    """
    # .....................................
    def test_file_invalid(self, invalid_phylip_alignment):
        """Tests that the invalid alignment files fail properly.

        Args:
            invalid_phylip_alignment (pytest.fixture): A parameterized pytest
                fixture providing invalid phylip alignment filenames.
        """
        with open(invalid_phylip_alignment) as in_phylip:
            with pytest.raises(dr.AlignmentIOError):
                dr.read_phylip_alignment_flo(in_phylip)

    # .....................................
    def test_file_valid(self, valid_phylip_alignment):
        """Tests that the valid alignment files do not fail.

        Args:
            valid_phylip_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid phylip alignment filenames.
        """
        with open(valid_phylip_alignment) as in_phylip:
            try:
                sequence_list = dr.read_phylip_alignment_flo(in_phylip)
                assert len(sequence_list) > 0
                for i in sequence_list:
                    assert isinstance(i, Sequence)
                    assert isinstance(i.name, string_formats)
                    assert i.seq is not None
            except Exception as e:
                print('Raised exception: {}'.format(str(e)))
                assert False

    # .....................................
    def test_stringio_invalid(self, invalid_phylip_alignment):
        """Test that invalid phylip alignment files cause a failure.

        Tests that the invalid alignment files fail properly when loaded into
        StringIO objects.

        Args:
            invalid_phylip_alignment (pytest.fixture): A parameterized pytest
                fixture providing invalid phylip alignment filenames.
        """
        with open(invalid_phylip_alignment) as in_phylip:
            phylip_stringio = StringIO()
            phylip_stringio.write(in_phylip.read())
            phylip_stringio.seek(0)
            with pytest.raises(dr.AlignmentIOError):
                dr.read_phylip_alignment_flo(phylip_stringio)

    # .....................................
    def test_stringio_valid(self, valid_phylip_alignment):
        """Tests that valid phylip alignment files are loaded properly.

        Tests that the valid alignment files do not fail when loaded into
        StringIO objects.

        Args:
            valid_phylip_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid phylip alignment filenames.
        """
        with open(valid_phylip_alignment) as in_phylip:
            phylip_stringio = StringIO()
            phylip_stringio.write(in_phylip.read())
            phylip_stringio.seek(0)
            try:
                sequence_list = dr.read_phylip_alignment_flo(
                    phylip_stringio)
                assert len(sequence_list) > 0
                for i in sequence_list:
                    assert isinstance(i, Sequence)
                    assert isinstance(i.name, string_formats)
                    assert i.seq is not None
            except Exception as e:
                print('Raised exception: {}'.format(str(e)))
                assert False


# .............................................................................
class Test_read_table_continuous_alignment_flo(object):
    """Test class for the read_table_continuous_alignment_flo method.
    """
    # .....................................
    def test_file_invalid(self, invalid_table_alignment):
        """Tests that the invalid alignment files fail properly.

        Args:
            invalid_table_alignment (pytest.fixture): A parameterized pytest
                fixture providing invalid table alignment filenames.
        """
        with open(invalid_table_alignment) as in_table:
            with pytest.raises(dr.AlignmentIOError):
                dr.read_table_alignment_flo(in_table)

    # .....................................
    def test_file_valid(self, valid_table_alignment):
        """Tests that the valid alignment files do not fail.

        Args:
            valid_table_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid table alignment filenames.
        """
        with open(valid_table_alignment) as in_table:
            try:
                sequence_list = dr.read_table_alignment_flo(in_table)
                assert len(sequence_list) > 0
                for i in sequence_list:
                    assert isinstance(i, Sequence)
                    assert isinstance(i.name, string_formats)
                    assert len(i.cont_values) > 0
            except Exception as e:
                print('Raised exception: {}'.format(str(e)))
                assert False

    # .....................................
    def test_stringio_invalid(self, invalid_table_alignment):
        """Tests that invalid table alignment files cause a failure.

        Tests that the invalid alignment files fail properly when loaded into
        StringIO objects.

        Args:
            invalid_table_alignment (pytest.fixture): A parameterized pytest
                fixture providing invalid table alignment filenames.
        """
        with open(invalid_table_alignment) as in_table:
            table_stringio = StringIO()
            table_stringio.write(in_table.read())
            table_stringio.seek(0)
            with pytest.raises(dr.AlignmentIOError):
                dr.read_table_alignment_flo(table_stringio)

    # .....................................
    def test_stringio_valid(self, valid_table_alignment):
        """Tests that valid table alignment files are loaded correctly.

        Tests that the valid alignment files do not fail when loaded into
        StringIO objects.

        Args:
            valid_table_alignment (pytest.fixture): A parameterized pytest
                fixture providing valid table alignment filenames.
        """
        with open(valid_table_alignment) as in_table:
            table_stringio = StringIO()
            table_stringio.write(in_table.read())
            table_stringio.seek(0)
            try:
                sequence_list = dr.read_table_alignment_flo(table_stringio)
                assert len(sequence_list) > 0
                for i in sequence_list:
                    assert isinstance(i, Sequence)
                    assert isinstance(i.name, string_formats)
                    assert len(i.cont_values) > 0
            except Exception as e:
                print('Raised exception: {}'.format(str(e)))
                assert False

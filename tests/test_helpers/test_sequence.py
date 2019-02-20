"""This module tests the sequence module.

Note:
    * This module uses pytest style tests for the sequence module.
"""
import pytest

from analyses.helpers import sequence


# .............................................................................
class Test_Sequence(object):
    """Test the Sequence class.
    """
    # .....................................
    def test_set_cont_values(self):
        """Test the set_cont_values method.
        """
        seq = sequence.Sequence(name='test')
        new_vals = [1, 2, 3, 4, 5]
        seq.set_cont_values(new_vals)
        assert seq.cont_values == new_vals

    # .....................................
    def test_set_qualstr(self):
        """Test the set_qualstr method.

        Note:
            * Only modifies qualarr if it is an empty list.
        """
        qual = 'abcdefg'
        seq1 = sequence.Sequence()
        seq2 = sequence.Sequence()
        seq2_qualarr = [1, 2, 3, 4, 5]
        seq2.set_qualarr(seq2_qualarr)

        # Set qualstr for sequences
        seq1.set_qualstr(qual)
        seq2.set_qualstr(qual)

        # qualstr should be the same for both
        assert seq1.qualstr == seq2.qualstr

        # qualarr should be different as seq2 should not have changed
        assert seq1.qualarr != seq2.qualarr
        assert seq2.qualarr == seq2_qualarr

    # .....................................
    def test_set_qualarr(self):
        """Test the set_qualarr function.

        Note:
            * Only modifies qualstr if it was an empty string.
        """
        qual = [1, 2, 3, 4, 5]
        seq2_qualstr = 'abcdefg'
        seq1 = sequence.Sequence()
        seq2 = sequence.Sequence()
        seq2.set_qualstr(seq2_qualstr)

        # Set qualarr for sequences
        seq1.set_qualarr(qual)
        seq2.set_qualarr(qual)

        # qualarr should be the same for both
        assert seq1.qualarr == seq2.qualarr

        # qualstr should be different as seq2 should not have changed
        assert seq1.qualstr != seq2.qualstr
        assert seq2.qualstr == seq2_qualstr

    # .....................................
    def test_get_fasta(self):
        """Test the get_fasta method.
        """
        name = 'test_name'
        seq_str = 'abcdefg'
        seq = sequence.Sequence(name=name, seq=seq_str)
        assert seq.get_fasta() == '>{}\n{}'.format(name, seq_str)

    # .....................................
    def test_get_fastq(self):
        """Test the get_fastq method.
        """
        name = 'test_name'
        seq_str = 'abcdefg'
        qual_str = 'abcdefg'
        seq = sequence.Sequence(name=name, seq=seq_str)
        seq.set_qualstr(qual_str)
        assert seq.get_fastq() == '@{}\n{}\n+\n{}'.format(name, seq_str,
                                                          qual_str)

    # .....................................
    def test_constructor(self):
        """Test Sequence construction.
        """
        _ = sequence.Sequence(name='test', seq='aaa')
        _ = sequence.Sequence(name='test')
        _ = sequence.Sequence(seq='aaa')
        _ = sequence.Sequence()

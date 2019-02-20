"""Module containing sequence class.
"""
import os
import sys


# .............................................................................
class Sequence(object):
    """Barebones class for sequences.

    This is a barebones class for sequences.  These can be aligned or not and
    can be any type of alphabet.

    Attributes:
        name (str): A name for this sequence.
        seq (str): A sequence string.
        qualstr (str): A string of sequence characters.
        qualarr (list): A list of offset code points.
        cont_values (list): A list of continuous values for a sequence.
    """
    # .....................................
    def __init__(self, name='', seq=''):
        """Constructor for Sequence.

        Args:
            name (str): A name for this sequence.
            seq (str): A sequence string.
        """
        self.name = name
        self.seq = seq
        self.qualstr = ""
        self.qualarr = []
        self.cont_values = []

    # .....................................
    def set_cont_values(self, values):
        """Set the continuous values for the sequence.

        Args:
            values (list): A list of values.
        """
        self.cont_values = values

    # .....................................
    def set_qualstr(self, qual):
        """Set the qualstr attribute.

        Args:
            qual (str): A new string to use for qualstr.

        Note:
            * An offset of 33 is assumed.

        Todo:
            * Should this reset both or neither?
        """
        self.qualstr = qual
        if len(self.qualarr) == 0:
            for j in self.qualstr:
                self.qualarr.append(ord(j)-33)

    # .....................................
    def set_qualarr(self, qual):
        """Set the qualarr attribute.

        Args:
            qual (:obj:`list` of :obj:`int`): A list of code point integers.

        Note:
            * An offset of 33 is assumed.

        Todo:
            * Should this reset both or neither?
        """
        self.qualarr = qual
        if len(self.qualstr) == 0:
            for j in self.qualarr:
                self.qualstr += chr(j+33)

    # .....................................
    def get_fasta(self):
        """Get a fasta string.
        """
        retstr = ">"
        retstr += self.name
        retstr += "\n"
        retstr += self.seq
        return retstr

    # .....................................
    def get_fastq(self):
        """Get a fastq string.
        """
        retstr = "@"
        retstr += self.name
        retstr += "\n"
        retstr += self.seq
        retstr += "\n+\n"
        retstr += self.qualstr
        return retstr

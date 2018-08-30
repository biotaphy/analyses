"""
@summary: Module containing sequence class
"""
import os
import sys


# .............................................................................
class Sequence(object):
    """
    @summary: This is a barebones class for sequences.  These can be aligned or
                not and can be any type of alphabet.
    """
    # .....................................
    def __init__(self, name="", seq=""):
        self.name = name
        self.seq = seq
        self.qualstr = ""
        self.qualarr = []
        self.cont_values = []

    # .....................................
    def set_cont_values(self, values):
        self.cont_values = values

    # .....................................
    def set_qualstr(self, qual):
        """
        @note: An offset of 33 is assumed
        @todo: Should this reset both or neither?
        """
        self.qualstr = qual
        if len(self.qualarr) == 0:
            for j in self.qualstr:
                self.qualarr.append(ord(j)-33)

    # .....................................
    def set_qualarr(self, qual):
        """
        @note: An offset of 33 is assumed
        @todo: Should this reset both or neither?
        """
        self.qualarr = qual
        if len(self.qualstr) == 0:
            for j in self.qualarr:
                self.qualstr += chr(j+33)

    # .....................................
    def get_fasta(self):
        retstr = ">"
        retstr += self.name
        retstr += "\n"
        retstr += self.seq
        return retstr

    # .....................................
    def get_fastq(self):
        retstr = "@"
        retstr += self.name
        retstr += "\n"
        retstr += self.seq
        retstr += "\n+\n"
        retstr += self.qualstr
        return retstr

"""Module containing tools for reading alignment files in various formats.

Todo:
    * FASTA reader?
"""
import csv
import json
import os
import sys

import numpy as np

from lmpy import Matrix

from analyses.helpers.sequence import Sequence


# .............................................................................
class AlignmentIOError(Exception):
    """Wrapper class for alignment errors.
    """
    pass


# .............................................................................
def create_sequence_list_from_dict(values_dict):
    """Creates a list of sequences from a dictionary

    Args:
        values_dict (dict) : A dictionary of taxon name keys and a list of
            values for each value.

    Note:
        * The dictionary should have structure::

            {
                "{taxon_name}" : [{values}]
            }


    Returns:
        A list of Sequence objects and None for headers.

    Raises:
        AlignmentIOError: If a dictionary value is not a list.
    """
    headers = None
    sequence_list = []
    for name, values in values_dict.items():
        if not isinstance(values, list):
            raise AlignmentIOError('Values must be a list')
        seq = Sequence(name=name)
        seq.set_cont_values(values)
        sequence_list.append(seq)
    return sequence_list, headers


# .............................................................................
def get_character_matrix_from_sequences_list(sequences, var_headers=None):
    """Converts a list of sequences into a character matrix.

    Args:
        sequences (:obj:`list` of :obj:`Sequence`): A list of Sequence objects
            to be converted.
        var_headers (:obj:`list` of headers, optional): If provided, uses these
            as variable headers for the columns in the matrix.

    Returns:
        Matrix: A matrix of sequence data.
    """
    if var_headers is not None:
        col_headers = var_headers
    else:
        col_headers = [
            'Column {}'.format(i) for i in range(
                len(sequences[0].cont_values))]
    data = np.zeros((len(sequences), len(col_headers)), dtype=float)
    row_headers = []
    i = 0
    for seq in sequences:
        row_headers.append(seq.name)
        data[i] = np.array(seq.cont_values)
        i += 1
    return Matrix(data, headers={'0': row_headers, '1': col_headers})


# .............................................................................
def load_alignment_from_filename(filename):
    """Attempts to load an alignment from a file path by guessing schema

    Args:
        filename (str): The file location containing the alignment

    Raises:
        RuntimeError: Raised with the method needed to load the alignment
            cannot be determined.

    Returns:
        tuple: Containing a list of sequences and headers
    """
    _, ext = os.path.splitext(filename)
    if ext == '.csv':
        load_method = read_csv_alignment_flo
    elif ext == '.json':
        load_method = read_json_alignment_flo
    elif ext == '.phylip':
        load_method = read_phylip_alignment_flo
    elif ext == '.tbl':
        load_method = read_table_alignment_flo
    else:
        raise RuntimeError(
            'Cannot determine load method for {} -- extension {}'.format(
                filename, ext))
    with open(filename) as align_file:
        ret = load_method(align_file)
    try:
        sequences, headers = ret
    except:
        sequences = ret
        headers = None
    return sequences, headers


# .............................................................................
def read_csv_alignment_flo(csv_flo):
    """Reads a CSV file-like object and return a list of sequences and headers.

    Args:
        csv_flo (file-like): A file-like object with CSV alignment data.

    Returns:
        A list of Sequence objects and headers.

    Raises:
        AlignmentIOError: If the number of columns is inconsistent across the
            sequences.
    """
    headers = None
    sequence_list = []

    has_header = csv.Sniffer().has_header(csv_flo.readline())
    csv_flo.seek(0)

    num_parts = None
    for line in csv_flo:
        parts = line.strip().split(',')
        if num_parts is None:
            num_parts = len(parts)
        else:
            if len(parts) != num_parts:
                raise AlignmentIOError('Number of columns is inconsistent')
        if has_header and headers is None:
            headers = parts[1:]
        else:
            name = parts[0]
            vals = [float(i) for i in parts[1:]]
            seq = Sequence(name=name)
            seq.set_cont_values(vals)
            sequence_list.append(seq)
    return sequence_list, headers


# .............................................................................
def read_json_alignment_flo(json_flo):
    """Read a JSON file-like object and return a list of sequences and headers.

    Args:
        json_flo (file-like): A file-like object with JSON alignment data.

    Note:
        * File should have structure::

            {
                "headers" : [{header_names}],
                "values" : [
                    {
                        "name" : "{taxon_name}",
                        "values" : [{values}]
                    }
                ]
            }

    Returns:
        A list of Sequence objects and headers.

    Raises:
        AlignmentIOError: If headers are provided but they are not a list.
    """
    json_vals = json.load(json_flo)

    if 'headers' in json_vals.keys():
        headers = json_vals['headers']
        if not isinstance(headers, list):
            raise AlignmentIOError(
                'If headers are provided, they must be a list')
    else:
        headers = None

    sequence_list = []
    for val_dict in json_vals['values']:
        name = val_dict['name']
        vals = [float(v) for v in val_dict['values']]
        seq = Sequence(name=name)
        seq.set_cont_values(vals)
        sequence_list.append(seq)
    return sequence_list, headers


# .............................................................................
def read_phylip_alignment_flo(phylip_flo):
    """Reads a phylip alignment file-like object and return the sequences.

    Args:
        phylip_flo (file-like): The phylip file-like object.

    Note:
        * We assume that the phylip files are extended and not strict (in terms
            of how many characters for taxon names).
        * The phylip file is in the format::
            numoftaxa numofsites
            seqlabel sequence
            seqlabel sequence

    Returns:
        A list of Sequence objects.

    Raises:
        AlignmentIOError: If there is a problem creating sequences.
    """
    seqlist = []
    # first line is the number of taxa and num of sites
    # we don't really even need to read this line,
    # so let's just skip it
    i = phylip_flo.readline()
    for i in phylip_flo:
        try:
            if len(i) > 2:
                spls = i.strip().split()
                name = spls[0].strip()
                seq = spls[1].strip()
                tseq = Sequence(name=name, seq=seq)
                seqlist.append(tseq)
        except Exception as e:
            raise AlignmentIOError(str(e))
    return seqlist


# .............................................................................
def read_table_alignment_flo(table_flo):
    """Reads a table from a file-like object.

    Args:
        table_flo (file-like): A file-like object containing table data.

    Returns:
        A list of Sequence objects.

    Raises:
        AlignmentIOError: If there is a problem creating sequences.
    """
    seqlist = []
    for i in table_flo:
        if len(i) > 2:
            try:
                spls = i.strip().split("\t")
                name = spls[0].strip()
                seq = spls[1].strip().split(" ")
                seq = [float(j) for j in seq]
                tseq = Sequence(name=name)
                tseq.set_cont_values(seq)
                seqlist.append(tseq)
            except Exception as e:
                raise AlignmentIOError(str(e))
    return seqlist

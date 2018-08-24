"""
@summary: This module contains functions for reading data files in various 
             formats.  This is for reading some common alignment files.  Each
             function returns a list of sequences
@note: If you want to read an alignment in fasta format, just use the fasta
          reader in the seq_reader.py
@todo: FASTA reader?
@todo: Use file-like objects
"""
import csv
import json
import os
import sys

from ancestral_reconstruction.helpers.sequence import Sequence

# .............................................................................
def create_sequence_list_from_dict(values_dict):
    """
    @summary: Creates a list of sequences from a dictionary
    @note: The dictionary should have structure:
              {
               "{taxon_name}" : [{values}]
              }
    """
    headers = None
    sequence_list = []
    for name, values in values_dict.iteritems():
       seq = Sequence(name=name)
       seq.set_cont_values(values)
       sequence_list.append(seq)
    return sequence_list, headers

# .............................................................................
def read_csv_alignment_flo(csv_flo):
    """
    @summary: Read a CSV file-like object and return a list of sequences and 
                 headers
    @param csv_flo: A file-like object with CSV alignment data
    """
    headers = None
    seqence_list = []
    
    has_header = csv.Sniffer().has_header(csv_flo.read(5000))
    csv_flo.seek(0)
    
    for line in csv_flo:
        parts = line.strip().split(',')
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
def read_json_alignment_file(json_filename):
    """
    @summary: Read a JSON file and return a list of sequences and headers
    @param json_filename: A file location of a JSON file with alignment data
    @note: File should have structure: 
               {
                "headers" : [{header_names}],
                "values" : [
                            {
                             "name" : "{taxon_name}",
                             "values" : [{values}]
                            }
                           ]
               }
    """
    with open(json_filename) as json_file:
        json_vals = json.load(json_file)
        
    if json_vals.has_key('headers'):
        headers = json_vals['headers']
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
def read_phylip_file(infilename):
    """
    @summary: This will read a phylip alignment file and return the list of 
                 sequences contained
    @param infilename: The phylip input file
    @note: We assume that the phylip files are extended and not strict (in 
              terms of how many characters for taxon names)
    @note: The phylip file is in the format:
              numoftaxa numofsites
              seqlabel sequence
              seqlabel sequence
    """
    infile = open(infilename,"r")
    seqlist = []
    # first line is the number of taxa and num of sites
    # we don't really even need to read this line, 
    # so let's just skip it
    i = infile.readline()
    for i in infile:
        if len(i) > 2:
            spls = i.strip().split()
            name = spls[0].strip()
            seq = spls[1].strip()
            tseq = Sequence(name=name,seq=seq)
            seqlist.append(tseq)
    infile.close()
    return seqlist

# .............................................................................
def read_phylip_cont_file(infile):
    """
    @summary: This will read a phylip alignment file with continuous characters 
                 and return the list of seqs
    @note: We assume that the phylip files are extended and not strict (in 
              terms of what type and how much white space and how many 
              characters for taxon names)
    @note: The phylip file is in the format:
              numoftaxa numofsites
              seqlabel contvalue contvalue
              seqlabel contvalue contvalue
    """
    seqlist = []
    # first line is the number of taxa and num of sites
    # we don't really even need to read this line, 
    # so let's just skip it
    i = infile.readline()
    for i in infile:
        if len(i) > 2:
            spls = i.strip().split("\t")
            name = spls[0].strip()
            seq = spls[1].strip().split(" ")
            seq = [float(j) for j in seq]
            tseq = Sequence(name=name)
            tseq.set_cont_values(seq)
            seqlist.append(tseq)
    return seqlist

# .............................................................................
def read_table_cont_file(infile):
    seqlist = []
    for i in infile:
        if len(i) > 2:
            spls = i.strip().split("\t")
            name = spls[0].strip()
            seq = spls[1].strip().split(" ")
            seq = [float(j) for j in seq]
            tseq = Sequence(name=name)
            tseq.set_cont_values(seq)
            seqlist.append(tseq)
    return seqlist




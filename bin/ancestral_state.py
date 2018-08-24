#!/usr/bin/env python
"""
@summary:
@todo: Rename without .py extension
@todo: Constants
@todo: Clean up help
"""
import argparse
import dendropy
import os

from ancestral_reconstruction.helpers.data_readers import (read_csv_alignment_file,
                              read_json_alignment_file, read_phylip_cont_file)
from ancestral_reconstruction.analysis.anc_dp import calculate_continuous_ancestral_states

# .............................................................................
if __name__ == '__main__':
   parser = argparse.ArgumentParser()

   parser.add_argument('in_tree_filename', type=str, 
                       help='Path to the tree file')
   parser.add_argument('in_tree_schema', type=str, 
                       help='The format of the tree', 
                       choices=['newick', 'nexml', 'nexus'])
   
   parser.add_argument('data_filename', type=str, 
                       help='Path to file with character state data')
   parser.add_argument('data_format', type=str, 
                       help='The format of the character data', 
                       choices=['csv', 'json', 'phylip'])
   parser.add_argument('-ch', '--column_headers', type=str, 
                  help='Path to file containing column headers, one per line')
   
   parser.add_argument('out_tree_filename', type=str, 
                       help='Path for new output tree')
   parser.add_argument('out_tree_schema', type=str, 
                       help='The format of the tree', 
                       choices=['newick', 'nexml', 'nexus'])

   args = parser.parse_args()
   
   # Check that input files exist
   if not os.path.exists(args.in_tree_filename):
      raise IOError('Input tree {} does not exist'.format(
                                                      args.in_tree_filename))
   if not os.path.exists(args.data_filename):
      raise IOError('Input data file {} does not exist'.format(
                                                      args.data_filename))
   
   # Read tree
   tree = dendropy.Tree.get(path=args.in_tree_filename, 
                                                   schema=args.in_tree_schema)

   # Read data
   if args.data_format == 'csv':
      sequences, headers = read_csv_alignment_file(args.data_filename)
   elif args.data_format == 'json':
      sequences, headers = read_json_alginment_file(args.data_filename)
   elif args.data_format == 'phylip':
      with open(args.data_filename) as in_file:
         sequences = read_phylip_cont_file(in_file)
      headers = None
   else:
      raise Exception, 'Unknown data format: {}'.format(args.data_format)
   
   # Run analysis
   calculate_continuous_ancestral_states(tree, sequences)
   
   # Write outputs
   with open(args.out_tree_filename, 'w') as out_file:
      out_file.write(tree.as_string(schema=args.out_tree_schema))


#    seqs = aln_reader.read_phylip_cont_file(open(sys.argv[2],"r"))
#    match_tips_and_cont_values(tree,seqs)
#    sqch = calc_cont_anc_states(tree)
#    outfile = open("contanc_dp.tre","w")
#    outfile.write(tree.as_string(schema="newick"))
#    outfile.close()

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

from ancestral_reconstruction.helpers import data_readers
from ancestral_reconstruction.analysis.anc_dp import calculate_continuous_ancestral_states

DESCRIPTION = """\
Generate ancestral state estimations for the nodes in a tree based on the 
variable values at the tips.
"""

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
    parser.add_argument(
        'out_characters_filename', type=str,
        help='A file location to write the character data as CSV')

    args = parser.parse_args()

    # Check that input files exist
    if not os.path.exists(args.in_tree_filename):
        raise IOError('Input tree {} does not exist'.format(
                                                      args.in_tree_filename))
    if not os.path.exists(args.data_filename):
        raise IOError('Input data file {} does not exist'.format(
                                                      args.data_filename))
    # Read data
    if args.data_format == 'csv':
        with open(args.data_filename) as in_file:
            sequences, headers = data_readers.read_csv_alignment_flo(in_file)
    elif args.data_format == 'json':
        with open(args.data_filename) as in_file:
            sequences, headers = data_readers.read_json_alginment_flo(in_file)
    elif args.data_format == 'phylip':
        with open(args.data_filename) as in_file:
            sequences = data_reders.read_phylip_alignment_flo(in_file)
        headers = None
    elif args.data_format == 'table':
        with open(args.data_filename) as in_file:
            sequences = data_readers.read_table_alignment_flo(in_file)
        headers = None
    else:
        raise Exception, 'Unknown data format: {}'.format(args.data_format)
   
    # Read tree
    tree = dendropy.Tree.get(path=args.in_tree_filename, 
                             schema=args.in_tree_schema)

    # Run analysis
    char_mtx = data_readers.get_character_matrix_from_sequences_list(sequences)
    (out_tree, anc_state) = calculate_continuous_ancestral_states(
        tree, char_mtx)

    # Write outputs
    with open(args.out_tree_filename, 'w') as out_file:
        out_file.write(out_tree.as_string(schema=args.out_tree_schema))

    with open(args.out_characters_filename, 'w') as out_file:
        anc_state.write_csv(out_file)

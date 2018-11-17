#!/usr/bin/env python
"""Tool for performing ancestral distribution computations

Todo:
    * Rename without .py extension
    * Constants
    * Clean up help
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

from ancestral_reconstruction.helpers import data_readers
from ancestral_reconstruction.analysis import anc_dp
from ancestral_reconstruction.lm_objects.tree import TreeWrapper

DESCRIPTION = """\
Generates ancestral distribution estimations based on the environmental 
distributions at the tips of the tree"""


# .............................................................................
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument('in_tree_filename', type=str, 
                        help='Path to the tree file')
    parser.add_argument('in_tree_schema', type=str, 
                        help='The format of the tree', 
                        choices=['newick', 'nexml', 'nexus'])
   
    parser.add_argument('data_filename', type=str, 
                        help='Path to file with character state data')
    parser.add_argument('data_format', type=str, 
                        help='The format of the character data', 
                        choices=['csv', 'json', 'phylip', 'table'])
    # out directory
    parser.add_argument('-o', '--output_directory', type=str, 
                        help='If provided, write the rates and plots here')

    # print tree
    # parser.add_argument('--print_tree', action='store_true', required=False,
    #                     help='Show the tree (requires ete3)')
    args = parser.parse_args()

    # Check that input files exist
    if not os.path.exists(args.in_tree_filename):
       raise IOError('Input tree {} does not exist'.format(
                                                       args.in_tree_filename))
    if not os.path.exists(args.data_filename):
       raise IOError('Input data file {} does not exist'.format(
                                                      args.data_filename))

    # Read the tree
    tree = TreeWrapper.get(path=args.in_tree_filename, 
                           schema=args.in_tree_schema)

    # Read data
    if args.data_format == 'csv':
        sequences, headers = data_readers.read_csv_alignment_flo(
            args.data_filename)
    elif args.data_format == 'json':
        sequences, headers = data_readers.read_json_alginment_flo(
            args.data_filename)
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

    # Get character matrix
    char_mtx = data_readers.get_character_matrix_from_sequences_list(
        sequences, var_headers=headers)
    # Run analysis
    tree, results = anc_dp.calculate_ancestral_distributions(tree, char_mtx)

    # Write outputs
    if args.output_directory is not None:
        # Write plots
        node_names = results.get_row_headers()
        num_cats = char_mtx.data.shape[1]
        low = 0
        high = num_cats
        
        x_grid = np.linspace(low, high, num_cats)
        for i in range(len(node_names)):
            plt.figure(figsize=(6, 4))
            plt.plot(x_grid, results.data[i, :, 0])
            plt.fill_between(x_grid, 0, results.data[i, :, 0], alpha=0.05)
            if np.any(results.data[i, :, 1] > 0.0):
                high_vals = results.data[i, :, 0] + results.data[i, :, 1]
                low_vals = results.data[i, :, 0] - results.data[i, :, 1]
                plt.plot(x_grid, high_vals, '--', alpha=0.55,)
                plt.plot(x_grid, low_vals, '--', alpha=0.55,)
            plt.grid(True)
            plt.savefig(os.path.join(args.output_directory, 
                                     '{}.png'.format(node_names[i])))
            plt.close()

        # Add rates to output
        # TODO: Calculate rates
        # plt.plot(x_grid, rates)
        # plt.ylabel('estimated rate')
        # plt.xlabel('state space')
        # plt.savefig(os.path.join(args.output_directory, 'rate_estimate.png'))

        # Write labeled tree
        tree.write(path=os.path.join(args.output_directory, 
                                     'labeled_tree.tre'), schema='newick')

        # Write results matrix (as csv)
        distributions_filename = os.path.join(args.output_directory, 
                                              'distributions.csv')
        with open(distributions_filename, 'w') as distributions_file:
            results.write_csv(distributions_file)
    
    # If we should print the tree, do it
    # if args.print_tree:
    #     print_tree_to_file(tree, args.output_directory)

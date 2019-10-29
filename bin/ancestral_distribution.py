#!/usr/bin/env python
"""Tool for performing ancestral distribution computations.

Todo:
    * Rename without .py extension.
    * Constants.
    * Clean up help.
"""
import argparse
import os

from lmpy import TreeWrapper

from analyses.ancestral_state import anc_dp
import analyses.tools.annotators as annotators
import analyses.tools.plots as tree_plots
from analyses.helpers import data_readers

DESCRIPTION = """\
Generates ancestral distribution estimations based on the environmental
 distributions at the tips of the tree"""


# .............................................................................
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        'in_tree_filename', type=str, help='Path to the tree file')
    parser.add_argument(
        'in_tree_schema', type=str, help='The format of the tree',
        choices=['newick', 'nexml', 'nexus'])

    parser.add_argument(
        'data_filename', type=str,
        help='Path to file with character state data')
    parser.add_argument(
        'data_format', type=str, help='The format of the character data',
        choices=['csv', 'json', 'phylip', 'table'])

    # Outputs
    # Annotated tree or trees
    # Plots
    # Matrix csv
    parser.add_argument(
        'out_tree_filename', type=str,
        help='Path to write the resulting annotated tree')
    parser.add_argument(
        'out_tree_schema', type=str,
        help='The format to use when writing the tree',
        choices=['newick', 'nexml', 'nexus'])
    parser.add_argument(
        '-l', '--annotate_labels', type=str,
        help='If provided, annotate the tree labels with this data column')
    parser.add_argument(
        '-p', '--plot_directory', type=str,
        help='If provided, write distribution plots to this directory')
    parser.add_argument(
        '-c', '--out_csv_filename', type=str,
        help='If provided, write the output character matrix CSV '
             'to this file location')

    args = parser.parse_args()

    # Check that input files exist
    if not os.path.exists(args.in_tree_filename):
        raise IOError(
            'Input tree {} does not exist'.format(args.in_tree_filename))
    if not os.path.exists(args.data_filename):
        raise IOError(
            'Input data file {} does not exist'.format(args.data_filename))

    # Read the tree
    tree = TreeWrapper.get(
        path=args.in_tree_filename, schema=args.in_tree_schema)

    # Read data
    if args.data_format == 'csv':
        with open(args.data_filename) as in_file:
            sequences, headers = data_readers.read_csv_alignment_flo(
                in_file)
    elif args.data_format == 'json':
        with open(args.data_filename) as in_file:
            sequences, headers = data_readers.read_json_alginment_flo(
                inf_file)
    elif args.data_format == 'phylip':
        with open(args.data_filename) as in_file:
            sequences = data_reders.read_phylip_alignment_flo(in_file)
        headers = None
    elif args.data_format == 'table':
        with open(args.data_filename) as in_file:
            sequences = data_readers.read_table_alignment_flo(in_file)
        headers = None
    else:
        raise Exception('Unknown data format: {}'.format(args.data_format))

    # Get the label annotation column, or None
    label_column = None
    if args.annotate_labels is not None:
        try:
            # Try looking for the string
            label_column = headers.index(args.annotate_labels)
        except:
            try:
                # Treat it as an integer
                label_column = int(args.annotate_labels)
            except:
                raise Exception(
                    'Could not find column to use for labels.  '
                    'Check the name to make sure it matches or use column'
                    ' index.')

    # Get character matrix
    char_mtx = data_readers.get_character_matrix_from_sequences_list(
        sequences, var_headers=headers)
    # Run analysis
    tree, results = anc_dp.calculate_ancestral_distributions(tree, char_mtx)

    # Should we annotate the tree labels?
    if label_column is not None:
        annotators.annotate_tree_with_label(
            tree, results, label_column=label_column)
    else:
        # Annotate tree
        annotators.add_all_annotations(tree, results, update=True)

    # Write the tree
    tree.write(path=args.out_tree_filename, schema=args.out_tree_schema)

    # CSV
    if args.out_csv_filename is not None:
        with open(args.out_csv_filename, 'w') as out_csv_f:
            results.write_csv(out_csv_f)
    # Plots
    if args.plot_directory is not None:
        tree_plots.create_distribution_plots(
            tree, results, args.plot_directory)

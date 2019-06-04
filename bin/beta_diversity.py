#!/usr/bin/env python
"""Tool for computing beta diversity and phylo beta diversity.

ARGS:
1) file location of the PAM (presence absence matrix)
2) file location of the tree associsated with PAM
3) phylo beta diversity family to use
4) file location to write output
5) (optional) the number of permuatations to perform 
6) (optional) alpha value to determine significance (defaul = 0.05)

Todo:
    * Rename without .py extension.
    * Constants.
    * Clean up help.
"""
import argparse
import os
import numpy as np

# import data_readers
# from tree import TreeWrapper
# import phylo_beta_diversity


from analyses.lm_objects.tree import TreeWrapper
from analyses.helpers import data_readers
from analyses.phylo_beta_diversity import phylo_beta_diversity

DESCRIPTION = """\
Computes phylogenetic & ecological beta diversity components for Sorensen and Jaccard Indices."""


# .............................................................................
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        'in_tree_filename', type=str, help='Path to the tree file')
    parser.add_argument(
        'in_tree_schema', type=str, help='The format of the tree',
        choices=['newick', 'nexml', 'nexus'])

    parser.add_argument(
        'pam_filename', type=str,
        help='Path to file with presence/ absence data (PAM)')
    parser.add_argument(
        'data_format', type=str, help='The format of the PAM',
        choices=['csv', 'json', 'phylip', 'table'])

    # Outputs
    # Annotated tree or trees
    # Plots
    # Matrix csv

# Might use this for nodal beta diversity 

    # parser.add_argument(
    #     'out_tree_filename', type=str,
    #     help='Path to write the resulting annotated tree')
    # parser.add_argument(
    #     'out_tree_schema', type=str,
    #     help='The format to use when writing the tree',
    #     choices=['newick', 'nexml', 'nexus'])
    # parser.add_argument(
    #     '-l', '--annotate_labels', type=str,
    #     help='If provided, annotate the tree labels with this data column')
    # parser.add_argument(
    #     '-p', '--plot_directory', type=str,
    #     help='If provided, write distribution plots to this directory')

    parser.add_argument(
        'family_name', type=str,
        help='Beta diversity family metric to calculate',
        choices=['sorensen', 'jaccard'])

    parser.add_argument(
        'out_foldername', type=str,
        help='Write the output of beta diversity calculations to this folder')

    parser.add_argument(
        '-n', '--number_permutations', 
        help='The number of permuatations to calculate')

    parser.add_argument(
        '-a', '--alpha',
        help='The alpha value to determine significance') # make default = 0.05

    args = parser.parse_args()

    # Number of iterations
    if args.number_permutations == None:
        nrand = 10
    else:
        nrand = int(args.number_permutations)
    print(nrand)

    # Check that input files exist
    if not os.path.exists(args.in_tree_filename):
        raise IOError(
            'Input tree {} does not exist'.format(args.in_tree_filename))
    if not os.path.exists(args.pam_filename):
        raise IOError(
            'Input data file {} does not exist'.format(args.pam_filename))

    # Read data
    if args.data_format == 'csv':
        with open(args.pam_filename) as in_file:
            sequences, headers = data_readers.read_csv_alignment_flo(
                in_file)
    elif args.data_format == 'json':
        with open(args.pam_filename) as in_file:
            sequences, headers = data_readers.read_json_alignment_flo(
                in_file)
    elif args.data_format == 'phylip':
        with open(args.pam_filename) as in_file:
            sequences = data_readers.read_phylip_alignment_flo(in_file) 
        headers = None
    elif args.data_format == 'table':
        with open(args.pam_filename) as in_file:
            sequences = data_readers.read_table_alignment_flo(in_file)
        headers = None
    else:
        raise Exception('Unknown data format: {}'.format(args.data_format))

    # Get the label annotation column, or None
    # label_column = None
    # if args.annotate_labels is not None:
    #     try:
    #         # Try looking for the string
    #         label_column = headers.index(args.annotate_labels)
    #     except:
    #         try:
    #             # Treat it as an integer
    #             label_column = int(args.annotate_labels)
    #         except:
    #             raise Exception(
    #                 'Could not find column to use for labels.  '
    #                 'Check the name to make sure it matches or use column'
    #                 ' index.')

    # Convert data to PAM format
    pam = data_readers.get_character_matrix_from_sequences_list(sequences, headers)

    # Read the tree
    tree = TreeWrapper.get(
        path=args.in_tree_filename, schema=args.in_tree_schema)

    # Run analysis
    if args.family_name == 'jaccard':
        results = phylo_beta_diversity.calculate_phylo_beta_diversity_jaccard(pam, tree)
        res_names = ['beta_jtu', 'phylo_beta_jtu', 'beta_jne', 'phylo_beta_jne', 'beta_jac', 'phylo_beta_jac']
    elif args.family_name == 'sorensen':
        results = phylo_beta_diversity.calculate_phylo_beta_diversity_sorensen(pam, tree)
        res_names = ['beta_sim', 'phylo_beta_sim', 'beta_sne', 'phylo_beta_sne', 'beta_sor', 'phylo_beta_sor']
    else:
        raise Exception('Could not find family name')

    # Should we annotate the tree labels?
    # if label_column is not None:
    #     annotators.annotate_tree_with_label(
    #         tree, results, label_column=label_column)
    # else:
    #     # Annotate tree
    #     annotators.add_all_annotations(tree, results, update=True)

    # Write the tree
    # tree.write(path=args.out_tree_filename, schema=args.out_tree_schema)

    # Write results to folder
    if not os.path.exists(args.out_foldername):
        os.makedirs(args.out_foldername)
        
    for table in range(len(results)):
        # print table, res_names[table], "\n", results[table].data, "\n"
        with open(os.path.join(args.out_foldername, '{}.csv'.format(res_names[table])), 'w') as out_csv_f:
            results[table].write_csv(out_csv_f)
       
    # Test randomization:
    # List of lists: 0:JTU; 1:JNE; 2:JAC.
    print (res_names[5])
    x = phylo_beta_diversity.calc_phylo_jac_distr(pam, tree, obs = results[5], metric = "jne", nrand = nrand)
    print(x.data)
    # print x[1].data[0]#[1][0]
    # print x[1].data[1][0]
    # for i in range(len(x)):
    #     print(x[i].data)

    # # Plots
    # if args.plot_directory is not None:
    #     tree_plots.create_distribution_plots(
    #         tree, results, args.plot_directory)

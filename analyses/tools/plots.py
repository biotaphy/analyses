"""Module containing tools for plotting tree nodes.

Todo:
    * Determine best path for this module.
    * Determine a good way to pass plot options.
"""
import os

import matplotlib.pyplot as plt
import numpy as np


plt.switch_backend('Agg')


# .............................................................................
def create_distribution_plots(lm_tree, node_matrix, output_directory):
    """Creates distribution plots for each of the nodes in the tree.

    Args:
        lm_tree (TreeWrapper): A Lifemapper tree object.
        node_matrix (Matrix): A Lifemapper Matrix object with rows matching the
            nodes in the provided tree.
        output_directory (str): A directory where the output plots should be
            written.
    """
    # If the output directory does not exist, create it
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    node_names = node_matrix.get_row_headers()
    num_cats = node_matrix.shape[1]
    low = 0
    high = num_cats

    x_grid = np.linspace(low, high, num_cats)
    for i in range(len(node_names)):
        # TODO(CJ) : Variable figure size
        plt.figure(figsize=(6, 4))
        # TODO(CJ) : Variable alpha
        plt.plot(x_grid, node_matrix[i, :, 0], alpha=0.05)
        if np.any(node_matrix[i, :, 1] > 0.0):
            # TODO(CJ) : Options for high and low value lines
            high_vals = node_matrix[i, :, 0] + node_matrix[i, :, 1]
            low_vals = node_matrix[i, :, 0] - node_matrix[i, :, 1]
            plt.plot(x_grid, high_vals, '--', alpha=0.55,)
            plt.plot(x_grid, low_vals, '--', alpha=0.55,)
        # TODO(CJ) : Does this print to screen?  Remove it?
        plt.grid(True)
        # TODO: Calculate rates
        # plt.plot(x_grid, rates)
        # plt.ylabel('estimated rate')
        # plt.xlabel('state space')
        # plt.savefig(os.path.join(args.output_directory, 'rate_estimate.png'))
        # TODO(CJ) : File name options?
        plt.savefig(os.path.join(
            output_directory, '{}.png'.format(node_names[i])))
        plt.close()

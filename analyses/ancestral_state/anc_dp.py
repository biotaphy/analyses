"""Module containing code for calculating ancestral states.
"""
import math
import sys

import numpy as np
import scipy.linalg as la

from lmpy import Matrix, TreeWrapper


# .............................................................................
def _get_node_label(node):
    """Returns the node label or taxon label if node is a tip.

    Args:
        node (Node): A tree node to get the label for.

    Returns:
        str: The node label or taxon label.
    """
    if node.label is not None:
        return node.label
    else:
        return node.taxon.label


# .............................................................................
def calculate_ancestral_distributions(tree, char_mtx):
    """Calculates ancestral distributions.

    Args:
        tree (Tree): A dendropy tree or TreeWrapper object.
        char_mtx (Matrix): A Matrix object with character information.  Each
            row should represent a tip in the tree and each column should be a
            bin to calculate ancestral distribution.

    Returns:
        A matrix of character data with the following dimensions:
            * rows: nodes / tips in the tree
            * columns: character variables
            * depth: first is the calculated value, second layer is standard
                error if desired
    """
    return calculate_continuous_ancestral_states(tree, char_mtx,
                                                 sum_to_one=True,
                                                 calc_std_err=True)


# .............................................................................
def calculate_continuous_ancestral_states(tree, char_mtx, sum_to_one=False,
                                          calc_std_err=False):
    """Calculates the continuous ancestral states for the nodes in a tree.

    Args:
        tree (Tree): A dendropy tree or TreeWrapper object.
        char_mtx (Matrix): A Matrix object with character information.  Each
            row should represent a tip in the tree and each column should be a
            variable to calculate ancestral state for.
        calc_std_err (:obj:`bool`, optional): If True, calculate standard error
            for each variable.  Defaults to False.
        sum_to_one (:obj:`bool`, optional): If True, standardize the character
            matrix so that the values in a row sum to one. Defaults to False.

    Returns:
        A matrix of character data with the following dimensions:
            * rows: nodes / tips in the tree
            * columns: character variables
            * depth: first is the calculated value, second layer is standard
                error if desired

    Todo:
        * Add function for consistent label handling.
    """
    # Wrap tree if dendropy tree
    if not isinstance(tree, TreeWrapper):
        tree = TreeWrapper.from_base_tree(tree)

    # Assign labels to nodes that don't have them
    tree.add_node_labels()

    # Synchronize tree and character data
    # Prune tree
    prune_taxa = []
    keep_taxon_labels = []
    init_row_headers = char_mtx.get_row_headers()
    for taxon in tree.taxon_namespace:
        label = taxon.label.replace(' ', '_')
        if label not in init_row_headers:
            prune_taxa.append(taxon)
            print(
                'Could not find {} in character matrix, pruning'.format(label))
        else:
            keep_taxon_labels.append(label)

    if len(keep_taxon_labels) == 0:
        raise Exception(
            'None of the tree tips were found in the character data')

    tree.prune_taxa(prune_taxa)
    tree.purge_taxon_namespace()

    # Prune character data
    keep_rows = []
    i = 0
    for label in init_row_headers:
        if label in keep_taxon_labels:
            keep_rows.append(i)
        else:
            print('Could not find {} in tree tips, pruning'.format(label))
        i += 1
    char_mtx = char_mtx.slice(keep_rows)

    # Standardize character matrix if requested
    tip_count, num_vars = char_mtx.shape
    if sum_to_one:
        for i in range(tip_count):
            sc = float(1.0) / np.sum(char_mtx[i])
            for j in range(num_vars):
                char_mtx[i, j] *= sc

    # Initialize data matrix
    num_nodes = len(tree.nodes())
    data_shape = (num_nodes, num_vars, 2 if calc_std_err else 1)
    data = np.zeros(data_shape, dtype=float)

    # Initialize headers
    row_headers = []

    tip_col_headers = char_mtx.get_column_headers()
    tip_row_headers = char_mtx.get_row_headers()
    tip_lookup = dict([
        (tip_row_headers[i].replace('_', ' '), i) for i in range(tip_count)])

    # Get the number of internal nodes in the tree
    internal_node_count = num_nodes - tip_count
    # Loop through the tree and set the matrix index for each node
    # Also set data values
    node_headers = []
    node_i = tip_count
    tip_i = 0
    node_index_lookup = {}
    for node in tree.nodes():
        label = _get_node_label(node)
        if len(node.child_nodes()) == 0:
            # Tip
            node_index_lookup[label] = tip_i
            row_headers.append(label)
            data[tip_i, :, 0] = char_mtx[tip_lookup[label]]
            tip_i += 1
        else:
            node_index_lookup[label] = node_i
            node_headers.append(label)
            # Internal node
            data[node_i, :, 0] = np.zeros((1, num_vars), dtype=float)
            node_i += 1

    # Row headers should be extended with node headers
    row_headers.extend(node_headers)

    # For each variable
    for x in range(num_vars):
        # Compute the ML estimate of the root
        full_mcp = np.zeros((internal_node_count, internal_node_count),
                            dtype=float)
        full_vcp = np.zeros(internal_node_count, dtype=float)

        for k in tree.postorder_edge_iter():
            i = k.head_node
            if len(i.child_nodes()) != 0:
                node_num_i = node_index_lookup[_get_node_label(i)] - tip_count
                for j in i.child_nodes():
                    tbl = 2./j.edge_length
                    full_mcp[node_num_i][node_num_i] += tbl
                    node_num_j = node_index_lookup[_get_node_label(j)]

                    if len(j.child_nodes()) == 0:
                        full_vcp[node_num_i] += (data[node_num_j, x, 0] * tbl)
                    else:
                        node_num_j -= tip_count
                        full_mcp[node_num_i][node_num_j] -= tbl
                        full_mcp[node_num_j][node_num_i] -= tbl
                        full_mcp[node_num_j][node_num_j] += tbl

        b = la.cho_factor(full_mcp)

        # these are the ML estimates for the ancestral states
        ml_est = la.cho_solve(b, full_vcp)
        sos = 0
        for k in tree.postorder_edge_iter():
            i = k.head_node
            node_num_i = node_index_lookup[_get_node_label(i)]
            if len(i.child_nodes()) != 0:
                data[node_num_i, x, 0] = ml_est[node_num_i - tip_count]

                if calc_std_err:
                    for j in i.child_nodes():
                        node_num_j = node_index_lookup[_get_node_label(j)]
                        temp = data[node_num_i, x, 0] - data[node_num_j, x, 0]
                        sos += temp * temp / j.edge_length

                    # nni is node_num_i adjusted for only nodes
                    nni = node_num_i - tip_count
                    qpq = full_mcp[nni][nni]
                    tm1 = np.delete(full_mcp, (nni), axis=0)
                    tm = np.delete(tm1, (nni), axis=1)
                    b = la.cho_factor(tm)
                    sol = la.cho_solve(b, tm1[:, nni])
                    temp_std_err = qpq - np.inner(tm1[:, nni], sol)
                    data[node_num_i, x, 1] = math.sqrt(2.0 * sos / (
                        (internal_node_count - 1) * temp_std_err))

    depth_headers = ['maximum_likelihood']
    if calc_std_err:
        depth_headers.append('standard_error')

    mtx_headers = {
        '0': row_headers,
        '1': tip_col_headers,
        '2': depth_headers
        }
    return tree, Matrix(data, headers=mtx_headers)

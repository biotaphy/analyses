"""
@summary: Module containing code for calculating ancestral states
"""
import dendropy as dp
import math
import numpy as np
import scipy.linalg as la
import sys

from ancestral_reconstruction.lm_objects.matrix import Matrix
from ancestral_reconstruction.lm_objects.tree import TreeWrapper


# .............................................................................
def _get_node_label(node):
    """
    @summary: Return the node label or the taxon label if node is a tip
    @param node: A tree node to get the label for
    """
    if node.label is not None:
        return node.label
    else:
        return node.taxon.label


# .............................................................................
def calculate_continuous_ancestral_states(tree, char_mtx, calc_std_err=False):
    """
    @summary: Calculates the continuous ancestral states for the nodes in a
                tree
    @param tree: A dendropy tree or TreeWrapper object
    @param char_mtx: A Matrix object with character information.  Each row
                        should represent a tip in the tree and each column
                        should be a variable to calculate ancestral state for
    @param calc_std_err: If True, calculate standard error for each variable
    @return: A matrix of character data with the following dimensions:
                rows - nodes / tips in the tree
                columns - character variables
                depth - first is the calculated value, second layer is
                            standard error if desired
    @todo: Standard error computations
    """
    # Wrap tree if dendropy tree
    if not isinstance(tree, TreeWrapper):
        tree = TreeWrapper.from_base_tree(tree)

    # Initialize data matrix
    num_nodes = len(tree.nodes())
    num_vars = char_mtx.data.shape[1]
    data_shape = (num_nodes, num_vars, 2 if calc_std_err else 1)
    data = np.zeros(data_shape, dtype=float)

    # Initialize headers
    row_headers = []

    # Assign labels to nodes that don't have them
    tree.add_node_labels()

    # TODO: Ensure that tree tips match provided character data and fill in
    #    data matrix for tips
    tip_col_headers = char_mtx.get_column_headers()
    tip_row_headers = char_mtx.get_row_headers()
    tip_count = len(tip_row_headers)
    tip_lookup = dict([(tip_row_headers[i], i) for i in range(tip_count)])

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
            if label not in tip_lookup.keys():
                # Cannot find data for tip, raise exception
                raise Exception(
                    'Could not find {} in character matrix'.format(label))
            else:
                data[tip_i, :, 0] = char_mtx.data[tip_lookup[label]]
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
        # TODO: Evaluate if this causes a memory problem
        full_mcp = np.zeros((internal_node_count, internal_node_count),
                            dtype=float)
        full_vcp = np.zeros(internal_node_count, dtype=float)
        # full_mcp = np.zeros((num_nodes, num_nodes), dtype=float)
        # full_vcp = np.zeros(num_nodes, dtype=float)

        for k in tree.postorder_edge_iter():
            i = k.head_node
            if len(i.child_nodes()) != 0:
                node_num_i = node_index_lookup[_get_node_label(i)] - tip_count
                for j in i.child_nodes():
                    tbl = 2./j.edge_length
                    full_mcp[node_num_i][node_num_i] += tbl
                    node_num_j = node_index_lookup[_get_node_label(j)]

                    if len(j.child_nodes()) == 0:
                        full_vcp[node_num_i] += (data[node_num_j,
                                                      x, 0] * tbl)
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
                # print(i.data['val'])
                # i.label = str(mle[nodenum[i]])
                # for j in i.child_nodes():
                #     node_num_j = _get_node_label(j)
                #     temp = (data[node_num_i, x, 0] - data[node_num_j, x, 0])
                #     sos += temp*temp / j.edge_length

    depth_headers = ['maximum_likelihood']
    if calc_std_err:
        depth_headers.append('standard_error')

    mtx_headers = {
        '0': row_headers,
        '1': tip_col_headers,
        '2': depth_headers
        }
    return tree, Matrix(data, headers=mtx_headers)


# .............................................................................
def calculate_continuous_ancestral_states_old(tree, sequences):
    """
    @summary: Calculates continuous ancestral states for a tree and a list of
                 sequences
    @note: sq_change which produces the same results as ML without SE
    @todo: Raise a different exception
    """
    match_tips_and_cont_values(tree, sequences)
    df = 0
    nodenum = {}
    count = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) == 0:
            i.data['val'] = float(i.data['cont_values'][0])
            i.data['valse'] = float(i.data['cont_values'][0])
        else:
            nodenum[i] = count
            count += 1
            df += 1
            i.data['val'] = 0.
            i.data['valse'] = 0.
    df -= 1
    # compute the mlest of the root
    fullMcp = np.zeros((df+1, df+1))
    fullVcp = np.zeros(df+1)
    count = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) != 0:
            nni = nodenum[i]
            for j in i.child_nodes():
                tbl = 2./j.edge_length
                fullMcp[nni][nni] += tbl
                if len(j.child_nodes()) == 0:
                    fullVcp[nni] += (j.data['val'] * tbl)
                else:
                    nnj = nodenum[j]
                    fullMcp[nni][nnj] -= tbl
                    fullMcp[nnj][nni] -= tbl
                    fullMcp[nnj][nnj] += tbl
            count += 1
    b = la.cho_factor(fullMcp)
    # these are the ML estimates for the ancestral states
    mle = la.cho_solve(b, fullVcp)
    sos = 0
    for k in tree.postorder_edge_iter():
        i = k.head_node
        if len(i.child_nodes()) != 0:
            i.data['val'] = mle[nodenum[i]]
            # print(i.data['val'])
            i.label = str(mle[nodenum[i]])
            for j in i.child_nodes():
                temp = (i.data['val'] - j.data['val'])
                sos += temp*temp / j.edge_length
    # print("Square Length: {}".format(sos))
    # calcSE
    """
    for i in tree.iternodes(order="postorder"):
        if i.istip == False:
            qpq = fullMcp[nodenum[i]][nodenum[i]]
            tm1 = np.delete(fullMcp,(nodenum[i]),axis=0)
            tm = np.delete(tm1,(nodenum[i]),axis=1)
            b = cho_factor(tm)
            sol = cho_solve(b,tm1[:,nodenum[i]])
            tempse = qpq - np.inner(tm1[:,nodenum[i]],sol)
            i.data['valse'] = math.sqrt(2*sos/(df*tempse))
    """
    return 0


# .............................................................................
def match_tips_and_cont_values(tree, seqs):
    """
    @summary: Match the tips in the tree with the values in the alignment
    @param tree: The tree to get the tips from
    @param seq: A list of Sequence objects to get the alignment values from
    @todo: Raise a different exception
    """
    for i in tree:
        i.data = {}
        if len(i.child_nodes()) == 0:
            test = False
            for j in seqs:
                if i.taxon.label == j.name:
                    test = True
                    i.data['cont_values'] = j.cont_values
                    break
            if not test:
                raise Exception('Could not find {} in cont_values'.format(
                    i.taxon.label))


# .............................................................................
# if __name__ == "__main__":
#     if len(sys.argv) != 3:
#         print("python "+sys.argv[0]+" newick.tre dataasphy")
#         sys.exit(0)
#     tree = dp.Tree.get(path=sys.argv[1], schema="newick")
#     seqs = aln_reader.read_phylip_cont_file(open(sys.argv[2], "r"))
#     match_tips_and_cont_values(tree,seqs)
#     sqch = calc_cont_anc_states(tree)
#     outfile = open("contanc_dp.tre","w")
#     outfile.write(tree.as_string(schema="newick"))
#     outfile.close()

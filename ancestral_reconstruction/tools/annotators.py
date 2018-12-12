"""Module containing tools for annotating trees with ancestral trait data
"""


# .............................................................................
def add_all_annotations(lm_tree, node_matrix, update=False):
    """Adds all annotations from the node matrix to the tree

    Args:
        lm_tree (TreeWrapper) : A Lifemapper tree object.
        node_matrix (Matrix) : A Lifemapper Matrix object with rows matching
            the nodes in the provided tree.
        update (bool) : Should any existing annotations be updated
    """
    annotations = {}
    for i, node_name in enumerate(node_matrix.get_row_headers()):
        anns = {}
        for j, col_name in enumerate(node_matrix.get_column_headers()):
            ann[col_name] = node_matrix.data[i, j]
        annotations[node_name] = anns
    lm_tree.annotate_tree(annotations, update=update)


# .............................................................................
def annotate_tree_with_label(lm_tree, node_matrix, label_column=0):
    """Annotates the tree by changing the label of each node

    Args:
        lm_tree (TreeWrapper) : A Lifemapper tree object.
        node_matrix (Matrix) : A Lifemapper Matrix object with rows matching
            the nodes in the provided tree.
        label_column (int) : The column in the matrix to use for node labels.
    """
    annotations = {}
    for i, node_name in enumerate(node_matrix.get_row_headers()):
        annotations[node_name] = node_matrix.data[i, label_column]
    lm_tree.annotate_tree(annotations, label_attribute=None)

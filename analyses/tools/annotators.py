"""Module containing tools for annotating trees with ancestral trait data.
"""


# .............................................................................
def add_all_annotations(lm_tree, node_matrix, update=False):
    """Adds all annotations from the node matrix to the tree.

    Args:
        lm_tree (TreeWrapper): A Lifemapper tree object.
        node_matrix (Matrix): A Lifemapper Matrix object with rows matching the
            nodes in the provided tree.
        update (:obj:`bool`, optional): Should any existing annotations be
            updated.  Defaults to False.
    """
    annotations = {}
    for i, node_name in enumerate(node_matrix.get_row_headers()):
        anns = {}
        for j, col_name in enumerate(node_matrix.get_column_headers()):
            anns[col_name] = node_matrix[i, j, 0]
        annotations[node_name] = anns
    lm_tree.annotate_tree(annotations, update=update)


# .............................................................................
def annotate_tree_with_label(lm_tree, node_matrix, label_column=0):
    """Annotates the tree by changing the label of each node.

    Args:
        lm_tree (TreeWrapper): A Lifemapper tree object.
        node_matrix (Matrix): A Lifemapper Matrix object with rows matching
            the nodes in the provided tree.
        label_column (:obj:`int`, optional): The column in the matrix to use
            for node labels.  Defaults to 0.
    """
    annotations = {}
    for i, node_name in enumerate(node_matrix.get_row_headers()):
        annotations[node_name] = node_matrix[i, label_column, 0]
    lm_tree.annotate_tree(annotations, label_attribute=None)

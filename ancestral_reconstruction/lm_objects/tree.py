"""
@summary: Module containing LmTree class
@author: CJ Grady
@version: 2.0
@status: alpha

@license: gpl2
@copyright: Copyright (C) 2018, University of Kansas Center for Research

          Lifemapper Project, lifemapper [at] ku [dot] edu,
          Biodiversity Institute,
          1345 Jayhawk Boulevard, Lawrence, Kansas, 66045, USA

          This program is free software; you can redistribute it and/or modify
          it under the terms of the GNU General Public License as published by
          the Free Software Foundation; either version 2 of the License, or (at
          your option) any later version.

          This program is distributed in the hope that it will be useful, but
          WITHOUT ANY WARRANTY; without even the implied warranty of
          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
          General Public License for more details.

          You should have received a copy of the GNU General Public License
          along with this program; if not, write to the Free Software
          Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
          02110-1301, USA.
@todo: Should we provide a method to collapse clades that only have one child?
@todo: Should pruning a tree collapse clades automatically?
@todo: Add method to remove annotations
@todo: Move label method out of internal functions
"""
import dendropy
import numpy as np

from ancestral_reconstruction.lm_objects.matrix import Matrix


# .............................................................................
# .                     Phylogenetic Tree Module Constants                    .
# .............................................................................
DEFAULT_TREE_SCHEMA = 'nexus'


# .............................
class PhyloTreeKeys(object):
    """
    @summary: Keys for phylogenetic trees
    """
    MTX_IDX = 'mx'  # The matrix index for this node
    SQUID = 'squid'  # This is the LM SQUID (species identifier) for the tip


# .............................................................................
class LmTreeException(Exception):
    """
    @summary: Wrapper around the base Exception class for tree related errors
    """
    pass


# .............................................................................
class TreeWrapper(dendropy.Tree):
    """
    @summary: Dendropy Tree wrapper that adds a little functionality and
                improves performance of some functions
    """
    # ..............................
    @classmethod
    def from_base_tree(cls, tree):
        """
        @summary: Creates a TreeWrapper object from a base dendropy.Tree
        """
        return cls.get(data=tree.as_string('nexus'), schema='nexus')

    # ..............................
    def add_node_labels(self, prefix=None, overwrite=False):
        """Add labels to the nodes in the tree

        Add labels to the unlabeled nodes in the tree.

        Args:
            prefix (:obj:`str`, optional): If provided, prefix the node labels
                with this string.
            overwrite (boolean): Indicates whether existing node labels should
                be overwritten or if they should be maintained.

        Note:
            This labels nodes the way that R does
        """
        self._label_tree_nodes(self.seed_node, len(self.get_labels()),
                               prefix=prefix, overwrite=overwrite)

    # ..............................
    def annotate_tree(self, attribute_name, annotation_pairs,
                      label_attribute='label', update=False):
        """
        @summary: Annotates the nodes of the tree
        @param attribute_name: The name of the annotation attribute to add
        @param annotation_pairs: A dictionary of label keys with annotation
                                    values
        @param label_attribute: If this is provided, use this annotation
                                    attribute as the key instead of the label
        """
        label_method = self._get_label_method(label_attribute)

        for taxon in self.taxon_namespace:
            try:
                # label = getattr(taxon, labelAttribute)
                label = label_method(taxon)

                if taxon.annotations.get_value(attribute_name) is not None:
                    if update:
                        # Remove existing values
                        for ann in taxon.annotations.findall(
                                                        name=attribute_name):
                            taxon.annotations.remove(ann)
                        # Set new value
                        taxon.annotations.add_new(attribute_name,
                                                  annotation_pairs[label])
                else:
                    taxon.annotations.add_new(attribute_name,
                                              annotation_pairs[label])
            except KeyError:
                # Pass if a label is not found in the dictionary,
                #     otherwise fail
                pass

    # ..............................
    def get_annotations(self, annotation_attribute):
        """
        @summary: Gets a list of (label, annotation) pairs
        @param annotation_attribute: The annotation attribute to retrieve
        """
        annotations = []
        for taxon in self.taxon_namespace:
            att = taxon.annotations.get_value(annotation_attribute)
            annotations.append((taxon.label, att))
        return annotations

    # ..............................
    def get_distance_matrix(self, label_attribute='label',
                            ordered_labels=None):
        """
        @summary: Get a Matrix object of phylogenetic distances between tips
                    using a lower memory footprint
        @param label_attribute: The attribute of the tips to use as labels for
                                    the matrix
        @param ordered_labels: If provided, use this order of labels
        """
        label_method = self._get_label_method(label_attribute)

        # Get list of labels
        if ordered_labels is None:
            ordered_labels = []
            for taxon in self.taxon_namespace:
                ordered_labels.append(label_method(taxon))

        label_lookup = dict([(
                ordered_labels[i], i) for i in range(len(ordered_labels))])

        dist_mtx = np.zeros((len(ordered_labels), len(ordered_labels)),
                            dtype=float)

        # Build path lookup dictionary
        path_lookups = {}
        for taxon in self.taxon_namespace:
            n = self.find_node_for_taxon(taxon)
            l = []
            while n is not None:
                if n.edge_length is not None:
                    # If this is still too big, drop id and do some logic with
                    #    lengths
                    l.append((id(n), n.edge_length))
                n = n.parent_node
            path_lookups[taxon.label] = l

        num_taxa = len(self.taxon_namespace)
        for i_1 in xrange(num_taxa - 1):
            taxon1 = self.taxon_namespace[i_1]

            label = label_method(taxon1)
            # Check for matrix index
            try:
                idx1 = label_lookup[label]

                # Build path to root for taxon 1
                # path_labels = []
                o_dist = 0.0
                t_path = path_lookups[taxon1.label]
                t_labels = []
                for label, p_dist in t_path:
                    o_dist += p_dist
                    t_labels.append(label)

                for i_2 in xrange(i_1, num_taxa):
                    taxon2 = self.taxon_namespace[i_2]

                    try:
                        idx2 = label_lookup[label_method(taxon2)]

                        # Initialize distance for these two taxa
                        dist = o_dist

                        # Loop through path back to root
                        t2_path = path_lookups[taxon2.label]

                        for label, p_dist in t2_path:
                            if label in t_labels:
                                dist -= p_dist
                            else:
                                dist += p_dist

                        # mrca = pdm.mrca(taxon1, taxon2)
                        # dist = pdm.patristic_distance(taxon1, taxon2)
                        dist_mtx[idx1, idx2] = dist
                        dist_mtx[idx2, idx1] = dist
                    except:
                        pass
            except:
                pass

        distance_matrix = Matrix(dist_mtx, headers={'0': ordered_labels,
                                                    '1': ordered_labels})
        return distance_matrix

    # ..............................
    def get_distance_matrix_dendropy(self, label_attribute='label',
                                     ordered_labels=None):
        """
        @summary: Get a Matrix object of phylogenetic distances between tips
        @param label_attribute: The attribute of the tips to use as labels for
                                    the matrix
        @param ordered_labels: If provided, use this order of labels
        """
        label_method = self._get_label_method(label_attribute)

        # Get list of labels
        if ordered_labels is None:
            ordered_labels = []
            for taxon in self.taxon_namespace:
                ordered_labels.append(label_method(taxon))

        label_lookup = dict([(
                ordered_labels[i], i) for i in range(len(ordered_labels))])

        dist_mtx = np.zeros((len(ordered_labels), len(ordered_labels)),
                            dtype=float)

        pdm = self.phylogenetic_distance_matrix()

        for taxon1 in self.taxon_namespace:
            label = label_method(taxon1)
            # Check for matrix index
            try:
                idx1 = label_lookup[label]

                for taxon2 in self.taxon_namespace:
                    try:
                        idx2 = label_lookup[label_method(taxon2)]
                        # mrca = pdm.mrca(taxon1, taxon2)
                        dist = pdm.patristic_distance(taxon1, taxon2)
                        dist_mtx[idx1, idx2] = dist
                    except:
                        pass
            except:
                pass

        distance_matrix = Matrix(dist_mtx, headers={'0': ordered_labels,
                                                    '1': ordered_labels})
        return distance_matrix

    # ..............................
    def get_labels(self):
        """
        @summary: Get tip labels for a clade
        @note: Bottom-up order
        """
        labels = []
        for taxon in self.taxon_namespace:
            labels.append(taxon.label)

        labels.reverse()
        return labels

    # ..............................
    def get_variance_covariance_matrix(self, label_attribute='label',
                                       ordered_labels=None):
        """
        @summary: Get a Matrix object of variance / covariance for tips in tree
        @param label_attribute: The attribute of the tips to use as labels for
                                    the matrix
        @param ordered_labels: If provided, use this order of labels
        @todo: Consider using the node labeling function to get the indices.
                    Haven't changed yet because this is producing the correct
                    result
        """
        if not self.has_branch_lengths():
            raise Exception('Cannot create VCV without branch lengths')

        label_method = self._get_label_method(label_attribute)

        # Get list of labels
        if ordered_labels is None:
            ordered_labels = []
            for taxon in self.taxon_namespace:
                ordered_labels.append(label_method(taxon))

        label_lookup = dict([(
                ordered_labels[i], i) for i in range(len(ordered_labels))])

        n = len(orderedLabels)
        vcv = np.zeros((n, n), dtype=float)

        edges = []
        for edge in self.postorder_edge_iter():
            edges.append(edge)

        edges.reverse()

        for edge in edges:
            try:
                el = edge.head_node.distance_from_root()
            except:
                el = None
            if el is not None:
                child_nodes = edge.head_node.child_nodes()
                if len(child_nodes) == 0:
                    idx = label_lookup[label_method(edge.head_node.taxon)]
                    vcv[idx, idx] = edge.head_node.distance_from_root()
                else:
                    left_child, right_child = edge.head_node.child_nodes()
                    left_tips = [
                        label_lookup[
                            label_method(tipNode.taxon)
                                    ] for tip_node in left_child.leaf_nodes()]
                    right_tips = [
                        label_lookup[
                            label_method(tipNode.taxon)
                                    ] for tip_node in right_child.leaf_nodes()]
                    # if len(leftTips) > 1 and len(rightTips) > 1:
                    for l in left_tips:
                        for r in right_tips:
                            vcv[l, r] = vcv[r, l] = el

            for node in self.leaf_nodes():
                idx = label_ookup[label_method(node.taxon)]
                vcv[idx, idx] = node.distance_from_root()

        vcv_matrix = Matrix(vcv, headers={'0': ordered_labels,
                                          '1': ordered_labels})
        return vcv_matrix

    # ..............................
    def has_branch_lengths(self):
        """
        @summary: Returns boolean indicating if the tree has branch lengths for
                      every clade
        """
        try:
            self.minmax_leaf_distance_from_root()
            return True
        except:
            return False

    # ..............................
    def has_polytomies(self):
        """
        @summary: Returns boolean indicating if the tree has polytomies
        """
        for n in self.nodes():
            if len(n.child_nodes()) > 2:
                return True
        return False

    # ..............................
    def is_binary(self):
        """
        @summary: Returns a boolean indicating if the tree is binary
        @note: Checks that every clade has either zero or two children
        """
        for n in self.nodes():
            if not len(n.child_nodes()) in [0, 2]:
                return False
        return True

    # ..............................
    def is_ultrametric(self, rel_tol=1e-03):
        """
        @summary: Check if the tree is ultrametric
        @param rel_tol: The relative tolerance to determine if the min and max
                        are equal.  We will say they are equal if they are
                        99.9%.
        @note: To be ultrametric, the branch length from root to tip must be
                   equal for all tips
        """
        try:
            min_bl, max_bl = self.minmax_leaf_distance_from_root()
            return np.isclose(min_bl, max_bl, rtol=rel_tol)
        except:
            pass
        return False

    # ..............................
    def prune_tips_without_attribute(self,
                                     search_attribute=PhyloTreeKeys.MTX_IDX):
        """
        @summary: Prunes the tree of any tips that don't have the specified
                      attribute
        """
        prune_taxa = []
        for taxon in self.taxon_namespace:
            val = None
            try:
                val = taxon.annotations.get_value(search_attribute)
            except:
                pass
            if val is None:
                prune_taxa.append(taxon)

        self.prune_taxa(prune_taxa)
        self.purge_taxon_namespace()

    # ..............................
    def _annotation_method(self, label_attribute):
        """
        @summary: Use the label attribute as the node label
        """
        def label_method(taxon):
            return taxon.annotations.get_value(label_attribute)
        return label_method

    # ..............................
    def _get_label_method(self, label_attribute):
        """
        @summary: Get the function to be used for retrieving labels
        """
        if label_attribute.lower() == 'label':
            return self._taxon_label_method
        else:
            return self._annotation_method(label_attribute)

    # ..............................
    def _label_tree_nodes(self, node, i, prefix=None, overwrite=False):
        """
        @summary: Private function to do the work when labeling nodes
        @note: Recursive
        """
        cn = node.child_nodes()

        # If we have children, label and recurse
        if len(cn) > 0:
            if node.label is None or overwrite:
                if prefix is not None:
                    node.label = '{}{}'.format(prefix, i)
                else:
                    node.label = str(i)
                i += 1
            # Loop through children and label nodes
            for child in cn:
                i = self._label_tree_nodes(child, i, prefix=prefix,
                                           overwrite=overwrite)
        # Return the current i value
        return i

    # ..............................
    def _taxon_label_method(self, taxon):
        """
        @summary: Use the taxon label for the label
        """
        return taxon.label

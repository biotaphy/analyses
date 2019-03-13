"""Module containing code for calculating phylogenetic beta diversity

Note:
    * Numpy and Scipy may have the methods needed for these computations.  Feel
        free to take advantage of those to libraries as they are already
        dependencies.
"""
import numpy as np

from analyses.lm_objects.matrix import Matrix


# .............................................................................
def get_species_index_lookup(pam):
    """Creates a lookup dictionary for species in a matrix

    Args:
        pam (:obj:`Matrix`): A Lifemapper Matrix object with presence absence
            values.

    Returns:
        A dictionary of species keys with matrix index values.
    """
    return dict(
        [(sp, idx) for (idx, sp) in enumerate(pam.get_column_headers())])


# .............................................................................
def calculate_phylo_beta_diversity_jaccard(pam, tree):
    """Calculates phylogenetic beta diversity for the jaccard index family.

    Args:
        pam (:obj:`Matrix`): A Lifemapper Matrix object with presence absence
            values (site rows by species columns).
        tree (:obj:`TreeWrapper`): A TreeWrapper object for a wrapped Dendropy
            phylogenetic tree.

    Returns:
        Phylogenetic beta diversity matrics (species by species)
            * beta_jtu: ADD DESCRIPTION
            * phylo_beta_jtu: ADD DESCRIPTION
            * beta_jne: ADD DESCRIPTION
            * phylo_beta_jne: ADD DESCRIPTION
            * beta_jac: ADD DESCRIPTION
            * phylo_beta_jac: ADD DESCRIPTION

    Note:
        * It looks like the scipy.spatial.distance.jaccard method may be useful
            here.

    Todo:
        * Fill in method documentation
        * Fill in method
    """
    # Get a lookup dictionary for the matrix index of each species in the PAM
    #    in case they are not in the same order as the taxa in the tree
    species_lookup = get_species_index_lookup(pam)

    # Build a header dictionary, all of the returned matricies will have the
    #    same headers, site rows by site columns.
    # Note: This will differ from the R method because each site will be
    #    present in both the rows and the columns.
    mtx_headers = {
        '0': pam.get_row_headers(),  # Row headers
        '1': pam.get_row_headers()  # Column headers
    }

    num_sites = pam.data.shape[0]  # Get the number of sites in the PAM

    # Note: For ease of development, use these numpy arrays for the
    #    computations.  They will be wrapped into a Matrix object when they are
    #    returned from the function.
    beta_jtu_data = np.zeros((num_sites, num_sites), dtype=np.float)
    phylo_beta_jtu_data = np.zeros((num_sites, num_sites), dtype=np.float)
    beta_jne_data = np.zeros((num_sites, num_sites), dtype=np.float)
    phylo_beta_jne_data = np.zeros((num_sites, num_sites), dtype=np.float)
    beta_jac_data = np.zeros((num_sites, num_sites), dtype=np.float)
    phylo_beta_jac_data = np.zeros((num_sites, num_sites), dtype=np.float)

    # TODO: Compute phylo beta diversity for jaccard index family

    return (
        Matrix(beta_jtu_data, headers=mtx_headers),
        Matrix(phylo_beta_jtu_data, headers=mtx_headers),
        Matrix(beta_jne_data, headers=mtx_headers),
        Matrix(phylo_beta_jne_data, headers=mtx_headers),
        Matrix(beta_jac_data, headers=mtx_headers),
        Matrix(phylo_beta_jac_data, headers=mtx_headers))


# .............................................................................
def calculate_phylo_beta_diversity_sorensen(pam, tree):
    """Calculates phylogenetic beta diversity for the sorensen index family.

    Args:
        pam (:obj:`Matrix`): A Lifemapper Matrix object with presence absence
            values.
        tree (:obj:`TreeWrapper`): A TreeWrapper object for a wrapped Dendropy
            phylogenetic tree.

    Returns:
        Phylogenetic beta diversity matrics (species by species)
            * beta_sim: ADD DESCRIPTION
            * phylo_beta_sim: ADD DESCRIPTION
            * beta_sne: ADD DESCRIPTION
            * phylo_beta_sne: ADD DESCRIPTION
            * beta_sor: ADD DESCRIPTION
            * phylo_beta_sor: ADD DESCRIPTION

    Todo:
        * Fill in method documentation
        * Fill in method
    """
    # Get a lookup dictionary for the matrix index of each species in the PAM
    #    in case they are not in the same order as the taxa in the tree
    species_lookup = get_species_index_lookup(pam)

    # Build a header dictionary, all of the returned matricies will have the
    #    same headers, site rows by site columns.
    # Note: This will differ from the R method because each site will be
    #    present in both the rows and the columns.
    mtx_headers = {
        '0': pam.get_row_headers(),  # Row headers
        '1': pam.get_row_headers()  # Column headers
    }

    num_sites = pam.data.shape[0]  # Get the number of sites in the PAM

    # Note: For ease of development, use these numpy arrays for the
    #    computations.  They will be wrapped into a Matrix object when they are
    #    returned from the function.
    beta_sim_data = np.zeros((num_sites, num_sites), dtype=np.float)
    phylo_beta_sim_data = np.zeros((num_sites, num_sites), dtype=np.float)
    beta_sne_data = np.zeros((num_sites, num_sites), dtype=np.float)
    phylo_beta_sne_data = np.zeros((num_sites, num_sites), dtype=np.float)
    beta_sor_data = np.zeros((num_sites, num_sites), dtype=np.float)
    phylo_beta_sor_data = np.zeros((num_sites, num_sites), dtype=np.float)

    # TODO: Compute phylo beta diversity for sorensen index family

    return (
        Matrix(beta_sim_data, headers=mtx_headers),
        Matrix(phylo_beta_sim_data, headers=mtx_headers),
        Matrix(beta_sne_data, headers=mtx_headers),
        Matrix(phylo_beta_sne_data, headers=mtx_headers),
        Matrix(beta_sor_data, headers=mtx_headers),
        Matrix(phylo_beta_sor_data, headers=mtx_headers))

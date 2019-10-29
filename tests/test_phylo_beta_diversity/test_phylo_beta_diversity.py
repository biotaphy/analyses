"""This module is used for testing phylogenetic beta diversity

Notes:
    * Uses pytest style testing.
"""
import numpy as np

from lmpy import Matrix, TreeWrapper

from analyses.phylo_beta_diversity import phylo_beta_diversity as pbd


# .............................................................................
class Test_phylo_beta_diversity_jaccard(object):
    """Test phylogenetic beta diversity using the jaccard index family
    """
    # .....................................
    def test_valid(self, valid_phylo_beta_diversity_package):
        """Test the method with valid data

        Note:
            * Test values were determined from example at
                https://rdrr.io/rforge/betapart/man/phylo.beta.pair.html
        """
        (pam_fn, tree_fn, test_beta_jac_fn, test_beta_jne_fn, test_beta_jtu_fn,
         _, _, _, test_phylo_beta_jac_fn, test_phylo_beta_jne_fn,
         test_phylo_beta_jtu_fn, _, _, _) = valid_phylo_beta_diversity_package

        with open(pam_fn) as in_f:
            pam = Matrix.load_csv(in_f, num_header_rows=1, num_header_cols=1)
        tree = TreeWrapper.from_filename(tree_fn)
        with open(test_beta_jac_fn) as in_f:
            test_beta_jac = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)
        with open(test_beta_jne_fn) as in_f:
            test_beta_jne = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)
        with open(test_beta_jtu_fn) as in_f:
            test_beta_jtu = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)
        with open(test_phylo_beta_jac_fn) as in_f:
            test_phylo_beta_jac = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)
        with open(test_phylo_beta_jne_fn) as in_f:
            test_phylo_beta_jne = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)
        with open(test_phylo_beta_jtu_fn) as in_f:
            test_phylo_beta_jtu = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)

        (beta_jtu, phylo_beta_jtu, beta_jne, phylo_beta_jne, beta_jac,
         phylo_beta_jac) = pbd.calculate_phylo_beta_diversity_jaccard(
             pam, tree)
        # Check matrix outputs to see if they are within tolerance
        assert np.allclose(beta_jtu, test_beta_jtu)
        assert np.allclose(phylo_beta_jtu, test_phylo_beta_jtu)
        assert np.allclose(beta_jne, test_beta_jne)
        assert np.allclose(phylo_beta_jne, test_phylo_beta_jne)
        assert np.allclose(beta_jac, test_beta_jac)
        assert np.allclose(phylo_beta_jac, test_phylo_beta_jac)

    # .....................................
    def test_extra_species_in_pam(self):
        pass

    # .....................................
    def test_extra_species_in_tree(self):
        pass


# .............................................................................
class Test_phylo_beta_diversity_sorensen(object):
    """Test phylogenetic beta diversity using the sorensen index family
    """
    # .....................................
    def test_valid(self, valid_phylo_beta_diversity_package):
        """Test the method with valid data

        Note:
            * Test values were determined from example at
                https://rdrr.io/rforge/betapart/man/phylo.beta.pair.html
        """
        (pam_fn, tree_fn, _, _, _, test_beta_sim_fn, test_beta_sne_fn,
         test_beta_sor_fn, _, _, _, test_phylo_beta_sim_fn,
         test_phylo_beta_sne_fn, test_phylo_beta_sor_fn
         ) = valid_phylo_beta_diversity_package

        with open(pam_fn) as in_f:
            pam = Matrix.load_csv(in_f, num_header_rows=1, num_header_cols=1)
        tree = TreeWrapper.from_filename(tree_fn)
        with open(test_beta_sim_fn) as in_f:
            test_beta_sim = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)
        with open(test_beta_sne_fn) as in_f:
            test_beta_sne = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)
        with open(test_beta_sor_fn) as in_f:
            test_beta_sor = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)
        with open(test_phylo_beta_sim_fn) as in_f:
            test_phylo_beta_sim = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)
        with open(test_phylo_beta_sne_fn) as in_f:
            test_phylo_beta_sne = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)
        with open(test_phylo_beta_sor_fn) as in_f:
            test_phylo_beta_sor = Matrix.load_csv(
                in_f, num_header_rows=1, num_header_cols=1)

        (beta_sim, phylo_beta_sim, beta_sne, phylo_beta_sne, beta_sor,
         phylo_beta_sor) = pbd.calculate_phylo_beta_diversity_sorensen(
             pam, tree)
        # Check matrix outputs to see if they are within tolerance
        assert np.allclose(beta_sim, test_beta_sim)
        assert np.allclose(phylo_beta_sim, test_phylo_beta_sim)
        assert np.allclose(beta_sne, test_beta_sne)
        assert np.allclose(phylo_beta_sne, test_phylo_beta_sne)
        assert np.allclose(beta_sor, test_beta_sor)
        assert np.allclose(phylo_beta_sor, test_phylo_beta_sor)

    # .....................................
    def test_extra_species_in_pam(self):
        pass

    # .....................................
    def test_extra_species_in_tree(self):
        pass

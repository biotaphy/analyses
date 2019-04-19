"""This module is used for testing phylogenetic beta diversity

Notes:
    * Uses pytest style testing.
"""
import numpy as np

from analyses.lm_objects.matrix import Matrix
from analyses.lm_objects.tree import TreeWrapper
from analyses.phylo_beta_diversity import phylo_beta_diversity as pbd


# .............................................................................
class Test_phylo_beta_diversity_jaccard(object):
    """Test phylogenetic beta diversity using the jaccard index family
    """
    # .....................................
    def test_valid(self):
        """Test the method with valid data

        Note:
            * Test values were determined from example at
                https://rdrr.io/rforge/betapart/man/phylo.beta.pair.html
        """
        site_headers = {
            '0': ['A', 'B', 'C', 'D', 'E', 'F'],
            '1': ['A', 'B', 'C', 'D', 'E', 'F']}
        pam = Matrix(
            np.array([
                [1, 1, 1, 0, 0, 0],
                [0, 1, 1, 1, 0, 0],
                [0, 0, 1, 1, 1, 0],
                [0, 0, 1, 1, 1, 1],
                [0, 0, 0, 1, 1, 1],
                [1, 0, 0, 1, 1, 1]]),
            headers={
                '0': ['A', 'B', 'C', 'D', 'E', 'F'],
                '1': ['sp1', 'sp2', 'sp3', 'sp4', 'sp5', 'sp6']})
        tree = TreeWrapper.get(
            data='(((sp1:1,sp2:1):5,(sp3:3,sp4:3):3):2,(sp5:7,sp6:7):1);',
            schema='newick')

        test_beta_jtu = Matrix(
            np.array([
                [1.0, 0.5, 0.8, 0.8, 1.0, 0.8],
                [0.5, 1.0, 0.5, 0.5, 0.8, 0.8],
                [0.8, 0.5, 1.0, 0.0, 0.5, 0.5],
                [0.8, 0.5, 0.0, 1.0, 0.0, 0.4],
                [1.0, 0.8, 0.5, 0.0, 1.0, 0.0],
                [0.8, 0.8, 0.5, 0.4, 0.0, 1.0]]), headers=site_headers)

        test_beta_jne = Matrix(
            np.array([
                [1.0, 0.0, 0.0, 0.0333333, 0.0, 0.0333333],
                [0.0, 1.0, 0.0, 0.1, 0.0, 0.0333333],
                [0.0, 0.0, 1.0, 0.25, 0.0, 0.1],
                [0.0333333, 0.1, 0.25, 1.0, 0.25, 0.0],
                [0.0, 0.0, 0.0, 0.25, 1.0, 0.25],
                [0.0333333, 0.0333333, 0.1, 0.0, 0.25, 1.0]]),
            headers=site_headers)

        test_beta_jac = Matrix(
            np.array([
                [1.0, 0.5, 0.8, 0.8333333, 1.0, 0.8333333],
                [0.5, 1.0, 0.5, 0.6, 0.8, 0.8333333],
                [0.8, 0.5, 1.0, 0.25, 0.5, 0.6],
                [0.8333333, 0.6, 0.25, 1.0, 0.25, 0.4],
                [1.0, 0.8, 0.5, 0.25, 1.0, 0.25],
                [0.8333333, 0.8333333, 0.6, 0.4, 0.25, 1.0]]),
            headers=site_headers)

        test_phylo_beta_jtu = Matrix(
            np.array([
                [1.0, 0.125, 0.6363636, 0.6363636, 0.8, 0.4210526],
                [0.125, 1.0, 0.5217391, 0.5217391, 0.6923077, 0.3809524],
                [0.6363636, 0.5217391, 1.0, 0.0, 0.2727273, 0.2727273],
                [0.6363636, 0.5217391, 0.0, 1.0, 0.0, 0.2068966],
                [0.8, 0.6923077, 0.2727273, 0.0, 1.0, 0.0],
                [0.4210526, 0.3809524, 0.2727273, 0.2068966, 0.0, 1.0]]),
            headers=site_headers)

        test_phylo_beta_jne = Matrix(
            np.array([
                [1.0, 0.09722222, 0.05594406, 0.12121212, 0.04848485,
                    0.24561404],
                [0.09722222, 1.0, 0.03826087, 0.13451087, 0.05769231,
                    0.22510823],
                [0.05594406, 0.03826087, 1.0, 0.26923077, 0.11188811,
                    0.22727273],
                [0.12121212, 0.13451087, 0.26923077, 1.0, 0.11538462,
                    0.07435345],
                [0.04848485, 0.05769231, 0.11188811, 0.11538462, 1.0,
                    0.20689655],
                [0.24561404, 0.22510823, 0.22727273, 0.07435345, 0.20689655,
                    1.0]]),
            headers=site_headers)  # 0.3826087 [b,c]

        test_phylo_beta_jac = Matrix(
            np.array([
                [1.0, 0.2222222, 0.6923077, 0.7575758, 0.8484848, 0.6666667],
                [0.2222222, 1.0, 0.56, 0.65625, 0.75, 0.6060606],
                [0.6923077, 0.56, 1.0, 0.2692308, 0.3846154, 0.5],
                [0.7575758, 0.65625, 0.2692308, 1.0, 0.1153846, 0.28125],
                [0.8484848, 0.75, 0.3846154, 0.1153846, 1.0, 0.2068966],
                [0.6666667, 0.6060606, 0.5, 0.28125, 0.2068966, 1.0]]),
            headers=site_headers)

        (beta_jtu, phylo_beta_jtu, beta_jne, phylo_beta_jne, beta_jac,
         phylo_beta_jac) = pbd.calculate_phylo_beta_diversity_jaccard(
             pam, tree)
        # Check matrix outputs to see if they are within tolerance
        assert np.all(np.isclose(beta_jtu.data, test_beta_jtu.data))
        assert np.all(
            np.isclose(phylo_beta_jtu.data, test_phylo_beta_jtu.data))
        assert np.all(np.isclose(beta_jne.data, test_beta_jne.data))
        assert np.all(
            np.isclose(phylo_beta_jne.data, test_phylo_beta_jne.data))
        assert np.all(np.isclose(beta_jac.data, test_beta_jac.data))
        assert np.all(
            np.isclose(phylo_beta_jac.data, test_phylo_beta_jac.data))

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
    def test_valid(self):
        """Test the method with valid data

        Note:
            * Test values were determined from example at
                https://rdrr.io/rforge/betapart/man/phylo.beta.pair.html
        """
        site_headers = {
            '0': ['A', 'B', 'C', 'D', 'E', 'F'],
            '1': ['A', 'B', 'C', 'D', 'E', 'F']}
        pam = Matrix(
            np.array([
                [1, 1, 1, 0, 0, 0],
                [0, 1, 1, 1, 0, 0],
                [0, 0, 1, 1, 1, 0],
                [0, 0, 1, 1, 1, 1],
                [0, 0, 0, 1, 1, 1],
                [1, 0, 0, 1, 1, 1]]),
            headers={
                '0': ['A', 'B', 'C', 'D', 'E', 'F'],
                '1': ['sp1', 'sp2', 'sp3', 'sp4', 'sp5', 'sp6']})
        tree = TreeWrapper.get(
            data='(((sp1:1,sp2:1):5,(sp3:3,sp4:3):3):2,(sp5:7,sp6:7):1);',
            schema='newick')

        test_beta_sim = Matrix(
            np.array([
                [1.0, 0.3333333, 0.6666667, 0.6666667, 1.0, 0.6666667],
                [0.3333333, 1.0, 0.3333333, 0.3333333, 0.6666667, 0.6666667],
                [0.6666667, 0.3333333, 1.0, 0.0, 0.3333333, 0.3333333],
                [0.6666667, 0.3333333, 0.0, 1.0, 0.0, 0.25],
                [1.0, 0.6666667, 0.3333333, 0.0, 1.0, 0.0],
                [0.6666667, 0.6666667, 0.3333333, 0.25, 0.0, 1.0]]),
            headers=site_headers)

        test_beta_sne = Matrix(
            np.array([
                [1.0, 0.0, 0.0, 0.04761905, 0.0, 0.04761905],
                [0.0, 1.0, 0.0, 0.0952381, 0.0, 0.04761905],
                [0.0, 0.0, 1.0, 0.14285714, 0.0, 0.0952381],
                [0.04761905, 0.0952381, 0.14285714, 1.0, 0.14285714, 0.0],
                [0.0, 0.0, 0.0, 0.14285714, 1.0, 0.14285714],
                [0.04761905, 0.04761905, 0.09523810, 0.0, 0.14285714, 1.0]]),
            headers=site_headers)

        test_beta_sor = Matrix(
            np.array([
                [1.0, 0.3333333, 0.6666667, 0.7142857, 1.0, 0.7142857],
                [0.3333333, 1.0, 0.3333333, 0.4285714, 0.6666667, 0.7142857],
                [0.6666667, 0.3333333, 1.0, 0.1428571, 0.3333333, 0.4285714],
                [0.7142857, 0.4285714, 0.1428571, 1.0, 0.1428571, 0.25],
                [1.0, 0.6666667, 0.3333333, 0.1428571, 1.0, 0.1428571],
                [0.7142857, 0.7142857, 0.4285714, 0.25, 0.1428571, 1.0]]),
            headers=site_headers)  # 0.1428571 -> 0.4285714

        test_phylo_beta_sim = Matrix(
            np.array([
                [1.0, 0.06666667, 0.46666667, 0.46666667, 0.66666667,
                    0.26666667],
                [0.06666667, 1.0, 0.35294118, 0.35294118, 0.52941176,
                    0.23529412],
                [0.46666667, 0.35294118, 1.0, 0.0, 0.157894, 0.15789474],
                [0.46666667, 0.35294118, 0.0, 1.0, 0.0, 0.11538462],
                [0.66666667, 0.52941176, 0.15789474, 0.0, 1.0, 0.0],
                [0.26666667, 0.23529412, 0.15789474, 0.11538462, 0.0, 1.0]]),
            headers=site_headers)  # phylo

        test_phylo_beta_sne = Matrix(
            np.array([
                [1.0, 0.05833333, 0.0627451, 0.14308943, 0.07017544,
                    0.23333333],
                [0.05833333, 1.0, 0.03594771, 0.13543092, 0.07058824,
                    0.19948849],
                [0.0627451, 0.03594771, 1.0, 0.15555556, 0.0802005, 0.1754386],
                [0.14308943, 0.13543092, 0.15555556, 1.0, 0.06122449,
                    0.04825175],
                [0.07017544, 0.07058824, 0.08020050, 0.06122449, 1.0,
                    0.11538462],
                [0.23333333, 0.19948849, 0.1754386, 0.04825175, 0.11538462,
                    1.0]]),
            headers=site_headers)  # phylo; 0.5833333 -> 0.05833333

        test_phylo_beta_sor = Matrix(
            np.array([
                [1.0, 0.125, 0.52941176, 0.60975610, 0.73684211, 0.5],
                [0.125, 1.0, 0.38888889, 0.48837209, 0.6, 0.43478261],
                [0.52941176, 0.388888889, 1.0, 0.15555556, 0.23809524,
                    0.3333333],
                [0.6097561, 0.48837209, 0.15555556, 1.0, 0.06122449,
                    0.16363636],
                [0.73684211, 0.6, 0.23809524, 0.06122449, 1.0, 0.11538462],
                [0.5, 0.43478261, 0.33333333, 0.16363636, 0.11538462, 1.0]]),
            headers=site_headers)  # phylo

        (beta_sim, phylo_beta_sim, beta_sne, phylo_beta_sne, beta_sor,
         phylo_beta_sor) = pbd.calculate_phylo_beta_diversity_sorensen(
             pam, tree)
        # Check matrix outputs to see if they are within tolerance
        assert np.all(np.isclose(beta_sim.data, test_beta_sim.data))
        assert np.all(
            np.isclose(phylo_beta_sim.data, test_phylo_beta_sim.data))
        assert np.all(np.isclose(beta_sne.data, test_beta_sne.data))
        assert np.all(
            np.isclose(phylo_beta_sne.data, test_phylo_beta_sne.data))
        assert np.all(np.isclose(beta_sor.data, test_beta_sor.data))
        assert np.all(
            np.isclose(phylo_beta_sor.data, test_phylo_beta_sor.data))

    # .....................................
    def test_extra_species_in_pam(self):
        pass

    # .....................................
    def test_extra_species_in_tree(self):
        pass

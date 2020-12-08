import os
import unittest
import pytest

from copy import deepcopy

from .testing_utils import get_test_data


class PhyloClusterTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.refpkg import ReferencePackage
        from treesapp.rel_evo_dist import RedTree
        from ete3 import Tree
        self.refpkg = ReferencePackage()
        self.refpkg.f__json = get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.refpkg.slurp()
        self.taxa_tree = self.refpkg.taxonomically_label_tree()

        self.red_tree = RedTree()
        self.red_tree.decorate_rel_dist(self.taxa_tree)

        self.mock_tree = Tree("(A:1,(B:1,(E:1,D:1):0.5):0.5);")
        x = 0
        for n in self.mock_tree.traverse("postorder"):  # type: Tree
            if not n.name:
                n.name = str(x)
                x += 1
        return

    def test_partition_nodes(self):
        from treesapp.phylo_cluster import PhyloClust
        p_clust = PhyloClust()
        p_clust.alpha = 2
        node_paritions = p_clust.partition_nodes(tree=deepcopy(self.mock_tree))
        self.assertEqual(4, len(node_paritions))
        return

    def test_define_tree_clusters(self):
        from treesapp.phylo_cluster import PhyloClust
        p_clust = PhyloClust()
        p_clust.alpha = 0.1
        cluster_map = p_clust.define_tree_clusters(self.taxa_tree)
        return

    def test_group_rel_dist(self):
        from treesapp.phylo_cluster import PhyloClust
        p_clust = PhyloClust()

        # Test grouping by species-level
        rel_dists = p_clust.group_rel_dists(tree=self.taxa_tree, hierarchy=self.refpkg.taxa_trie)
        for group, dists in rel_dists.items():
            self.assertEqual('s', group[0])
            self.assertTrue(1 > max(dists))

        # Test grouping by a deeper taxonomic relationship, family
        rel_dists = p_clust.group_rel_dists(self.taxa_tree, self.refpkg.taxa_trie, group_rank="family", norm=False)
        return

    def test_cluster_phylogeny(self):
        from treesapp.phylo_cluster import cluster_phylogeny
        with pytest.raises(SystemExit):
            cluster_phylogeny([])
        return


if __name__ == '__main__':
    unittest.main()

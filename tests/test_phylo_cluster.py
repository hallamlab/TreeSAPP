import os
import unittest
import pytest

from copy import deepcopy
from ete3 import Tree, TreeNode

from .testing_utils import get_test_data


class PhyloClusterTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.refpkg import ReferencePackage
        from treesapp.rel_evo_dist import RedTree
        self.refpkg = ReferencePackage()
        self.refpkg.f__json = get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.refpkg.slurp()
        self.taxa_tree = self.refpkg.taxonomically_label_tree()

        self.red_tree = RedTree()
        self.red_tree.decorate_rel_dist(self.taxa_tree)

        self.mock_tree = Tree("(A:1,(B:0.1,(E:0.08,D:0.02):0.2):0.2);")
        x = 0
        for n in self.mock_tree.traverse("postorder"):  # type: Tree
            if not n.name:
                n.name = str(x)
                x += 1
        return

    def test_partition_nodes(self):
        from treesapp.phylo_cluster import PhyloClust
        p_clust = PhyloClust()
        p_clust.alpha = 0.1
        node_paritions = p_clust.partition_nodes(tree=deepcopy(self.mock_tree))
        # Ensure all the tree's nodes are present in the partitions
        count = 0
        for st in node_paritions.values():
            for _ in st.traverse():
                count += 1
        self.assertEqual(count, (len(self.mock_tree)*2)-1)
        self.assertEqual(4, len(node_paritions))

        # Test each leaf in a separate cluster
        p_clust.alpha = 0.2
        self.assertEqual(3, len(p_clust.partition_nodes(tree=deepcopy(self.mock_tree))))
        return

    def test_build_edge_node_index(self):
        from treesapp.phylo_cluster import PhyloClust
        p_clust = PhyloClust()
        node_edge_map = p_clust.build_edge_node_index(self.mock_tree)
        self.assertEqual(7, len(node_edge_map))
        self.assertEqual(0, node_edge_map['A'])
        self.assertEqual(6, node_edge_map['2'])
        return

    def test_split_node_partition_edges(self):
        from treesapp.phylo_cluster import PhyloClust
        p_clust = PhyloClust()
        p_clust.alpha = 0.1
        # Test empty
        edge_clusters = p_clust.split_node_partition_edges({})
        self.assertEqual(0, len(edge_clusters))

        a, b, c = TreeNode(name='A', dist=0.2), TreeNode(name='B', dist=0.4), TreeNode(name='0', dist=0.1)
        mock_partitions = {1: c, 2: b}
        c.children = [a]
        a.up = c
        self.assertEqual(2, len(mock_partitions))
        edge_clusters = p_clust.split_node_partition_edges(mock_partitions)
        self.assertEqual(3, len(edge_clusters))
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

    def test_define_tree_clusters(self):
        from treesapp.phylo_cluster import PhyloClust
        p_clust = PhyloClust()
        p_clust.alpha = 0.1
        cluster_map = p_clust.define_tree_clusters(tree=deepcopy(self.taxa_tree))
        self.assertTrue(len(self.taxa_tree) <= len(cluster_map) <= sum([1 for _ in self.taxa_tree.traverse()]))

        # Ensure each edge number appears in only one cluster
        edge_names = []
        for l in cluster_map.values():
            edge_names += l
        self.assertEqual(len(edge_names), len(set(edge_names)))
        return

    def test_cluster_phylogeny(self):
        from treesapp.phylo_cluster import cluster_phylogeny
        with pytest.raises(SystemExit):
            cluster_phylogeny([])
        return


if __name__ == '__main__':
    unittest.main()

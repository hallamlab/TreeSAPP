import os
import re
import unittest
import pytest
import shutil

from copy import deepcopy
from ete3 import Tree, TreeNode

from .testing_utils import get_test_data


class PhyloClusterTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.refpkg import ReferencePackage
        from treesapp.rel_evo_dist import RedTree
        self.tmp_dir = os.path.join(os.getcwd(), "tests", "phyloclust_test_dir")
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.mkdir(self.tmp_dir)

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

    def tearDown(self) -> None:
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
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
        self.assertEqual(3, len(node_paritions))

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
        p_clust.tax_rank = "family"
        p_clust.group_rel_dists(self.taxa_tree, self.refpkg.taxa_trie)
        return

    def test_match_edges_to_clusters(self):
        from treesapp.phylo_cluster import PhyloClust
        p_clust = PhyloClust()
        p_clust.alpha = 0.1
        p_clust.match_edges_to_clusters(tree=deepcopy(self.taxa_tree))

        # Ensure each edge number appears in only one cluster
        self.assertEqual(len(p_clust._edges_to_cluster_index), (len(self.taxa_tree)*2)-1)
        return

    def test_set_pquery_sample_name(self):
        from treesapp.phylo_cluster import PhyloClust
        from treesapp.phylo_seq import PQuery
        p_clust = PhyloClust()
        test_pq = PQuery()
        test_pq.seq_name = '3300019225.a:Ga0179949_11324471|McrA|1_124'
        p_clust.sample_re = re.compile(r'^(\d+)\.a:.*')
        p_clust.set_pquery_sample_name(test_pq, "test")
        self.assertEqual("3300019225", getattr(test_pq, "sample_name"))
        return

    def test_define_tree_clusters(self):
        from treesapp.phylo_cluster import PhyloClust
        p_clust = PhyloClust()
        p_clust.alpha = 0.05
        p_clust.define_tree_clusters(tree=deepcopy(self.taxa_tree))
        self.assertTrue(len(self.taxa_tree) <= len(p_clust.cluster_index) <= sum([1 for _ in self.taxa_tree.traverse()]))
        return

    def test_cluster_phylogeny(self):
        from treesapp.phylo_cluster import cluster_phylogeny
        # Should fail when multiple refpkgs provided
        with pytest.raises(SystemExit):
            cluster_phylogeny(["--refpkg_path",
                               self.refpkg.f__json, get_test_data(os.path.join("refpkgs", "McrB_build.pkl")),
                               "--jplace", get_test_data("epa_result.jplace")])

        # Test input is a single JPlace file
        cluster_phylogeny(["--refpkg_path", self.refpkg.f__json,
                           "--jplace", get_test_data("epa_result.jplace"),
                           "--output", self.tmp_dir])

        # Test input is a treesapp assign output directory
        cluster_phylogeny(["--refpkg_path", self.refpkg.f__json,
                           "--assign_output", get_test_data("test_output_TarA"),
                           "--output", self.tmp_dir,
                           "--alpha", str(0.4)])
        return


if __name__ == '__main__':
    unittest.main()

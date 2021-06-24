import os
import re
import unittest
import pytest
import shutil

from copy import deepcopy
from ete3 import Tree, TreeNode

from .testing_utils import get_test_data, random_ete_tree


class PhyloClusterTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.refpkg import ReferencePackage
        from treesapp.rel_evo_dist import RedTree
        self.tmp_dir = os.path.join(os.getcwd(), "tests", "phyloclust_test_dir")
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.mkdir(self.tmp_dir)

        self.refpkg = ReferencePackage()
        self.refpkg.f__pkl = get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
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

    def test_load_args(self):
        from treesapp.phylo_cluster import PhyloClust
        p_clust = PhyloClust()
        cli_args = ["-o", self.tmp_dir,
                    "--assign_output", get_test_data("test_output_TarA"), get_test_data("marker_test_results"),
                    "--refpkg_path", self.refpkg.f__pkl,
                    "--mode", "ref_guided",
                    "--tax_rank", "genus"]
        p_clust.load_args(p_clust.get_arguments(cli_args))
        self.assertEqual(os.path.join(self.tmp_dir, "intermediates", "load_pqueries/"),
                         p_clust.stage_output_dir)
        self.assertTrue(os.path.isdir(p_clust.stage_output_dir))
        self.assertTrue("mmseqs" in p_clust.executables)
        return

    def test_load_sample_placement_files(self):
        from treesapp import phylo_cluster
        p_clust = phylo_cluster.PhyloClust()
        # Simulate failure by de novo clustering with JPlace file
        with pytest.raises(SystemExit):
            p_clust.load_args(p_clust.arg_parser.parse_args(["--refpkg_path", self.refpkg.f__pkl,
                                                             "--mode", "de_novo",
                                                             "--jplace", get_test_data("epa_result.jplace")]))
            p_clust.load_sample_placement_files()
        return

    def test_partition_nodes(self):
        from treesapp.phylo_cluster import PhyloClust
        p_clust = PhyloClust()
        p_clust.alpha = 0.1
        node_paritions = p_clust.partition_nodes(tree=deepcopy(self.mock_tree), alpha=p_clust.alpha)
        # Ensure all the tree's nodes are present in the partitions
        count = 0
        for st in node_paritions.values():
            for _ in st.traverse():
                count += 1
        self.assertEqual(count, (len(self.mock_tree)*2)-1)
        self.assertEqual(3, len(node_paritions))

        # Test each leaf in a separate cluster
        p_clust.alpha = 0.2
        self.assertEqual(3, len(p_clust.partition_nodes(tree=deepcopy(self.mock_tree), alpha=p_clust.alpha)))
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
        edge_clusters = p_clust.split_node_partition_edges(p_clust.alpha, {})
        self.assertEqual(0, len(edge_clusters))

        a, b, c = TreeNode(name='A', dist=0.2), TreeNode(name='B', dist=0.4), TreeNode(name='0', dist=0.1)
        mock_partitions = {1: c, 2: b}
        c.children = [a]
        a.up = c
        self.assertEqual(2, len(mock_partitions))
        edge_clusters = p_clust.split_node_partition_edges(p_clust.alpha, mock_partitions)
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
        p_clust.normalize = True
        p_clust.tax_rank = "family"
        rel_dists = p_clust.group_rel_dists(self.taxa_tree, self.refpkg.taxa_trie)
        self.assertEqual(5, len(rel_dists))
        self.assertEqual(403, len(sum(rel_dists.values(), [])))
        return

    def test_calculate_distance_threshold(self):
        from treesapp import phylo_cluster
        p_clust = phylo_cluster.PhyloClust()
        p_clust.clustering_mode = "de_novo"
        alpha = p_clust.calculate_distance_threshold(taxa_tree=self.taxa_tree,
                                                     taxonomy=self.refpkg.taxa_trie)
        self.assertTrue(0.04 < round(alpha, 3) <= 0.08)

        # Define a subtree for a group of related taxa
        def clade_member(node: TreeNode) -> bool:
            t = getattr(node, "taxon")
            if "g__Methanosphaera" in [r.prefix_taxon() for r in t.lineage()]:
                return True
            return False

        subtree_leaves = self.taxa_tree.get_leaves(is_leaf_fn=clade_member)
        small_tree = subtree_leaves[0].get_common_ancestor(subtree_leaves)
        self.assertEqual(4, len(small_tree))

        self.assertEqual(-1, p_clust.calculate_distance_threshold(small_tree,
                                                                  self.refpkg.taxa_trie,
                                                                  override_rank="class"))
        p_clust.calculate_distance_threshold(small_tree, self.refpkg.taxa_trie, override_rank="species")

        p_clust.ref_pkg = self.refpkg
        p_clust.clustering_mode = "local"
        alpha = p_clust.calculate_distance_threshold(small_tree, self.refpkg.taxa_trie)
        self.assertTrue(0.90 < alpha < 1.0)
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

        p_clust.define_tree_clusters(tree=deepcopy(self.taxa_tree), override_alpha=0.1)
        self.assertEqual(0.05, p_clust.alpha)
        return

    def test_infer_cluster_phylogeny(self):
        from treesapp import phylo_cluster
        from treesapp import fasta
        p_clust = phylo_cluster.PhyloClust()
        p_clust.load_args(p_clust.get_arguments(["--refpkg_path", self.refpkg.f__pkl,
                                                 "--assign_output", get_test_data("test_output_TarA"),
                                                 "--output", self.tmp_dir,
                                                 "--alpha", str(0.4)]))
        q_seqs = fasta.FASTA(get_test_data("McrA_eval.faa"))
        q_seqs.load_fasta()
        _fa_input, ete_tree = phylo_cluster.infer_cluster_phylogeny(fa_input=q_seqs.file,
                                                                    executables=p_clust.executables,
                                                                    output_dir=p_clust.stage_output_dir)
        # Test for node names
        ex_name = "KKH90701"
        self.assertTrue(ex_name in ete_tree.get_leaf_names())
        self.assertEqual(os.path.join(self.tmp_dir, 'intermediates', 'load_pqueries/'), p_clust.stage_output_dir)

        # Test for number of nodes
        self.assertEqual(q_seqs.n_seqs(), len(ete_tree))
        return

    def test_assign_pqueries_to_leaf_clusters(self):
        from treesapp import phylo_cluster
        from treesapp import file_parsers
        from treesapp import assign
        from treesapp import seq_clustering
        p_clust = phylo_cluster.PhyloClust()
        pqueries = file_parsers.load_classified_sequences_from_assign_output(assign_output_dir=get_test_data("test_output_TarA"),
                                                                             assigner_cls=assign.Assigner(),
                                                                             refpkg_name="McrA")["McrA"]
        pquery_names = []
        for pq in pqueries:
            pquery_names.append(pq.place_name)
            setattr(pq, "sample_name", "test")

        i = 0
        spread = 3
        precluster_map = {}
        while pquery_names:
            p_otu = phylo_cluster.PhylOTU(name=str(i))
            p_otu.tree_node = random_ete_tree(pquery_names[0:spread])
            p_clust.cluster_index[str(i)] = p_otu
            if i % 2 == 0:
                mock_cluster = seq_clustering.Cluster(pquery_names[0])
                mock_cluster.members = pqueries[i*spread:(i*spread)+spread]
                precluster_map[i] = {str(i): mock_cluster}
            i += 1
            pquery_names = pquery_names[spread:]
        self.assertEqual([0, 2], list(precluster_map.keys()))
        self.assertEqual(4, len(p_clust.cluster_index))

        p_clust.assign_pqueries_to_leaf_clusters(pqueries=pqueries,
                                                 cluster_map=precluster_map)
        self.assertEqual({'0': 3, '1': 3, '2': 3, '3': 2},
                         p_clust.sample_mat["test"])

        return

    def test_format_precluster_map(self):
        from treesapp import phylo_cluster
        precluster_map = {0: "Cluster", 1: "Cluster", 2: "Cluster"}
        phylo_groups = phylo_cluster.format_precluster_map(cluster_method="align", precluster_map=precluster_map)
        self.assertEqual(1, len(phylo_groups))
        self.assertEqual(3, len(phylo_groups[0]))
        return

    def test_de_novo_phylo_clusters(self):
        from treesapp import phylo_cluster
        from treesapp import rel_evo_dist
        p_clust = phylo_cluster.PhyloClust()
        p_clust.load_args(p_clust.get_arguments(["--refpkg_path", self.refpkg.f__pkl,
                                                 "--assign_output", get_test_data("test_output_TarA"),
                                                 "--output", self.tmp_dir,
                                                 "--alpha", str(0.4)]))
        p_clust.load_sample_placement_files()
        # Test with placement space clustering
        taxa_tree = p_clust.ref_pkg.taxonomically_label_tree()
        red_tree = rel_evo_dist.RedTree()
        red_tree.decorate_rel_dist(taxa_tree)

        phylo_cluster.de_novo_phylo_clusters(p_clust, taxa_tree, cluster_method="psc")
        self.assertEqual(11, len(p_clust.clustered_pqueries))
        self.assertEqual(len(p_clust.clustered_pqueries),
                         len(set([pq.p_otu for pq in p_clust.clustered_pqueries])))
        # Test with pairwise sequence clustering
        p_clust.cluster_index.clear()
        p_clust.clustered_pqueries.clear()
        p_clust._edges_to_cluster_index.clear()
        phylo_cluster.de_novo_phylo_clusters(p_clust, taxa_tree, cluster_method="align", drep_id=0.8)
        self.assertEqual(10, len(set([pq.p_otu for pq in p_clust.clustered_pqueries])))

        phylo_cluster.de_novo_phylo_clusters(p_clust, taxa_tree, cluster_method="psc", phylo_group=None)

        return

    def test_cluster_by_local_alignment(self):
        from treesapp import phylo_cluster
        p_clust = phylo_cluster.PhyloClust()
        p_clust.load_args(p_clust.get_arguments(["--refpkg_path", self.refpkg.f__pkl,
                                                 "--assign_output", get_test_data("test_output_TarA"),
                                                 "--output", self.tmp_dir]))
        p_clust.load_sample_placement_files()
        phylo_cluster.cluster_by_local_alignment(phylo_clust=p_clust)
        self.assertEqual(11, len(p_clust.clustered_pqueries))
        self.assertEqual(10, len(p_clust.cluster_index))
        self.assertEqual(10, len(p_clust.sample_mat["test_output_TarA"]))

        return

    def test_subtree_finder(self):
        from treesapp import phylo_cluster
        target_leaves = {'9_McrA', '6_McrA', '7_McrA',
                         '37_McrA', '49_McrA', '27_McrA',
                         '191_McrA', '192_McrA', '194_McrA'}
        subtree_size = 3
        sub_root = phylo_cluster.subtree_finder(ref_tree=self.refpkg.get_ete_tree(),
                                                leaf_nodes=target_leaves,
                                                tree_size=subtree_size)
        self.assertTrue(set(sub_root.get_leaf_names()).issubset(target_leaves))
        self.assertEqual(subtree_size, len(sub_root.get_leaf_names()))
        return

    def test_select_subtree_sequences(self):
        from treesapp import phylo_cluster
        from treesapp.seq_clustering import Cluster
        from treesapp.phylo_seq import PQuery, PhyloPlace
        mock_cluster = Cluster("TarA")
        test_placement_edges = [62, 60, 65, 63]
        mock_placements = []
        while test_placement_edges:
            pplace = PhyloPlace()
            pplace.edge_num = test_placement_edges.pop()
            mock_placements.append(pplace)
        while mock_placements:
            pq = PQuery()
            pq.consensus_placement = mock_placements.pop()
            mock_cluster.members.append(pq)

        n = 8
        ref_fasta = phylo_cluster.select_subtree_sequences(ref_pkg=self.refpkg,
                                                           clusters=[mock_cluster],
                                                           subtree_size=n)
        self.assertEqual(n+1, ref_fasta.n_seqs())
        return

    def test_get_outgroup(self):
        from treesapp import phylo_cluster
        outgroup = phylo_cluster.get_outgroup(tree=self.mock_tree, target="D")
        self.assertIsInstance(outgroup, TreeNode)
        self.assertEqual('A', outgroup.name)
        return

    def test_cluster_phylogeny(self):
        from treesapp.phylo_cluster import cluster_phylogeny
        # Should fail when multiple refpkgs provided
        with pytest.raises(SystemExit):
            cluster_phylogeny(["--refpkg_path",
                               self.refpkg.f__pkl, get_test_data(os.path.join("refpkgs", "McrB_build.pkl")),
                               "--jplace", get_test_data("epa_result.jplace")])

        # TODO: Remove once PSC with ref-guided is implemented
        cluster_phylogeny(["--refpkg_path", self.refpkg.f__pkl,
                           "--jplace", get_test_data("epa_result.jplace"),
                           "--output", self.tmp_dir,
                           "--mode", "ref_guided"])

        # Test input is a treesapp assign output directory
        cluster_phylogeny(["--refpkg_path", self.refpkg.f__pkl,
                           "--assign_output", get_test_data("test_output_TarA"),
                           "--output", self.tmp_dir,
                           "--alpha", str(0.4),
                           "--mode", "de_novo",
                           "--delete"])
        self.assertFalse(os.path.isdir(os.path.join(self.tmp_dir, "intermediates")))

        # Test input is a single JPlace file
        cluster_phylogeny(["--refpkg_path", self.refpkg.f__pkl,
                           "--jplace", get_test_data("epa_result.jplace"),
                           "--output", self.tmp_dir])

        self.assertTrue(os.path.isfile(os.path.join(self.tmp_dir, "final_outputs", "phylotu_taxa.tsv")))
        self.assertTrue(os.path.isfile(os.path.join(self.tmp_dir, "final_outputs", "phylotu_matrix.tsv")))
        self.assertTrue(os.path.isfile(os.path.join(self.tmp_dir, "final_outputs", "phylotu_pquery_assignments.tsv")))
        with open(os.path.join(self.tmp_dir, "final_outputs", "phylotu_pquery_assignments.tsv")) as asns:
            lines = asns.readlines()
        self.assertEqual(13, len(lines))

        # Test local alignment clustering
        cluster_phylogeny(["--refpkg_path", self.refpkg.f__pkl,
                           "--assign_output", get_test_data("test_output_TarA"),
                           "--output", self.tmp_dir,
                           "--mode", "local",
                           "--delete"])
        return


if __name__ == '__main__':
    unittest.main()

import unittest
import pytest
import os
from .testing_utils import get_test_data


@pytest.fixture(scope="class")
def test_pquery(request):
    from treesapp.jplace_utils import jplace_parser, demultiplex_pqueries
    from treesapp.entish import map_internal_nodes_leaves
    jplace_data = jplace_parser(get_test_data("epa_result.jplace"))
    pqueries = demultiplex_pqueries(jplace_data)
    request.cls.db = pqueries.pop(9)
    request.cls.db.node_map = map_internal_nodes_leaves(jplace_data.tree)
    request.cls.db.abundance = 2.5
    return


@pytest.mark.usefixtures("test_pquery")
class PhyloSeqTests(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.phylo_seq import PQuery
        from treesapp.refpkg import ReferencePackage
        # A generic placement dictionary parsed from a JPlace file
        self.placement_dict = {'p': [[489, -50.7, 0.6, 0.859, 1.227],
                                     [1, -50.8, 0.4, 0.1, 1.1]],
                               'n': ['seq_test_1']}
        self.field_order = ['edge_num', 'likelihood', 'like_weight_ratio', 'distal_length', 'pendant_length']

        # Prepare the various PQuery instances for tests
        self.pquery_test_1 = PQuery()
        self.pquery_test_1.placements = self.placement_dict
        self.pquery_test_2 = self.db
        self.pquery_test_3 = PQuery()
        self.pquery_test_3.placements = {'p': [[245, -1000.2, 0.98, 0.1, 0.1]],
                                         'n': ['seq_test_3']}
        self.pquery_test_4 = PQuery()
        self.pquery_test_4.placements = {'p': [[454, -56278.9228361584, 0.3661573002, 0.0153492795, 0.0220553182],
                                               [452, -56278.3740983653, 0.6338426998, 0.0076035529, 0.0233278435]],
                                         'n': ['PKL62129.1 [Methanomicrobiales archaeon HGW-Methanomicrobiales-2]']}
        self.pquery_test_5 = PQuery()
        self.pquery_test_5.placements = {'p': [[397, -56289.2785, 0.3114285092, 0.0267577444, 0.068344758],
                                               [395, -56288.7317004902, 0.0327271581, 0.0100713972, 0.0666057481],
                                               [393, -56289.2785075033, 0.0189423138, 0.0153064082, 0.0729607682],
                                               [388, -56285.7632899202, 0.6369020189, 0.0103030924, 0.0686955093]],
                                         'n': ['AFD09581.1']}
        for pq in [self.pquery_test_1, self.pquery_test_2, self.pquery_test_3, self.pquery_test_4, self.pquery_test_5]:
            pq.node_map = self.db.node_map

        # Reference packages
        puha_rp = ReferencePackage()
        puha_rp.f__pkl = get_test_data(os.path.join("refpkgs", "PuhA_build.pkl"))
        puha_rp.slurp()
        node_map = puha_rp.get_internal_node_leaf_map()

        # PQueries
        pquery_1 = PQuery()
        pquery_2 = PQuery()
        pquery_1.seq_name, pquery_1.start, pquery_1.end, pquery_1.ref_name, pquery_1.node_map = "seq1", 1, 112, "PuhA", node_map
        pquery_2.seq_name, pquery_2.start, pquery_2.end, pquery_2.ref_name, pquery_2.node_map = "seq2", 3, 184, "PuhA", node_map
        self.pqueries = {"PuhA": [pquery_1, pquery_2]}

        # Prepare the reference package
        self.refpkg = ReferencePackage(refpkg_name="McrA")
        self.refpkg.f__pkl = get_test_data("refpkgs/McrA_build.pkl")
        self.refpkg.slurp()
        return

    def test_abundify_tree_saps(self):
        from treesapp import phylo_seq
        # Set up the input data
        abund_dict = {"seq1|PuhA|1_112": 100,
                      "seq2|PuhA|3_184": 120,
                      "seq2|NxrA|2_210": 80}
        pquery_1, pquery_2 = self.pqueries["PuhA"]

        # Test when no names map
        phylo_seq.quantify_pquery_instances(tree_saps=self.pqueries, abundance_dict=abund_dict)
        self.assertEqual(0.0, pquery_2.abundance)

        # Set the place name so the names now match
        pquery_1.place_name = "{}|{}|{}_{}".format(pquery_1.seq_name, pquery_1.ref_name, pquery_1.start, pquery_1.end)
        phylo_seq.quantify_pquery_instances(tree_saps=self.pqueries, abundance_dict=abund_dict)
        self.assertEqual(100, pquery_1.abundance)

        return

    def test_phyloplace_str(self):
        from treesapp import phylo_seq
        pplace = phylo_seq.PhyloPlace(placement={"n": "test_seq_1", "p": [["13", "0.01", "0.02"]]},
                                      field_positions=["edge_num", "distal_length", "pendant_length"])
        self.assertEqual("Placement of sequence 'test_seq_1' on edge 13", str(pplace))
        return

    def test_calc_mean_tip_length(self):
        placement = self.pquery_test_2.placements[0]
        placement.calc_mean_tip_length(internal_leaf_node_map=self.pquery_test_2.node_map, ref_tree=self.refpkg.tree)
        self.assertTrue(0.0 < placement.mean_tip_length)
        return

    def test_children_lineage(self):
        self.pquery_test_2.process_max_weight_placement(self.refpkg.taxonomically_label_tree())
        with pytest.raises(SystemExit):
            self.pquery_test_2.children_lineage(leaves_taxa_map={'1': 'd__Archaea', '2': 'd__Archaea',
                                                                 '100': 'd__Archaea; p__Euryarchaeota;'
                                                                        ' c__Methanobacteria; o__Methanobacteriales;'
                                                                        ' f__Methanobacteriaceae; g__Methanobrevibacter'})
        return

    def test_process_max_weight_placement(self):
        from treesapp.phylo_seq import split_placements, PhyloPlace
        # Trigger a failure
        self.pquery_test_1.placements = split_placements({'p': [[511, -50.7, 0.7, 0.859, 1.227],
                                                                [2, -50.8, 0.3, 0.1, 1.1]],
                                                          'n': ['seq_test_1']})
        with pytest.raises(SystemExit):
            self.pquery_test_1.process_max_weight_placement(self.refpkg.taxonomically_label_tree())

        # Use valid placements to succeed
        self.pquery_test_1.placements = split_placements(self.placement_dict)
        self.pquery_test_1.process_max_weight_placement(self.refpkg.taxonomically_label_tree())
        self.assertEqual(2, len(self.pquery_test_1.placements))
        self.assertIsInstance(self.pquery_test_1.consensus_placement, PhyloPlace)
        self.assertEqual(0.6, self.pquery_test_1.consensus_placement.like_weight_ratio)
        return

    def test_calculate_consensus_placement(self):
        from treesapp.phylo_seq import split_placements
        # Exit if the placements dictionary of placement strings hasn't been converted to PhyloPlace instances
        with pytest.raises(AttributeError):
            self.pquery_test_1.calculate_consensus_placement(self.refpkg.taxonomically_label_tree())
        # Now format and test for a proper PQuery with multiple placements near the root
        self.pquery_test_1.placements = split_placements(self.placement_dict)
        self.pquery_test_1.calculate_consensus_placement(self.refpkg.taxonomically_label_tree())
        self.assertEqual("r__Root", self.pquery_test_1.lct)
        self.assertEqual(1, len(self.pquery_test_1.placements))
        self.assertEqual(1.1, self.pquery_test_1.consensus_placement.pendant_length)
        self.assertEqual(498, self.pquery_test_1.consensus_placement.edge_num)

        # Test for a PQuery with just a single placement
        self.pquery_test_3.placements = split_placements(self.pquery_test_3.placements)
        self.pquery_test_3.calculate_consensus_placement(self.refpkg.taxonomically_label_tree())
        self.assertEqual(245, self.pquery_test_3.consensus_placement.edge_num)

        # Test with multiple placements where the LCA is one of the placements
        self.pquery_test_4.placements = split_placements(self.pquery_test_4.placements)
        self.pquery_test_4.calculate_consensus_placement(self.refpkg.taxonomically_label_tree())
        self.assertEqual(454, self.pquery_test_4.consensus_placement.edge_num)

        # Test with multiple placements and the LCA is none of the placements
        self.pquery_test_5.placements = split_placements(self.pquery_test_5.placements)
        self.pquery_test_5.calculate_consensus_placement(self.refpkg.taxonomically_label_tree(), min_aelw=0.7)
        self.assertEqual(397, self.pquery_test_5.consensus_placement.edge_num)

        # Ensure none of the likelihood weight ratios exceed 1
        for pq in [self.pquery_test_1, self.pquery_test_3, self.pquery_test_4, self.pquery_test_5]:
            self.assertTrue(pq.consensus_placement.like_weight_ratio <= 1.0)
        return

    def test_assignments_to_pqueries(self):
        from treesapp.phylo_seq import assignments_to_pqueries
        assignment_lines = [['test_TarA.1', 'scaffold_5431_c1_4', 'DsrAB', '79', '161', 'r__Root', '0.0', '282', '7.3e-10', '1.0', '3.085', '0.599,1.952,0.534'],
                            ['test_TarA.1', 'scaffold_59587_c1_1', 'DsrAB', '1', '182', 'r__Root; d__Bacteria; p__Proteobacteria; c__Deltaproteobacteria', '0.0', '979', '2.1e-63', '1.0', '0.143', '0.022,0.121,0.0']]
        pqueries = assignments_to_pqueries(assignment_lines)
        self.assertTrue("DsrAB" in pqueries)
        self.assertEqual(['scaffold_5431_c1_4|DsrAB|79_161', 'scaffold_59587_c1_1|DsrAB|1_182'],
                         [x.place_name for x in pqueries["DsrAB"]])
        self.assertEqual(['scaffold_5431_c1_4', 'scaffold_59587_c1_1'],
                         [x.seq_name for x in pqueries["DsrAB"]])

        # Fail due to bad line format
        with pytest.raises(SystemExit):
            assignments_to_pqueries(["Taxonomy\tAbundance\tiNode\tE-value\tLWR\tEvoDist\tDistances".split("\t")])
        return

    def test_phylo_place(self):
        from treesapp.phylo_seq import PhyloPlace
        bad_dict = {}
        for k, v in self.placement_dict.items():
            bad_dict[k] = v.copy()
        bad_dict['n'].append("second seq name")
        with pytest.raises(SystemExit):
            PhyloPlace(bad_dict)
        pplace = PhyloPlace(self.placement_dict)
        self.assertEqual(pplace.like_weight_ratio, 0.6)
        self.assertEqual(pplace.distal_length, 0.859)
        return

    def test_name_placed_sequence(self):
        # Test when the seq_name attribute is empty

        # Test when the seq_name attribute is set - nothing should be changed

        return

    def test_sort_centroids_from_clusters(self):
        from treesapp import phylo_seq
        with pytest.raises(AttributeError):
            phylo_seq.sort_centroids_from_clusters(list(range(8)), [0]*2 + [1]*6)
        pq_1 = phylo_seq.PQuery()
        pq_2 = phylo_seq.PQuery()
        final_clusters = phylo_seq.sort_centroids_from_clusters(pqueries=[pq_1, pq_2],
                                                                cluster_indices=[0, 0])
        self.assertIsInstance(final_clusters, list)
        self.assertEqual(1, len(final_clusters))
        return

    def test_cluster_pquery_distances(self):
        from treesapp import phylo_seq
        from treesapp import file_parsers
        from treesapp.seq_clustering import Cluster
        from treesapp.assign import Assigner
        refpkg_pquery_map = file_parsers.load_classified_sequences_from_assign_output(get_test_data("test_output_TarA"),
                                                                                      Assigner(),
                                                                                      refpkg_name="McrA")
        test_pqueries = refpkg_pquery_map["McrA"]
        self.assertEqual(11, len(test_pqueries))

        pquery_clusters = phylo_seq.cluster_pquery_placement_space_distances(test_pqueries)
        self.assertEqual(11, len(pquery_clusters))
        self.assertIsInstance(pquery_clusters['0'], Cluster)

        phylo_seq.cluster_pquery_placement_space_distances(test_pqueries, min_cluster_size=2)
        self.assertEqual(1, max([len(i.members) for i in pquery_clusters.values()]))
        self.assertEqual(11, len(pquery_clusters))
        return

    def test_split_placements(self):
        from treesapp.phylo_seq import split_placements
        self.pquery_test_1.placements = split_placements(self.placement_dict)
        self.pquery_test_1.name_placed_sequence()
        self.assertEqual(2, len(self.pquery_test_1.placements))
        self.assertEqual("seq_test_1", self.pquery_test_1.place_name)
        return

    def test_sum_abundances_per_node(self):
        self.pquery_test_2.consensus_placement = self.pquery_test_2.placements[0]
        # Test the exception
        self.pquery_test_2.abundance = None
        self.assertEqual(0, self.pquery_test_2.sum_abundances_per_node(leaf_abundance_sums={})['121_McrA'])

        self.pquery_test_2.abundance = 10
        leaf_abundance_sums = self.pquery_test_2.sum_abundances_per_node(leaf_abundance_sums={})
        self.assertEqual(5, len(leaf_abundance_sums))
        self.assertEqual(2, leaf_abundance_sums['121_McrA'])
        return


if __name__ == '__main__':
    unittest.main()

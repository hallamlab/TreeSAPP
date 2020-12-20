import unittest
import pytest
from .testing_utils import get_test_data


@pytest.fixture(scope="class")
def test_pquery(request):
    from treesapp.jplace_utils import jplace_parser, demultiplex_pqueries
    from treesapp.entish import map_internal_nodes_leaves
    jplace_data = jplace_parser(get_test_data("epa_result.jplace"))
    pqueries = demultiplex_pqueries(jplace_data)
    request.cls.db = pqueries.pop(0)
    request.cls.db.node_map = map_internal_nodes_leaves(jplace_data.tree)
    request.cls.db.abundance = 2.5
    return


@pytest.mark.usefixtures("test_pquery")
class PhyloSeqtests(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.phylo_seq import PQuery
        from treesapp.refpkg import ReferencePackage
        # A generic placement dictionary parsed from a JPlace file
        self.placement_dict = {'p': [[489, -50.7, 0.7, 0.859, 1.227],
                                     [1, -50.8, 0.3, 0.1, 1.1]],
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

        # Prepare the reference package
        self.refpkg = ReferencePackage(refpkg_name="McrA")
        self.refpkg.f__json = get_test_data("refpkgs/McrA_build.pkl")
        self.refpkg.slurp()
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
        self.assertEqual(0.7, self.pquery_test_1.consensus_placement.like_weight_ratio)
        return

    def test_calculate_consensus_placement(self):
        from treesapp.phylo_seq import split_placements
        # Exit if the placements dictionary of placement strings hasn't been converted to PhyloPlace instances
        with pytest.raises(AttributeError):
            self.pquery_test_1.calculate_consensus_placement(self.refpkg.taxonomically_label_tree())
        # Now format and test for a proper PQuery with multiple placements near the root
        self.pquery_test_1.placements = split_placements(self.placement_dict)
        self.pquery_test_1.calculate_consensus_placement(self.refpkg.taxonomically_label_tree(), min_aelw=0.6)
        self.assertEqual("r__Root", self.pquery_test_1.lct)
        self.assertEqual(2, len(self.pquery_test_1.placements))
        self.assertEqual(490, self.pquery_test_1.consensus_placement.edge_num)

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
        self.pquery_test_5.calculate_consensus_placement(self.refpkg.taxonomically_label_tree())
        self.assertEqual(398, self.pquery_test_5.consensus_placement.edge_num)

        # Ensure none of the likelihood weight ratios exceed 1
        for pq in [self.pquery_test_1, self.pquery_test_3, self.pquery_test_4, self.pquery_test_5]:
            self.assertTrue(pq.consensus_placement.like_weight_ratio <= 1.0)
        return

    def test_assignments_to_treesaps(self):
        from treesapp.phylo_seq import assignments_to_treesaps
        assignment_lines = [['test_TarA.1', 'scaffold_5431_c1_4', 'DsrAB', '79', '161', 'r__Root', '0.0', '282', '7.3e-10', '1.0', '3.085', '0.599,1.952,0.534'],
                            ['test_TarA.1', 'scaffold_59587_c1_1', 'DsrAB', '1', '182', 'r__Root; d__Bacteria; p__Proteobacteria; c__Deltaproteobacteria', '0.0', '979', '2.1e-63', '1.0', '0.143', '0.022,0.121,0.0']]
        pqueries = assignments_to_treesaps(assignment_lines)
        self.assertTrue("DsrAB" in pqueries)
        self.assertEqual(['scaffold_5431_c1_4|DsrAB|79_161', 'scaffold_59587_c1_1|DsrAB|1_182'],
                         [x.place_name for x in pqueries["DsrAB"]])
        self.assertEqual(['scaffold_5431_c1_4', 'scaffold_59587_c1_1'],
                         [x.seq_name for x in pqueries["DsrAB"]])

        # Fail due to bad line format
        with pytest.raises(SystemExit):
            assignments_to_treesaps(["Taxonomy\tAbundance\tiNode\tE-value\tLWR\tEvoDist\tDistances".split("\t")])
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
        self.assertEqual(pplace.like_weight_ratio, 0.7)
        self.assertEqual(pplace.distal_length, 0.859)
        return

    def test_name_placed_sequence(self):
        # Test when the seq_name attribute is empty

        # Test when the seq_name attribute is set - nothing should be changed

        return

    def test_split_placements(self):
        from treesapp.phylo_seq import split_placements
        self.pquery_test_1.placements = split_placements(self.placement_dict)
        self.pquery_test_1.name_placed_sequence()
        self.assertEqual(2, len(self.pquery_test_1.placements))
        self.assertEqual("seq_test_1", self.pquery_test_1.place_name)
        return

    def test_sum_rpkms_per_node(self):
        self.pquery_test_2.consensus_placement = self.pquery_test_2.placements[0]
        # Test the exception
        self.pquery_test_2.abundance = None
        self.assertEqual(0, self.pquery_test_2.sum_rpkms_per_node(leaf_rpkm_sums={})['21_McrA'])

        self.pquery_test_2.abundance = 12
        leaf_rpkm_sums = self.pquery_test_2.sum_rpkms_per_node(leaf_rpkm_sums={})
        self.assertEqual(12, len(leaf_rpkm_sums))
        self.assertEqual(1, leaf_rpkm_sums['21_McrA'])
        return


if __name__ == '__main__':
    unittest.main()

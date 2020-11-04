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
        self.placement_dict = {'p': [[33, -50.7, 0.7, 0.859, 1.227],
                                     [2, -50.8, 0.3, 0.1, 1.1]],
                               'n': ['seq_test_1']}
        self.field_order = ['edge_num', 'likelihood', 'like_weight_ratio', 'distal_length', 'pendant_length']
        self.pquery_test_1 = PQuery()
        self.pquery_test_1.placements = self.placement_dict
        self.pquery_test_2 = self.db
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
        with pytest.raises(AttributeError):
            self.pquery_test_1.calculate_consensus_placement(self.refpkg.taxonomically_label_tree())
        self.pquery_test_1.placements = split_placements(self.placement_dict)
        self.pquery_test_1.calculate_consensus_placement(self.refpkg.taxonomically_label_tree())
        self.assertEqual("r__Root", self.pquery_test_1.lct)
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

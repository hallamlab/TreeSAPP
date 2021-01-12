import os
import shutil
import unittest
import pytest

from collections import namedtuple
from .testing_utils import get_test_data


class AssignerTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp import refpkg
        from treesapp import phylo_seq
        from .testing_utils import get_test_data
        self.output_dir = "./tests/assigner_test/"
        self.mock_args = namedtuple('args', ['reads', 'reverse', 'paired', 'rpkm',
                                             'molecule', 'overwrite', 'stage', 'output'])
        # Reference packages
        puha_rp = refpkg.ReferencePackage()
        puha_rp.f__json = get_test_data(os.path.join("refpkgs", "PuhA_build.pkl"))
        puha_rp.slurp()
        node_map = puha_rp.get_internal_node_leaf_map()

        # PQueries
        pquery_1 = phylo_seq.PQuery()
        pquery_2 = phylo_seq.PQuery()
        pquery_1.seq_name, pquery_1.start, pquery_1.end, pquery_1.ref_name, pquery_1.node_map = "seq1", 1, 112, "PuhA", node_map
        pquery_2.seq_name, pquery_2.start, pquery_2.end, pquery_2.ref_name, pquery_2.node_map = "seq2", 3, 184, "PuhA", node_map
        self.pqueries = {"PuhA": [pquery_1, pquery_2]}

        self.refpkg_dict = {"PuhA": puha_rp}

        # Ensure the output directory exists
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        return

    def tearDown(self) -> None:
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

    def test_abundify_tree_saps(self):
        from treesapp import assign
        # Set up the input data
        abund_dict = {"seq1|PuhA|1_112": 100,
                      "seq2|PuhA|3_184": 120,
                      "seq2|NxrA|2_210": 80}
        pquery_1, pquery_2 = self.pqueries["PuhA"]

        # Test when no names map
        assign.abundify_tree_saps(tree_saps=self.pqueries, abundance_dict=abund_dict)
        self.assertEqual(0.0, pquery_2.abundance)

        # Set the place name so the names now match
        pquery_1.place_name = "{}|{}|{}_{}".format(pquery_1.seq_name, pquery_1.ref_name, pquery_1.start, pquery_1.end)
        assign.abundify_tree_saps(tree_saps=self.pqueries, abundance_dict=abund_dict)
        self.assertEqual(100, pquery_1.abundance)

        return

    def test_select_consensus_placements(self):
        from treesapp.assign import select_query_placements
        with pytest.raises(ValueError):
            select_query_placements(self.pqueries, self.refpkg_dict, "red")
        return

    def test_determine_confident_lineage(self):
        from treesapp.assign import determine_confident_lineage
        from treesapp.phylo_seq import PhyloPlace
        # Create test data
        pquery_1, pquery_2 = self.pqueries["PuhA"]

        # test when mode is 'aelw' and the lineage is not adjusted
        with pytest.raises(SystemExit):
            determine_confident_lineage(tree_saps=self.pqueries, refpkg_dict=self.refpkg_dict, mode="aelw")

        # Instantiate the consensus placement attributes
        pquery_1.consensus_placement = PhyloPlace()
        pquery_1.consensus_placement.edge_num = 46
        pquery_1.lct = "r__Root; d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales; f__Rhodobacteraceae; g__Rhodobacter; s__Rhodobacter sphaeroides"
        pquery_1.avg_evo_dist = 3.5
        pquery_2.consensus_placement = PhyloPlace()
        pquery_2.consensus_placement.edge_num = 1
        pquery_2.lct = "d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales"
        determine_confident_lineage(tree_saps=self.pqueries, refpkg_dict=self.refpkg_dict, mode="aelw")
        self.assertEqual(pquery_1.lct, pquery_1.recommended_lineage)
        self.assertEqual(pquery_2.recommended_lineage, "r__Root; d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales")

        # test when mode is 'max_lwr' and linear model is used
        determine_confident_lineage(tree_saps=self.pqueries, refpkg_dict=self.refpkg_dict, mode="max_lwr")
        self.assertEqual(pquery_1.recommended_lineage, "r__Root")
        return

    def test_get_info(self):
        from treesapp import assign
        ts_assigner = assign.Assigner()
        info_string = ts_assigner.get_info()
        self.assertTrue(info_string)
        return

    def test_decide_stage(self):
        from treesapp import assign
        ts_assigner = assign.Assigner()
        ts_assigner.output_dir = self.output_dir
        ts_assigner.var_output_dir = self.output_dir
        ts_assigner.final_output_dir = self.output_dir
        # Test normal condition for protein
        ts_assigner.decide_stage(self.mock_args(molecule="prot", overwrite=False, stage="continue",
                                                output=self.output_dir, reads="", reverse="", rpkm=False, paired="pe"))
        self.assertTrue(ts_assigner.stage_status("clean"))
        self.assertFalse(ts_assigner.stage_status("orf-call"))
        self.assertFalse(ts_assigner.stage_status("abundance"))

        # Test a nucleotide input but with no rpkm
        ts_assigner.change_stage_status("orf-call", True)
        ts_assigner.change_stage_status("abundance", True)
        ts_assigner.decide_stage(self.mock_args(molecule="dna", overwrite=False, stage="continue",
                                                output=self.output_dir, reads="", reverse="", rpkm=False, paired="pe"))
        self.assertTrue(ts_assigner.stage_status("orf-call"))
        self.assertFalse(ts_assigner.stage_status("abundance"))

        # Test an RPKM run with bad fastq files
        with pytest.raises(SystemExit):
            ts_assigner.decide_stage(
                self.mock_args(molecule="dna", overwrite=False, stage="continue", output=self.output_dir,
                               reads="test_TarA.1.fq", reverse="", rpkm=True, paired="pe"))

        # Test a successful RPKM run
        ts_assigner.decide_stage(self.mock_args(molecule="dna", overwrite=False,
                                                stage="continue", output=self.output_dir,
                                                reads=get_test_data("test_TarA.1.fq"),
                                                reverse=get_test_data("test_TarA.2.fq"),
                                                rpkm=True, paired="pe"))
        self.assertTrue(ts_assigner.stage_status("abundance"))
        return


if __name__ == '__main__':
    unittest.main()

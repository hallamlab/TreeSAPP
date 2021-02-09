import os
import shutil
import unittest
import pytest

from .testing_utils import get_test_data


class AssignerTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp import refpkg
        from treesapp import phylo_seq
        from treesapp import treesapp_args
        from .testing_utils import get_test_data
        self.output_dir = "./tests/assigner_test/"

        # Simulate the treesapp assign arguments
        self.mock_parser = treesapp_args.TreeSAPPArgumentParser()
        treesapp_args.add_classify_arguments(self.mock_parser)

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

        # Test normal condition for protein
        args = self.mock_parser.parse_args(["--fastx_input", get_test_data("marker_test_suite.faa"),
                                            "--molecule", "prot",
                                            "--refpkg_dir", self.output_dir,
                                            "--output", self.output_dir])
        ts_assigner.furnish_with_arguments(args)
        ts_assigner.check_previous_output(overwrite=False)
        ts_assigner.check_classify_arguments(args)
        ts_assigner.decide_stage(args)
        self.assertTrue(ts_assigner.stage_status("clean"))
        self.assertFalse(ts_assigner.stage_status("orf-call"))
        self.assertFalse(ts_assigner.stage_status("abundance"))

        # Test a nucleotide input but with no abundance
        ts_assigner.change_stage_status("orf-call", True)
        ts_assigner.change_stage_status("abundance", True)
        args = self.mock_parser.parse_args(["--fastx_input", get_test_data("marker_test_suite.faa"),
                                            "--molecule", "dna",
                                            "--refpkg_dir", self.output_dir,
                                            "--output", self.output_dir])
        ts_assigner.furnish_with_arguments(args)
        ts_assigner.decide_stage(args)
        self.assertTrue(ts_assigner.stage_status("orf-call"))
        self.assertFalse(ts_assigner.stage_status("abundance"))

        # Test an RPKM run with bad fastq files
        with pytest.raises(SystemExit):
            args = self.mock_parser.parse_args(["--fastx_input", get_test_data("marker_test_suite.faa"),
                                                "--molecule", "dna",
                                                "--refpkg_dir", self.output_dir,
                                                "--output", self.output_dir,
                                                "--reads", "test_TarA.1.fq", "--rel_abund"])
            ts_assigner.furnish_with_arguments(args)
            ts_assigner.decide_stage(args)

        # Test a successful RPKM run
        args = self.mock_parser.parse_args(["--fastx_input", get_test_data("marker_test_suite.faa"),
                                            "--molecule", "dna",
                                            "--refpkg_dir", self.output_dir,
                                            "--output", self.output_dir,
                                            "--reads", get_test_data("test_TarA.1.fq"),
                                            "--reverse", get_test_data("test_TarA.2.fq"),
                                            "--rel_abund", "--pairing", "pe"])
        ts_assigner.furnish_with_arguments(args)
        ts_assigner.decide_stage(args)
        self.assertTrue(ts_assigner.stage_status("abundance"))
        return

    def test_fetch_hmmsearch_outputs(self):
        from treesapp import assign
        ts_assigner = assign.Assigner()
        ts_assigner.output_dir = get_test_data("marker_test_results")
        ts_assigner.var_output_dir = os.path.join(ts_assigner.output_dir, "intermediates")
        for stage_order, stage in ts_assigner.stages.items():
            stage.dir_path = ts_assigner.var_output_dir + os.sep + stage.name + os.sep

        # Fail due to wrong stage
        ts_assigner.current_stage = ts_assigner.stages[0]
        with pytest.raises(SystemExit):
            ts_assigner.fetch_hmmsearch_outputs({"McrA"})

        ts_assigner.current_stage = ts_assigner.stage_lookup("search")
        ts_assigner.stage_output_dir = ts_assigner.current_stage.dir_path
        # Reference package targets don't match
        self.assertEqual([], ts_assigner.fetch_hmmsearch_outputs({"DsrAB"}))

        # Intended functionality test
        hmm_domtbls = ts_assigner.fetch_hmmsearch_outputs({"McrA", "McrB"})
        self.assertEqual(2, len(hmm_domtbls))

        return

    def test_gather_split_msa(self):
        from treesapp import assign
        test_refpkg_names = ["McrA", "McrB", "DsrAB"]
        aln_path = os.path.join("marker_test_results", "intermediates", "align") + os.sep

        split_msa_files = assign.gather_split_msa(align_dir=get_test_data(aln_path),
                                                  refpkg_names=test_refpkg_names)

        self.assertEqual(["McrA", "McrB"], list(split_msa_files.keys()))
        self.assertEqual("McrA_hmm_purified_group0-BMGE_references.mfa",
                         os.path.basename(split_msa_files["McrA"][0].ref))
        self.assertEqual("McrA_hmm_purified_group0-BMGE_queries.mfa",
                         os.path.basename(split_msa_files["McrA"][0].query))
        return


if __name__ == '__main__':
    unittest.main()

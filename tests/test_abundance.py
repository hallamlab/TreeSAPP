import os
import shutil
import pytest
import unittest

from .testing_utils import get_test_data


class AbundanceTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp import treesapp_args
        from treesapp.classy import Abundance
        self.mock_abund = Abundance()
        ts_assign_output = get_test_data("test_output_TarA/")
        self.mock_abund_dir = os.path.join("tests", "abundance_test_input")
        if os.path.isdir(self.mock_abund_dir):
            shutil.rmtree(self.mock_abund_dir)
        shutil.copytree(ts_assign_output, self.mock_abund_dir)
        mock_cli_args = ["--treesapp_output", self.mock_abund_dir,
                         "--reads", get_test_data("test_TarA.1.fq"),
                         "--reverse", get_test_data("test_TarA.2.fq"),
                         "--overwrite"]

        self.abundance_parser = treesapp_args.TreeSAPPArgumentParser()
        treesapp_args.add_abundance_arguments(self.abundance_parser)

        self.mock_args = self.abundance_parser.parse_args(mock_cli_args)
        self.mock_abund.furnish_with_arguments(self.mock_args)

        return

    def tearDown(self) -> None:
        if os.path.isdir(self.mock_abund_dir):
            shutil.rmtree(self.mock_abund_dir)
        return

    def test_generate_simplebar(self):
        from treesapp.abundance import generate_simplebar
        from treesapp.phylo_seq import PQuery, PhyloPlace
        test_file = os.path.join(self.mock_abund_dir, "test_McrA_abundance_simplebar.txt")
        refpkg_name = "McrA"
        mock_node_map = {28: ['198_McrA'],
                         170: ['123_McrA', '114_McrA', '118_McrA', '121_McrA', '116_McrA'],
                         394: ['36_McrA']}
        leaves = sum(mock_node_map.values(), [])
        mock_placement_edges = [28, 170, 394]
        mock_pqueries = []
        # Instantiate the PQuery instances
        for name in ["seq1", "seq2", "seq3"]:
            pq = PQuery()
            pq.seq_name = name
            pq.ref_name = refpkg_name
            pq.consensus_placement = PhyloPlace()
            pq.consensus_placement.edge_num = mock_placement_edges.pop()
            pq.node_map = mock_node_map
            pq.abundance = 1000.0
            mock_pqueries.append(pq)

        generate_simplebar(target_marker=refpkg_name,
                           itol_bar_file=test_file,
                           pqueries=mock_pqueries)

        # Evaluate the output file
        self.assertTrue(os.path.isfile(test_file))
        abundance_sum = 0.0
        with open(test_file) as simplebar:
            for line in simplebar:
                try:
                    leaf, tpm = line.split(',')
                except ValueError:
                    continue
                if leaf in leaves:
                    abundance_sum += float(tpm)
        self.assertEqual(3000.0, abundance_sum)
        return

    def test_strip_file_to_sample_name(self):
        from treesapp.classy import Abundance
        mock_abund = Abundance()
        mock_abund.strip_file_to_sample_name("path/to/test_TarA.1.fastq")
        self.assertEqual("test_TarA", mock_abund.sample_prefix)

        mock_abund.strip_file_to_sample_name("path/to/test_sample1.1.fq")
        self.assertEqual("test_sample1", mock_abund.sample_prefix)

        mock_abund.strip_file_to_sample_name("path/to/test_sample_R1.fq")
        self.assertEqual("test_sample", mock_abund.sample_prefix)

        mock_abund.strip_file_to_sample_name("SRR7188253_1.fastq.gz")
        self.assertEqual("SRR7188253", mock_abund.sample_prefix)
        return

    def test_parse_args(self):
        from treesapp import treesapp_args
        mock_cli_args = ["--treesapp_output", self.mock_abund_dir,
                         "--reads", get_test_data("test_TarA.1.fq")]

        parser = treesapp_args.TreeSAPPArgumentParser()
        treesapp_args.add_abundance_arguments(parser)
        args = parser.parse_args(mock_cli_args)
        self.assertEqual("pe", args.pairing)
        self.assertEqual("tpm", args.metric)
        return

    def test_check_arguments(self):
        from treesapp import treesapp_args
        # Set up the parser
        parser = treesapp_args.TreeSAPPArgumentParser()
        treesapp_args.add_abundance_arguments(parser)

        # Test fail if path to reads doesn't exist
        bad_fwd_args = parser.parse_args(["--treesapp_output", self.mock_abund_dir, "--reads", "test_TarA.1.fq"])
        bad_rev_args = parser.parse_args(["--treesapp_output", self.mock_abund_dir,
                                          "--reads", get_test_data("test_TarA.1.fq"),
                                          "--reverse", "path/to/bad/rev.fq"])
        with pytest.raises(SystemExit):
            self.mock_abund.check_arguments(bad_fwd_args)
            self.mock_abund.check_arguments(bad_rev_args)

        self.mock_abund.check_arguments(self.mock_args)
        self.assertEqual("tests/abundance_test_input/intermediates/orf-call/PitchLake_TarA_Mcr_Dsr_ORFs.fna",
                         os.path.relpath(start=os.getcwd(), path=self.mock_abund.ref_nuc_seqs))
        self.assertEqual("tests/abundance_test_input/final_outputs/" + self.mock_abund.classification_tbl_name,
                         os.path.relpath(start=os.getcwd(), path=self.mock_abund.classifications))

        return

    def test_decide_stage(self):
        from treesapp.classy import Abundance
        mock_abund = Abundance()
        mock_abund.furnish_with_arguments(self.mock_args)

        # Test with report=='update'
        mock_abund.check_arguments(self.mock_args)
        mock_abund.decide_stage(self.mock_args)
        self.assertTrue(mock_abund.stage_status("align_map"))
        self.assertTrue(mock_abund.stage_status("summarise"))
        self.assertTrue(mock_abund.stage_status("sam_sum"))
        open(os.path.join(self.mock_abund_dir, "intermediates", "align_map", mock_abund.sample_prefix) + ".sam",
             'w').close()

        # Simulate a re-run when all the FASTQ files in args.reads and args.reverse are the same and SAM exists
        rerun_args = self.abundance_parser.parse_args(["--treesapp_output", self.mock_abund_dir,
                                                       "--reads", get_test_data("test_TarA.1.fq"),
                                                       "--reverse", get_test_data("test_TarA.2.fq"),
                                                       "--report", "update"])
        mock_abund.check_arguments(rerun_args)
        mock_abund.decide_stage(rerun_args)
        self.assertFalse(mock_abund.stage_status("align_map"))
        self.assertTrue(mock_abund.stage_status("sam_sum"))
        self.assertTrue(mock_abund.stage_status("summarise"))

        # Test when report=='append'
        rerun_args = self.abundance_parser.parse_args(["--treesapp_output", self.mock_abund_dir,
                                                       "--reads", get_test_data("SRR3669912_1.fastq"),
                                                       "--reverse", get_test_data("SRR3669912_2.fastq"),
                                                       "--report", "append"])
        mock_abund.check_arguments(rerun_args)
        mock_abund.decide_stage(rerun_args)
        self.assertTrue(mock_abund.stage_status("align_map"))
        self.assertTrue(mock_abund.stage_status("sam_sum"))
        self.assertTrue(mock_abund.stage_status("summarise"))

        return


if __name__ == '__main__':
    unittest.main()

import pytest
import unittest


class TreesappTester(unittest.TestCase):
    def test_assign_prot(self):
        ref_pkgs = ["M0701", "M0702", "S0001"]
        from treesapp import commands
        from treesapp import file_parsers
        from . import testing_utils as utils
        assign_commands_list = ["--fastx_input", utils.get_test_data("marker_test_suite.faa"),
                                "--targets", ','.join(ref_pkgs),
                                "--num_procs", str(2),
                                "-m", "prot",
                                "--output", "./TreeSAPP_assign/",
                                "--stringency", "relaxed",
                                "--trim_align", "--no_svm"]
        commands.assign(assign_commands_list)
        lines = file_parsers.read_marker_classification_table("./TreeSAPP_assign/final_outputs/marker_contig_map.tsv")
        self.assertEqual(15, len(lines))

    def test_assign_dna(self):
        ref_pkgs = ["M0701", "M0702"]
        from treesapp import commands
        from treesapp import file_parsers
        from . import testing_utils as utils
        assign_commands_list = ["--fastx_input", utils.get_test_data("marker_test_suite.fna"),
                                "--targets", ','.join(ref_pkgs),
                                "--num_procs", str(2),
                                "-m", "dna",
                                "--output", "./TreeSAPP_assign/",
                                "--stringency", "strict",
                                "--trim_align", "--no_svm", "--overwrite"]
        commands.assign(assign_commands_list)
        lines = file_parsers.read_marker_classification_table("./TreeSAPP_assign/final_outputs/marker_contig_map.tsv")
        self.assertEqual(6, len(lines))

    def test_create(self):
        from commands import create
        from . import testing_utils as utils
        create_commands_list = ["--fastx_input", utils.get_test_data("create_test.faa"),
                                "--refpkg_name", "Crt",
                                "--identity", "0.95",
                                "--bootstraps", str(0),
                                "--molecule", "prot",
                                "--trim_align", "--cluster", "--fast", "--headless",
                                "--screen", "Bacteria,Archaea", "--filter", "Archaea",
                                "--min_taxonomic_rank", 'p',
                                "--output", "./TreeSAPP_create",
                                "--num_proc", str(2)]
        create(create_commands_list)
        return

    def test_evaluate(self):
        from commands import evaluate
        from . import testing_utils as utils
        test_ranks = ["Class", "Species"]
        evaluate_command_list = ["--fastx_input", utils.get_test_data("McrA_eval.faa"),
                                 "--reference_marker", "McrA",
                                 "-o", "./TreeSAPP_evaluate",
                                 "--trim_align",
                                 "-m", "prot",
                                 "-t", ','.join(test_ranks),
                                 "-n", str(2)]
        evaluate(evaluate_command_list)
        return

    # TODO: Write tests for abundance, purity, layer, update and train


if __name__ == '__main__':
    unittest.main()

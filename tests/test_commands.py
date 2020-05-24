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
                                "--trim_align", "--no_svm", "--overwrite", "--delete"]
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
                                "--trim_align", "--no_svm", "--overwrite", "--delete"]
        commands.assign(assign_commands_list)
        lines = file_parsers.read_marker_classification_table("./TreeSAPP_assign/final_outputs/marker_contig_map.tsv")
        self.assertEqual(6, len(lines))

    def test_create(self):
        from commands import create
        from refpkg import ReferencePackage
        from . import testing_utils as utils
        create_commands_list = ["--fastx_input", utils.get_test_data("create_test.faa"),
                                "--accession2taxid", utils.get_test_data("create_test.accession2taxid"),
                                "--refpkg_name", "Crt",
                                "--identity", "0.95",
                                "--bootstraps", str(0),
                                "--molecule", "prot",
                                "--screen", "Bacteria,Archaea", "--filter", "Archaea",
                                "--min_taxonomic_rank", 'p',
                                "--output", "./TreeSAPP_create",
                                "--num_proc", str(2),
                                "--trim_align", "--cluster", "--fast", "--headless", "--overwrite", "--delete"]
        create(create_commands_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_create/final_outputs/Crt_build.json"
        test_refpkg.slurp()
        test_refpkg.validate()
        self.assertEqual(68, test_refpkg.num_seqs)
        return

    def test_evaluate(self):
        from commands import evaluate
        from . import testing_utils as utils
        import os
        test_ranks = ["Class"]
        evaluate_command_list = ["--fastx_input", utils.get_test_data("McrA_eval.faa"),
                                 "--refpkg_path", utils.get_treesapp_file(os.path.join("treesapp", "data",
                                                                                       "McrA_build.json")),
                                 "--accession2lin", utils.get_test_data("McrA_eval_accession_id_lineage_map.tsv"),
                                 "-o", "./TreeSAPP_evaluate",
                                 "-m", "prot",
                                 "-t", ','.join(test_ranks),
                                 "-n", str(2),
                                 "--trim_align", "--overwrite"]
        evaluate(evaluate_command_list)
        return

    def test_abundance(self):
        import os
        from commands import abundance
        from file_parsers import read_marker_classification_table
        from . import testing_utils as utils
        classification_table = os.path.join(utils.get_example_output(), "final_outputs", "marker_contig_map.tsv")
        pre_lines = read_marker_classification_table(utils.get_test_data(classification_table))
        abundance_command_list = ["--treesapp_output", utils.get_test_data("test_output_TarA/"),
                                  "--reads", utils.get_test_data("test_TarA.1.fq"),
                                  "--reverse", utils.get_test_data("test_TarA.2.fq"),
                                  "--pairing", "pe",
                                  "--num_procs", str(2),
                                  "--delete"]
        abundance(abundance_command_list)
        post_lines = read_marker_classification_table(utils.get_test_data(classification_table))
        self.assertEqual(len(pre_lines), len(post_lines))
        return

    def test_purity(self):
        from commands import purity
        import os
        from . import testing_utils as utils
        purity_command_list = ["--fastx_input", utils.get_treesapp_file("dev_utils/TIGRFAM_seed_named.faa"),
                               "--extra_info", utils.get_treesapp_file("dev_utils/TIGRFAM_info.tsv"),
                               "--output", "./TreeSAPP_purity",
                               "--refpkg_path", os.path.join(utils.get_treesapp_path(),
                                                             "treesapp", "data", "McrA_build.json"),
                               "--trim_align", "--molecule", "prot", "-n", str(2)]
        purity(purity_command_list)
        return

    def test_layer(self):
        import os
        from commands import layer
        from . import testing_utils as utils
        from file_parsers import read_marker_classification_table
        classification_table = os.path.join(utils.get_example_output(), "final_outputs", "marker_contig_map.tsv")
        pre_lines = read_marker_classification_table(utils.get_test_data(classification_table))
        layer_command_list = ["--treesapp_output", utils.get_test_data("test_output_TarA/")]
        layer(layer_command_list)
        post_lines = read_marker_classification_table(utils.get_test_data(classification_table))
        self.assertEqual(len(pre_lines), len(post_lines))
        return

    def test_update(self):
        from commands import update
        from . import testing_utils as utils
        update_command_list = ["--fastx_input", utils.get_test_data("McrA_eval.faa"),
                               "--refpkg_name", "McrA",
                               "--treesapp_output", utils.get_test_data("test_output_TarA/"),
                               "--identity", str(0.97),
                               "--output", "./TreeSAPP_update",
                               "--num_proc", str(2),
                               "--molecule", "prot",
                               "-b", str(0),
                               "--trim_align", "--cluster", "--fast", "--headless", "--overwrite", "--delete"]
        update(update_command_list)
        return

    def test_train(self):
        import os
        from commands import train
        from . import testing_utils as utils
        train_command_list = ["--fastx_input", utils.get_test_data("McrA_eval.faa"),
                              "--output", "./TreeSAPP_train",
                              "--refpkg_path", "McrA",
                              "--pkg_path", os.path.join(utils.get_treesapp_path(), "treesapp", "data"),
                              "--num_proc", str(2),
                              "--molecule", "prot",
                              "--trim_align", "--delete", "--overwrite"]
        train(train_command_list)
        return


if __name__ == '__main__':
    unittest.main()

    # Clean up output directories
    import os
    output_prefix = "TreeSAPP_"
    for test_name in ["assign", "train", "update", "evaluate", "create"]:
        output_dir = output_prefix + test_name
        if os.path.isdir(output_dir):
            os.rmdir(output_dir)
import unittest
import os


class TreesappTester(unittest.TestCase):
    def test_assign_prot(self):
        ref_pkgs = ["McrA", "M0702", "S0001"]
        from treesapp.commands import assign
        from treesapp.file_parsers import read_classification_table
        from .testing_utils import get_test_data
        assign_commands_list = ["--fastx_input", get_test_data("marker_test_suite.faa"),
                                "--targets", ','.join(ref_pkgs),
                                "--num_procs", str(2),
                                "-m", "prot",
                                "--output", "./TreeSAPP_assign/",
                                "--stringency", "relaxed",
                                "--trim_align", "--overwrite", "--delete"]
        assign(assign_commands_list)
        lines = read_classification_table("./TreeSAPP_assign/final_outputs/marker_contig_map.tsv")
        self.assertEqual(15, len(lines))

    def test_assign_dna(self):
        ref_pkgs = ["M0701", "M0702"]
        from treesapp.commands import assign
        from treesapp.file_parsers import read_classification_table
        from .testing_utils import get_test_data
        assign_commands_list = ["--fastx_input", get_test_data("marker_test_suite.fna"),
                                "--targets", ','.join(ref_pkgs),
                                "--num_procs", str(2),
                                "-m", "dna",
                                "--output", "./TreeSAPP_assign/",
                                "--stringency", "strict",
                                "--trim_align", "--overwrite", "--delete"]
        assign(assign_commands_list)
        lines = read_classification_table("./TreeSAPP_assign/final_outputs/marker_contig_map.tsv")
        self.assertEqual(6, len(lines))

    def test_create(self):
        from treesapp.commands import create
        from treesapp.refpkg import ReferencePackage
        from .testing_utils import get_test_data
        create_commands_list = ["--fastx_input", get_test_data("create_test.faa"),
                                "--accession2taxid", get_test_data("create_test.accession2taxid"),
                                "--refpkg_name", "Crt",
                                "--similarity", "0.95",
                                "--bootstraps", str(0),
                                "--molecule", "prot",
                                "--screen", "Bacteria,Archaea", "--filter", "Archaea",
                                "--min_taxonomic_rank", 'p',
                                "--output", "./TreeSAPP_create",
                                "--num_proc", str(2),
                                "--trim_align", "--cluster", "--fast", "--headless", "--overwrite", "--delete"]
        create(create_commands_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_create/final_outputs/Crt_build.pkl"
        test_refpkg.slurp()
        test_refpkg.validate()
        self.assertEqual(68, test_refpkg.num_seqs)
        return

    def test_create_eggnog(self):
        from commands import create
        from refpkg import ReferencePackage
        from .testing_utils import get_test_data
        create_commands_list = ["--fastx_input", get_test_data("ENOG4111FIN.txt"),
                                "--output", "./TreeSAPP_create_PuhA",
                                "--refpkg_name", "PuhA",
                                "--similarity", "0.97",
                                "--bootstraps", str(0),
                                "--molecule", "prot",
                                "--screen", "Bacteria,Archaea",
                                "--num_proc", str(2),
                                "--min_taxonomic_rank", 'p',
                                "--trim_align", "--outdet_align", "--cluster", "--fast", "--headless",
                                "--overwrite", "--delete"]
        create(create_commands_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_create_PuhA/final_outputs/PuhA_build.pkl"
        test_refpkg.slurp()
        test_refpkg.validate()
        self.assertEqual(44, test_refpkg.num_seqs)
        return

    def test_evaluate(self):
        from commands import evaluate
        from .testing_utils import get_test_data
        evaluate_command_list = ["--fastx_input", get_test_data("McrA_eval.faa"),
                                 "--refpkg_path", get_test_data(os.path.join("refpkgs", "McrA_build.pkl")),
                                 "--accession2lin", get_test_data("McrA_eval_accession_id_lineage_map.tsv"),
                                 "-o", "./TreeSAPP_evaluate",
                                 "-m", "prot",
                                 "--taxonomic_ranks", "class", "order",
                                 "-n", str(2),
                                 "--trim_align", "--overwrite", "--delete"]
        evaluate(evaluate_command_list)
        return

    def test_abundance(self):
        from commands import abundance
        from file_parsers import read_classification_table
        from .testing_utils import get_test_data, get_example_output
        classification_table = os.path.join(get_example_output(), "final_outputs", "marker_contig_map.tsv")
        pre_lines = read_classification_table(get_test_data(classification_table))
        abundance_command_list = ["--treesapp_output", get_test_data("test_output_TarA/"),
                                  "--reads", get_test_data("test_TarA.1.fq"),
                                  "--reverse", get_test_data("test_TarA.2.fq"),
                                  "--pairing", "pe",
                                  "--num_procs", str(2),
                                  "--delete"]
        abundance(abundance_command_list)
        post_lines = read_classification_table(get_test_data(classification_table))
        self.assertEqual(len(pre_lines), len(post_lines))
        return

    def test_purity(self):
        from commands import purity
        from .testing_utils import get_treesapp_file, get_test_data
        purity_command_list = ["--fastx_input", get_treesapp_file("dev_utils/TIGRFAM_seed_named.faa"),
                               "--extra_info", get_treesapp_file("dev_utils/TIGRFAM_info.tsv"),
                               "--output", "./TreeSAPP_purity",
                               "--refpkg_path", get_test_data(os.path.join("refpkgs", "McrA_build.pkl")),
                               "--trim_align", "--molecule", "prot", "-n", str(2)]
        purity(purity_command_list)
        return

    def test_layer(self):
        from commands import layer
        from .testing_utils import get_test_data, get_example_output
        from file_parsers import read_classification_table
        classification_table = os.path.join(get_example_output(), "final_outputs", "marker_contig_map.tsv")
        pre_lines = read_classification_table(get_test_data(classification_table))
        layer_command_list = ["--treesapp_output", get_test_data("test_output_TarA/")]
        layer_command_list += ["--refpkg_dir", get_test_data("refpkgs/")]
        layer(layer_command_list)
        post_lines = read_classification_table(get_test_data(classification_table))
        self.assertEqual(len(pre_lines), len(post_lines))
        return

    def test_update(self):
        from commands import update
        from refpkg import ReferencePackage
        from .testing_utils import get_test_data
        update_command_list = ["--fastx_input", get_test_data("Photosynthesis/PuhA/ENOG4111FIN_PF03967_seed.faa"),
                               "--refpkg_path", get_test_data(os.path.join("refpkgs", "PuhA_build.json")),
                               "--treesapp_output", get_test_data("assign_PF03967_seed/"),
                               "--similarity", str(0.97),
                               "--output", "./TreeSAPP_update",
                               "--num_proc", str(2),
                               "--molecule", "prot",
                               "-b", str(0),
                               "--trim_align", "--cluster", "--fast", "--headless",
                               "--overwrite", "--delete", "--skip_assign", "--resolve"]
        update(update_command_list)
        # Test 2 with a seq2lineage table
        update_command_list += ["--seqs2lineage", get_test_data("PF03967_seed_seqs2lineage.txt")]
        update(update_command_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_update/final_outputs/PuhA_build.json"
        test_refpkg.slurp()
        test_refpkg.validate()
        self.assertEqual(46, test_refpkg.num_seqs)
        return

    def test_train(self):
        from commands import train
        from .testing_utils import get_test_data
        train_command_list = ["--fastx_input",  get_test_data("ENOG4111FIN.txt"),
                              "--output", "./TreeSAPP_train",
                              "--refpkg_path", get_test_data(os.path.join("refpkgs", "PuhA_build.pkl")),
                              "--accession2lin", get_test_data("ENOG4111FIN_accession_id_lineage_map.tsv"),
                              "--num_proc", str(4),
                              "--molecule", "prot",
                              "--svm_kernel", "rbf",
                              "--trim_align", "--delete", "--overwrite"]
        train(train_command_list)
        return

    def test_classifier_trainer(self):
        from classifier_trainer import classifier_trainer as ct_main
        from .testing_utils import get_test_data
        command_list = ["--fastx_input", get_test_data("McrA_eval.faa"),
                        "--output", "./TreeSAPP_svc_train",
                        "--refpkg_path", get_test_data(os.path.join("refpkgs", "McrA_build.pkl")),
                        "--accession2lin", get_test_data("McrA_eval_accession_id_lineage_map.tsv"),
                        "--num_proc", str(2),
                        "--molecule", "prot",
                        "--trim_align", "--delete", "--overwrite"]
        ct_main(command_list)

    def test_package(self):
        from commands import package
        from .testing_utils import get_test_data
        command_list = ["view",
                        "prefix",
                        "--refpkg_path", get_test_data(os.path.join("refpkgs", "McrA_build.pkl")),
                        "--output", "./TreeSAPP_package"]
        package(command_list)
        return


if __name__ == '__main__':
    unittest.main()

    # Clean up output directories
    output_prefix = "TreeSAPP_"
    for test_name in ["assign", "train", "update", "evaluate", "create"]:
        output_dir = output_prefix + test_name
        if os.path.isdir(output_dir):
            os.rmdir(output_dir)

import unittest
import pytest
import os
from shutil import rmtree


class TreesappTester(unittest.TestCase):
    def setUp(self) -> None:
        from .testing_utils import get_test_data
        self.num_procs = 4
        # FASTA files
        self.aa_test_fa = get_test_data("marker_test_suite.faa")
        self.nt_test_fa = get_test_data("marker_test_suite.fna")

        # Reference package pickles
        self.refpkg_dir = get_test_data(os.path.join("refpkgs"))
        self.mcra_pkl = get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.mcrb_pkl = get_test_data(os.path.join("refpkgs", "McrB_build.pkl"))
        self.puha_pkl = get_test_data(os.path.join("refpkgs", "PuhA_build.pkl"))

        # Output examples
        self.ts_assign_output = get_test_data("test_output_TarA/")

        return

    def tearDown(self) -> None:
        # Clean up output directories
        output_prefix = os.path.join(os.path.abspath("./"), "TreeSAPP_")
        for test_name in ["assign", "train", "update", "evaluate", "create", "package", "purity", "MCC", "colour"]:
            output_dir = output_prefix + test_name
            if os.path.isdir(output_dir):
                rmtree(output_dir)
        return

    def test_abundance(self):
        from treesapp.commands import abundance
        from treesapp.file_parsers import read_classification_table
        from .testing_utils import get_test_data
        classification_table = os.path.join(self.ts_assign_output, "final_outputs", "marker_contig_map.tsv")
        pre_lines = read_classification_table(get_test_data(classification_table))
        abundance_command_list = ["--treesapp_output", self.ts_assign_output,
                                  "--reads", get_test_data("test_TarA.1.fq"),
                                  "--reverse", get_test_data("test_TarA.2.fq"),
                                  "--pairing", "pe",
                                  "--num_procs", str(self.num_procs),
                                  "--delete"]
        abundance(abundance_command_list)
        post_lines = read_classification_table(get_test_data(classification_table))
        self.assertEqual(len(pre_lines), len(post_lines))
        return

    def test_assign_prot(self):
        ref_pkgs = ["McrA", "M0702", "S0001"]
        from treesapp.commands import assign
        from treesapp.file_parsers import read_classification_table
        assign_commands_list = ["--fastx_input", self.aa_test_fa,
                                "--targets", ','.join(ref_pkgs),
                                "--refpkg_dir", self.refpkg_dir,
                                "--num_procs", str(self.num_procs),
                                "-m", "prot",
                                "--output", "./TreeSAPP_assign/",
                                "--stringency", "relaxed",
                                "--placement_summary", "max_lwr",
                                "--trim_align", "--overwrite", "--delete", "--svm"]
        assign(assign_commands_list)
        lines = read_classification_table("./TreeSAPP_assign/final_outputs/marker_contig_map.tsv")
        self.assertEqual(15, len(lines))
        return

    def test_assign_dna(self):
        ref_pkgs = ["M0701", "M0702"]
        from treesapp.commands import assign
        from treesapp.file_parsers import read_classification_table
        assign_commands_list = ["--fastx_input", self.nt_test_fa,
                                "--targets", ','.join(ref_pkgs),
                                "--num_procs", str(self.num_procs),
                                "-m", "dna",
                                "--output", "./TreeSAPP_assign/",
                                "--stringency", "strict",
                                "--trim_align", "--overwrite", "--delete"]
        assign(assign_commands_list)
        lines = read_classification_table("./TreeSAPP_assign/final_outputs/marker_contig_map.tsv")
        self.assertEqual(6, len(lines))

    def test_colour(self):
        from treesapp.commands import colour
        colour_commands = ["-r", self.mcra_pkl, self.mcrb_pkl,
                           "-l", "family",
                           "-o", "./TreeSAPP_colour",
                           "--filter", "Methanococcales",
                           "--min_proportion", str(0.01),
                           "--no_polyphyletic"]
        colour(colour_commands)
        self.assertEqual(True, True)

    def test_colour_phenotypes(self):
        from treesapp.commands import colour
        from .testing_utils import get_test_data
        colour_commands = ["-r", self.mcra_pkl,
                           "-o", "./TreeSAPP_colour",
                           "--palette", "viridis",
                           "--taxa_map", get_test_data("Mcr_taxonomy-phenotype_map.tsv")]
        colour(colour_commands)
        self.assertTrue(True)

    def test_create(self):
        from treesapp.commands import create
        from treesapp.refpkg import ReferencePackage
        from .testing_utils import get_test_data
        create_commands_list = ["--fastx_input", get_test_data("create_test.faa"),
                                "--accession2taxid", get_test_data("create_test.accession2taxid"),
                                "--seqs2lineage", get_test_data("create_test_accession2lin.tsv"),
                                "--refpkg_name", "Crt",
                                "--similarity", "0.95",
                                "--bootstraps", str(0),
                                "--molecule", "prot",
                                "--screen", "Bacteria", "--filter", "Archaea",
                                "--min_taxonomic_rank", 'p',
                                "--output", "./TreeSAPP_create",
                                "--num_procs", str(self.num_procs),
                                "--trim_align", "--cluster", "--fast", "--headless", "--overwrite", "--delete"]
        create(create_commands_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_create/final_outputs/Crt_build.pkl"
        test_refpkg.slurp()
        self.assertTrue(test_refpkg.validate())
        self.assertEqual(69, test_refpkg.num_seqs)
        self.assertEqual(2, len(test_refpkg.pfit))
        self.assertTrue(test_refpkg.pfit[0] < 0)
        return

    def test_create_eggnog(self):
        from treesapp.commands import create
        from treesapp.refpkg import ReferencePackage
        from .testing_utils import get_test_data
        create_commands_list = ["--fastx_input", get_test_data("ENOG4111FIN.txt"),
                                "--output", "./TreeSAPP_create",
                                "--refpkg_name", "PuhA",
                                "--similarity", "0.90",
                                "--profile", get_test_data("PuhA_search.hmm"),
                                "--molecule", "prot",
                                "--screen", "Bacteria,Archaea",
                                "--num_procs", str(self.num_procs),
                                "--raxml_model", "LG+F+R4",
                                "--min_taxonomic_rank", 'p',
                                "--stage", "support",
                                "--trim_align", "--outdet_align", "--cluster", "--headless",
                                "--overwrite", "--delete"]
        create(create_commands_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_create/final_outputs/PuhA_build.pkl"
        test_refpkg.slurp()
        self.assertTrue(test_refpkg.validate())
        self.assertEqual(39, test_refpkg.num_seqs)
        self.assertEqual("LG+F+R4", test_refpkg.sub_model)
        return

    def test_create_accession2lin(self):
        from treesapp.commands import create
        from treesapp.refpkg import ReferencePackage
        from .testing_utils import get_test_data
        cmd_list = ["--fastx_input", get_test_data("PF01280_test.fasta"),
                    "--accession2lin", get_test_data("PF01280_test.tsv"),
                    "--output", "./TreeSAPP_create",
                    "--refpkg_name", "PF01280",
                    "--similarity", "0.95",
                    "--bootstraps", str(0),
                    "--molecule", "prot",
                    "--screen", "Archaea", "--filter", "Bacteria,Eukaryota",
                    "--min_taxonomic_rank", 'g',
                    "--num_procs", str(self.num_procs),
                    "--stage", "evaluate",
                    "--trim_align", "--cluster", "--fast", "--headless", "--overwrite", "--delete"]
        create(cmd_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_create/final_outputs/PF01280_build.pkl"
        test_refpkg.slurp()
        self.assertTrue(test_refpkg.validate())
        self.assertEqual(51, test_refpkg.num_seqs)
        return

    def test_evaluate(self):
        from treesapp.commands import evaluate
        from .testing_utils import get_test_data
        evaluate_command_list = ["--fastx_input", get_test_data("McrA_eval.faa"),
                                 "--refpkg_path", self.mcra_pkl,
                                 "--accession2lin", get_test_data("McrA_eval_accession_id_lineage_map.tsv"),
                                 "-o", "./TreeSAPP_evaluate",
                                 "-m", "prot",
                                 "--taxonomic_ranks", "class", "order",
                                 "-n", str(self.num_procs),
                                 "--trim_align", "--overwrite", "--delete"]
        evaluate(evaluate_command_list)
        self.assertEqual(True, True)
        return

    def test_info(self):
        from treesapp.commands import info
        with pytest.raises(SystemExit):
            info(["--verbose", "--refpkg_dir", "./"])
        info([])
        return

    def test_layer(self):
        from treesapp.commands import layer
        from .testing_utils import get_test_data
        from treesapp.file_parsers import read_classification_table
        # Layering annotations from multiple reference packages
        original_table = os.path.join(self.ts_assign_output, "final_outputs", "marker_contig_map.tsv")
        layered_table = os.path.join(self.ts_assign_output, "final_outputs",
                                     "extra_annotated_marker_contig_map.tsv")
        pre_lines = read_classification_table(get_test_data(original_table))
        layer_command_list = ["--treesapp_output", self.ts_assign_output]
        layer_command_list += ["--refpkg_dir", self.refpkg_dir]
        layer(layer_command_list)
        post_lines = read_classification_table(get_test_data(layered_table))
        self.assertEqual(len(pre_lines), len(post_lines))

        # With a different reference package, XmoA, just to be sure
        classification_table = os.path.join(get_test_data("p_amoA_FunGene9.5_isolates_assign/"),
                                            "final_outputs", "marker_contig_map.tsv")
        pre_lines = read_classification_table(get_test_data(classification_table))
        layer_command_list = ["--colours_style", get_test_data("XmoA_Function.txt"),
                              "--treesapp_output", get_test_data("p_amoA_FunGene9.5_isolates_assign/"),
                              "--refpkg_dir", self.refpkg_dir]
        layer(layer_command_list)
        post_lines = read_classification_table(get_test_data(classification_table))
        self.assertEqual(len(pre_lines), len(post_lines))
        return

    def test_purity(self):
        from treesapp.commands import purity
        from .testing_utils import get_treesapp_file
        purity_command_list = ["--fastx_input", get_treesapp_file("dev_utils/TIGRFAM_seed_named.faa"),
                               "--extra_info", get_treesapp_file("dev_utils/TIGRFAM_info.tsv"),
                               "--output", "./TreeSAPP_purity",
                               "--refpkg_path", self.mcra_pkl,
                               "--num_procs", str(self.num_procs),
                               "--trim_align", "--molecule", "prot"]
        purity(purity_command_list)
        self.assertEqual(True, True)
        return

    def test_package(self):
        from treesapp.commands import package
        from treesapp.refpkg import ReferencePackage
        view_command_list = ["view", "lineage_ids",
                             "--refpkg_path", self.mcra_pkl,
                             "--output", "./TreeSAPP_package"]
        package(view_command_list)
        edit_command_list = ["edit",
                             "f__msa", self.aa_test_fa,
                             "--refpkg_path", self.mcra_pkl,
                             "--output", "./TreeSAPP_package"]
        package(edit_command_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_package/McrA_build.pkl"
        test_refpkg.slurp()
        self.assertFalse(test_refpkg.validate())
        return

    def test_mcc_calculator(self):
        from treesapp import MCC_calculator
        from .testing_utils import get_test_data
        cmd = ["--fastx_input", get_test_data("EggNOG_McrA.faa"),
               "--annot_map", get_test_data("EggNOG_McrA_annot_map.tsv"),
               "--output", "./TreeSAPP_MCC",
               "--refpkg_dir", self.refpkg_dir,
               "--targets", "McrA",
               "--molecule", "prot",
               "--tool", "treesapp",
               "--num_procs", str(self.num_procs),
               "--delete", "--svm", "--overwrite"]
        MCC_calculator.mcc_calculator(cmd)
        self.assertEqual(True, True)
        return

    def test_train(self):
        import csv
        from treesapp.commands import train
        from .testing_utils import get_test_data
        output_dir_path = "./TreeSAPP_train"
        max_ex = 50
        train_command_list = ["--fastx_input", get_test_data("ENOG4111FIN.txt"),
                              "--annot_map", get_test_data("ENOG4111FIN_annot_map.tsv"),
                              "--output", output_dir_path,
                              "--refpkg_path", self.puha_pkl,
                              "--accession2lin", get_test_data("ENOG4111FIN_accession_id_lineage_map.tsv"),
                              "--max_examples", str(max_ex),
                              "--num_proc", str(self.num_procs),
                              "--molecule", "prot",
                              "--svm_kernel", "rbf",
                              "--classifier", "bin",
                              "--trim_align", "--delete", "--overwrite"]
        train(train_command_list)
        rank_list = []
        with open(os.path.join(output_dir_path, "final_outputs", "placement_info.tsv")) as placement_tbl:
            csv_handler = csv.reader(placement_tbl, delimiter="\t")
            next(csv_handler)
            for fields in csv_handler:
                rank_list.append(fields[1])

        self.assertEqual(0, len({"class", "order", "family", "genus", "species"}.difference(set(rank_list))))
        self.assertEqual(max_ex, len(rank_list))
        return

    def test_update_resolve(self):
        from treesapp.commands import update
        from treesapp.refpkg import ReferencePackage
        from .testing_utils import get_test_data
        update_command_list = ["--fastx_input", get_test_data("Photosynthesis/PuhA/ENOG4111FIN_PF03967_seed.faa"),
                               "--refpkg_path", self.puha_pkl,
                               "--treesapp_output", get_test_data("assign_SwissProt_PuhA/"),
                               "--output", "./TreeSAPP_update",
                               "--num_proc", str(self.num_procs),
                               "--molecule", "prot",
                               "--trim_align", "--cluster", "--fast", "--headless",
                               "--overwrite", "--delete", "--skip_assign", "--resolve"]
        update(update_command_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_update/final_outputs/PuhA_build.pkl"
        test_refpkg.slurp()
        self.assertTrue(test_refpkg.validate())
        self.assertEqual(49, test_refpkg.num_seqs)
        return

    def test_update_seqs2lineage(self):
        from treesapp.commands import update
        from treesapp.refpkg import ReferencePackage
        from .testing_utils import get_test_data
        update_command_list = ["--fastx_input", get_test_data("Photosynthesis/PuhA/ENOG4111FIN_PF03967_seed.faa"),
                               "--refpkg_path", self.puha_pkl,
                               "--treesapp_output", get_test_data("assign_SwissProt_PuhA/"),
                               "--seqs2lineage", get_test_data("SwissProt_PuhA_seqs2lineage.txt"),
                               "--output", "./TreeSAPP_update",
                               "--num_procs", str(self.num_procs),
                               "--molecule", "prot",
                               "-b", str(0),
                               "--trim_align", "--cluster", "--fast", "--headless",
                               "--overwrite", "--delete", "--skip_assign"]
        update(update_command_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_update/final_outputs/PuhA_build.pkl"
        test_refpkg.slurp()
        self.assertTrue(test_refpkg.validate())
        self.assertEqual(49, test_refpkg.num_seqs)
        return

    # def test_tmp(self):
    #     from treesapp.commands import create
    #     base_dir = "/home/connor/Bioinformatics/Hallam_projects/RefPkgs/"
    #     cmd = "".format(base_dir)
    #     create(cmd.split())
    #     return


if __name__ == '__main__':
    unittest.main()

import unittest
import os
from shutil import rmtree


class TreesappTester(unittest.TestCase):
    def tearDown(self) -> None:
        # Clean up output directories
        output_prefix = os.path.join(os.path.abspath("./"), "TreeSAPP_")
        for test_name in ["assign", "train", "update", "evaluate", "create", "package", "purity", "MCC"]:
            output_dir = output_prefix + test_name
            if os.path.isdir(output_dir):
                rmtree(output_dir)
        return

    def test_assign_prot(self):
        ref_pkgs = ["McrA", "M0702", "S0001"]
        from treesapp.commands import assign
        from treesapp.file_parsers import read_classification_table
        from .testing_utils import get_test_data
        assign_commands_list = ["--fastx_input", get_test_data("marker_test_suite.faa"),
                                "--targets", ','.join(ref_pkgs),
                                "--refpkg_dir", get_test_data(os.path.join("refpkgs")),
                                "--num_procs", str(2),
                                "-m", "prot",
                                "--output", "./TreeSAPP_assign/",
                                "--stringency", "relaxed",
                                "--trim_align", "--overwrite", "--delete", "--svm"]
        assign(assign_commands_list)
        lines = read_classification_table("./TreeSAPP_assign/final_outputs/marker_contig_map.tsv")
        self.assertEqual(16, len(lines))
        return

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
                                "--screen", "Bacteria", "--filter", "Archaea",
                                "--min_taxonomic_rank", 'p',
                                "--output", "./TreeSAPP_create",
                                "--num_proc", str(2),
                                "--trim_align", "--cluster", "--fast", "--headless", "--overwrite", "--delete"]
        create(create_commands_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_create/final_outputs/Crt_build.pkl"
        test_refpkg.slurp()
        test_refpkg.validate()
        self.assertEqual(69, test_refpkg.num_seqs)
        self.assertEqual(2, len(test_refpkg.pfit))
        self.assertTrue(test_refpkg.pfit[0] < 0)
        return

    def test_create_eggnog(self):
        from treesapp.commands import create
        from refpkg import ReferencePackage
        from .testing_utils import get_test_data
        create_commands_list = ["--fastx_input", get_test_data("ENOG4111FIN.txt"),
                                "--output", "./TreeSAPP_create_PuhA",
                                "--refpkg_name", "PuhA",
                                "--similarity", "0.90",
                                "--profile", get_test_data("PuhA_search.hmm"),
                                "--molecule", "prot",
                                "--screen", "Bacteria,Archaea",
                                "--num_proc", str(2),
                                "--min_taxonomic_rank", 'p',
                                "--stage", "support",
                                "--trim_align", "--outdet_align", "--cluster", "--headless",
                                "--overwrite", "--delete"]
        create(create_commands_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_create_PuhA/final_outputs/PuhA_build.pkl"
        test_refpkg.slurp()
        test_refpkg.validate()
        self.assertEqual(39, test_refpkg.num_seqs)
        return

    def test_create_accession2lin(self):
        from treesapp.commands import create
        from refpkg import ReferencePackage
        from .testing_utils import get_test_data
        cmd_list = ["--fastx_input", get_test_data("PF01280_test.fasta"),
                    "--accession2lin", get_test_data("PF01280_test.tsv"),
                    "--output", "./TreeSAPP_create_PF01280",
                    "--refpkg_name", "PF01280",
                    "--similarity", "0.95",
                    "--bootstraps", str(0),
                    "--molecule", "prot",
                    "--screen", "Archaea", "--filter", "Bacteria,Eukaryota",
                    "--min_taxonomic_rank", 'g',
                    "--num_proc", str(2),
                    "--stage", "evaluate",
                    "--trim_align", "--cluster", "--fast", "--headless", "--overwrite", "--delete"]
        create(cmd_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_create_PF01280/final_outputs/PF01280_build.pkl"
        test_refpkg.slurp()
        test_refpkg.validate()
        self.assertEqual(52, test_refpkg.num_seqs)
        return

    def test_evaluate(self):
        from treesapp.commands import evaluate
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
        self.assertEqual(True, True)
        return

    def test_abundance(self):
        from treesapp.commands import abundance
        from file_parsers import read_classification_table
        from .testing_utils import get_test_data
        classification_table = os.path.join(get_test_data("test_output_TarA/"),
                                            "final_outputs", "marker_contig_map.tsv")
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
        from treesapp.commands import purity
        from .testing_utils import get_treesapp_file, get_test_data
        purity_command_list = ["--fastx_input", get_treesapp_file("dev_utils/TIGRFAM_seed_named.faa"),
                               "--extra_info", get_treesapp_file("dev_utils/TIGRFAM_info.tsv"),
                               "--output", "./TreeSAPP_purity",
                               "--refpkg_path", get_test_data(os.path.join("refpkgs", "McrA_build.pkl")),
                               "--trim_align", "--molecule", "prot", "-n", str(2)]
        purity(purity_command_list)
        self.assertEqual(True, True)
        return

    def test_layer(self):
        from treesapp.commands import layer
        from .testing_utils import get_test_data
        from file_parsers import read_classification_table
        original_table = os.path.join(get_test_data("test_output_TarA/"), "final_outputs", "marker_contig_map.tsv")
        layered_table = os.path.join(get_test_data("test_output_TarA/"), "final_outputs",
                                     "extra_annotated_marker_contig_map.tsv")
        pre_lines = read_classification_table(get_test_data(original_table))
        layer_command_list = ["--treesapp_output", get_test_data("test_output_TarA/")]
        layer_command_list += ["--refpkg_dir", get_test_data("refpkgs/")]
        layer(layer_command_list)
        post_lines = read_classification_table(get_test_data(layered_table))
        self.assertEqual(len(pre_lines), len(post_lines))
        return

    def test_xmoa_layer(self):
        from treesapp.commands import layer
        from .testing_utils import get_test_data
        from file_parsers import read_classification_table
        classification_table = os.path.join(get_test_data("p_amoA_FunGene9.5_isolates_assign/"),
                                            "final_outputs", "marker_contig_map.tsv")
        pre_lines = read_classification_table(get_test_data(classification_table))
        layer_command_list = ["--colours_style", get_test_data("XmoA_Function.txt"),
                              "--treesapp_output", get_test_data("p_amoA_FunGene9.5_isolates_assign/"),
                              "--refpkg_dir", get_test_data("refpkgs/")]
        layer(layer_command_list)
        post_lines = read_classification_table(get_test_data(classification_table))
        self.assertEqual(len(pre_lines), len(post_lines))
        return

    def test_update_resolve(self):
        from treesapp.commands import update
        from refpkg import ReferencePackage
        from .testing_utils import get_test_data
        update_command_list = ["--fastx_input", get_test_data("Photosynthesis/PuhA/ENOG4111FIN_PF03967_seed.faa"),
                               "--refpkg_path", get_test_data(os.path.join("refpkgs", "PuhA_build.pkl")),
                               "--treesapp_output", get_test_data("assign_SwissProt_PuhA/"),
                               "--output", "./TreeSAPP_update",
                               "--num_proc", str(4),
                               "--molecule", "prot",
                               "-b", str(0),
                               "--trim_align", "--cluster", "--fast", "--headless",
                               "--overwrite", "--delete", "--skip_assign", "--resolve"]
        update(update_command_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_update/final_outputs/PuhA_build.pkl"
        test_refpkg.slurp()
        test_refpkg.validate()
        self.assertEqual(49, test_refpkg.num_seqs)
        return

    def test_update_seqs2lineage(self):
        from treesapp.commands import update
        from refpkg import ReferencePackage
        from .testing_utils import get_test_data
        update_command_list = ["--fastx_input", get_test_data("Photosynthesis/PuhA/ENOG4111FIN_PF03967_seed.faa"),
                               "--refpkg_path", get_test_data(os.path.join("refpkgs", "PuhA_build.pkl")),
                               "--treesapp_output", get_test_data("assign_SwissProt_PuhA/"),
                               "--seqs2lineage", get_test_data("SwissProt_PuhA_seqs2lineage.txt"),
                               "--output", "./TreeSAPP_update",
                               "--num_proc", str(2),
                               "--molecule", "prot",
                               "-b", str(0),
                               "--trim_align", "--cluster", "--fast", "--headless",
                               "--overwrite", "--delete", "--skip_assign"]
        update(update_command_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__json = "./TreeSAPP_update/final_outputs/PuhA_build.pkl"
        test_refpkg.slurp()
        test_refpkg.validate()
        self.assertEqual(49, test_refpkg.num_seqs)
        return

    def test_train(self):
        import csv
        from treesapp.commands import train
        from .testing_utils import get_test_data
        output_dir_path = "./TreeSAPP_train"
        train_command_list = ["--fastx_input", get_test_data("ENOG4111FIN.txt"),
                              "--annot_map", get_test_data("ENOG4111FIN_annot_map.tsv"),
                              "--output", output_dir_path,
                              "--refpkg_path", get_test_data(os.path.join("refpkgs", "PuhA_build.pkl")),
                              "--accession2lin", get_test_data("ENOG4111FIN_accession_id_lineage_map.tsv"),
                              "--num_proc", str(4),
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
        self.assertEqual(120, len(rank_list))
        return

    def test_package(self):
        from treesapp.commands import package
        from .testing_utils import get_test_data
        command_list = ["view",
                        "lineage_ids",
                        "--refpkg_path", get_test_data(os.path.join("refpkgs", "McrA_build.pkl")),
                        "--output", "./TreeSAPP_package"]
        package(command_list)
        self.assertEqual(True, True)
        return

    def test_mcc_calculator(self):
        from treesapp import MCC_calculator
        from .testing_utils import get_test_data
        cmd = ["--fastx_input", get_test_data("EggNOG_McrA.faa"),
               "--annot_map", get_test_data("EggNOG_McrA_annot_map.tsv"),
               "--output", "./TreeSAPP_MCC",
               "--targets", "McrA",
               "--molecule", "prot",
               "--tool", "treesapp",
               "--num_procs", str(4),
               "--delete", "--svm", "--overwrite"]
        MCC_calculator.mcc_calculator(cmd)
        self.assertEqual(True, True)
        return

    def test_tmp(self):
        from treesapp.commands import create
        base_dir = "/home/connor/Bioinformatics/Hallam_projects/RefPkgs/"
        cmd = "-n 4 -m prot --headless --overwrite --trim_align --delete --fast --cluster --screen Bacteria,Archaea --min_taxonomic_rank g" \
              " -i {0}/Translation/PF00900/PF00900_ncbi_ML60_uclust99.faa " \
              "-c PF00900 -p 0.95 " \
              "-o {0}/Translation/PF00900/seed_refpkg/ " \
              "--accession2taxid /mnt/sdb/Hallam_projects/Hallam_Databases/raw/Taxonomy/prot.accession2taxid".format(base_dir)
              # "--accession2lin {0}/Translation/PF00900/PF00900_accession2lineage_map.tsv".format(base_dir)
        create(cmd.split())
        return


if __name__ == '__main__':
    unittest.main()

import unittest
import pytest
import os
from shutil import rmtree, copyfile


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
        self.xmoa_pkl = get_test_data(os.path.join("refpkgs", "XmoA_build.pkl"))

        # Output examples
        self.ts_assign_output = get_test_data("test_output_TarA/")

        return

    def tearDown(self) -> None:
        # Clean up output directories
        output_prefix = os.path.join(os.path.abspath("./"), "TreeSAPP_")
        for test_name in ["assign", "train", "update", "evaluate", "create", "package", "purity", "colour"]:
            output_dir = output_prefix + test_name
            if os.path.isdir(output_dir):
                rmtree(output_dir)
        return

    def test_abundance(self):
        from treesapp.abundance import abundance
        from treesapp.file_parsers import read_classification_table
        from treesapp.phylo_seq import assignments_to_pqueries
        from .testing_utils import get_test_data
        classification_table = os.path.join(self.ts_assign_output, "final_outputs", "classifications.tsv")
        # Copy the classification table to replace after overwrite
        copyfile(classification_table, os.path.join(self.ts_assign_output, "tmp.tsv"))
        pre_lines = read_classification_table(get_test_data(classification_table))
        abundance_command_list = ["--treesapp_output", self.ts_assign_output,
                                  "--reads", get_test_data("test_TarA.1.fq"),
                                  "--reverse", get_test_data("test_TarA.2.fq"),
                                  "--pairing", "pe",
                                  "--report", "update",
                                  "--metric", "fpkm",
                                  "--num_procs", str(self.num_procs),
                                  "--overwrite"]
        abund_dict = abundance(abundance_command_list)["test_TarA"]
        post_lines = read_classification_table(get_test_data(classification_table))
        # Ensure no lines were added or removed
        self.assertEqual(len(pre_lines), len(post_lines))
        # Ensure the name of the sample is substituted for the sample ID
        self.assertEqual({"test_TarA"}, set([line[0] for line in post_lines]))
        # Ensure the value
        self.assertEqual(148, round(sum(abund_dict.values())))

        # Test an interleaved fastq file with TPM relative abundance measure
        abundance_command_list = ["--treesapp_output", self.ts_assign_output,
                                  "--reads", get_test_data("test_TarA.fastq.gz"),
                                  "--pairing", "pe",
                                  "--metric", "tpm",
                                  "--report", "append",
                                  "--num_procs", str(self.num_procs),
                                  "--delete", "--overwrite"]
        abundance(abundance_command_list)

        pqueries = assignments_to_pqueries(read_classification_table(get_test_data(classification_table)))
        abundant_pqueries = []
        for rp in pqueries:
            abundant_pqueries += [pq for pq in pqueries[rp] if pq.abundance > 0.0]
        self.assertEqual(33850, round(sum([pq.abundance for pq in abundant_pqueries])))

        # Replace the classification table
        copyfile(os.path.join(self.ts_assign_output, "tmp.tsv"), classification_table)
        os.remove(os.path.join(self.ts_assign_output, "tmp.tsv"))
        return

    def test_assign(self):
        ref_pkgs = ["McrA", "McrB", "DsrAB"]
        from treesapp import assign
        from .testing_utils import get_test_data
        from treesapp.file_parsers import read_classification_table
        from treesapp.phylo_seq import assignments_to_pqueries
        assign_commands_list = ["--fastx_input", self.aa_test_fa,
                                "--targets", ','.join(ref_pkgs),
                                "--refpkg_dir", self.refpkg_dir,
                                "--num_procs", str(self.num_procs),
                                "--hmm_coverage", str(20),
                                "--output", "./TreeSAPP_assign/",
                                "--stringency", "relaxed",
                                "--placement_summary", "max_lwr",
                                "--trim_align", "--overwrite", "--delete", "--svm"]
        assign.assign(assign_commands_list)
        lines = read_classification_table("./TreeSAPP_assign/final_outputs/classifications.tsv")
        self.assertEqual(15, len(lines))
        self.assertTrue(os.path.isfile("./TreeSAPP_assign/final_outputs/marker_test_suite_classified.faa"))

        # Test nucleotide sequence input
        ref_pkgs = ["McrA", "McrB"]
        assign_commands_list = ["--fastx_input", self.nt_test_fa,
                                "--targets", ','.join(ref_pkgs),
                                "--num_procs", str(self.num_procs),
                                "--output", "./TreeSAPP_assign/",
                                "--stringency", "strict",
                                "--trim_align", "--overwrite", "--delete"]
        assign.assign(assign_commands_list)

        assign_commands_list += ["--rel_abund",
                                 "--reads", get_test_data("SRR3669912_1.fastq"),
                                 "--reverse", get_test_data("SRR3669912_2.fastq")]
        assign.assign(assign_commands_list)
        assignments_tbl = "./TreeSAPP_assign/final_outputs/classifications.tsv"
        lines = read_classification_table(assignments_tbl)
        self.assertEqual(6, len(lines))
        classified_seqs = set()
        pqueries = assignments_to_pqueries(lines)
        for rp_pqs in pqueries.values():
            classified_seqs.update({pq.seq_name for pq in rp_pqs if pq.abundance >= 0.0})
        self.assertTrue("LYOS01000003.1:168824-170509_1" in classified_seqs)
        self.assertTrue(os.path.isfile("./TreeSAPP_assign/final_outputs/marker_test_suite_classified.fna"))
        return

    def test_colour(self):
        from treesapp.commands import colour
        # Test creating iTOL files for a taxonomic rank
        colour_commands = ["-r", self.mcra_pkl, self.mcrb_pkl,
                           "-l", "family",
                           "-o", "./TreeSAPP_colour",
                           "--filter", "Methanococcales",
                           "--min_proportion", str(0.01),
                           "--no_polyphyletic"]
        colour(colour_commands)
        self.assertTrue(os.path.isfile(os.path.join("TreeSAPP_colour", "McrA_family_colour_strip.txt")))

        # Test a failure when the desired attribute doesn't exist in the reference package(s)
        with pytest.raises(AssertionError):
            colour(["-r", self.mcra_pkl, self.mcrb_pkl,
                    "-o", "./TreeSAPP_colour",
                    "--palette", "viridis",
                    "--attribute", "Pathway"])

        # Test creating strip and styles file for an annotated feature
        colour_commands = ["-r", self.xmoa_pkl,
                           "-o", "./TreeSAPP_colour",
                           "--palette", "viridis",
                           "--attribute", "Paralog"]
        colour(colour_commands)
        itol_strip_file = os.path.join("TreeSAPP_colour", "XmoA_Paralog_colour_strip.txt")
        self.assertTrue(os.path.isfile(itol_strip_file))
        with open(itol_strip_file) as strip_dat:
            self.assertEqual(25, len(strip_dat.readlines()))
        return

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
                                "--screen", "Bacteria", "--filter", "Archaea",
                                "--min_taxonomic_rank", 'p',
                                "--output", "./TreeSAPP_create",
                                "--num_procs", str(self.num_procs),
                                "--max_examples", str(100),
                                "--trim_align", "--cluster", "--fast", "--headless",
                                "--overwrite", "--delete", "--deduplicate"]
        create(create_commands_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__pkl = "./TreeSAPP_create/final_outputs/Crt_build.pkl"
        test_refpkg.slurp()
        self.assertTrue(test_refpkg.validate())
        self.assertEqual(77, test_refpkg.num_seqs)
        self.assertEqual(2, len(test_refpkg.pfit))
        self.assertTrue(test_refpkg.pfit[0] < 0)

        # Test with different header format (EggNOG), infer phylogeny with RAxML-NG and skip training
        create_commands_list = ["--fastx_input", get_test_data("ENOG4111FIN.txt"),
                                "--output", "./TreeSAPP_create",
                                "--refpkg_name", "PuhA",
                                "--similarity", "0.90",
                                "--bootstraps", str(4),
                                "--profile", get_test_data("PuhA_search.hmm"),
                                "--screen", "Bacteria,Archaea",
                                "--num_procs", str(self.num_procs),
                                "--raxml_model", "LG+F+R4",
                                "--min_taxonomic_rank", 'p',
                                "--stage", "support",
                                "--trim_align", "--outdet_align", "--cluster", "--headless",
                                "--overwrite", "--delete"]
        create(create_commands_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__pkl = "./TreeSAPP_create/final_outputs/PuhA_build.pkl"
        test_refpkg.slurp()
        self.assertTrue(test_refpkg.validate())
        self.assertEqual(39, test_refpkg.num_seqs)
        self.assertEqual("LG+F+R4", test_refpkg.sub_model)

        # Test with another header format and the lineages are provided in an accession2lin table
        cmd_list = ["--fastx_input", get_test_data("PF01280_test.fasta"),
                    "--accession2lin", get_test_data("PF01280_test.tsv"),
                    "--output", "./TreeSAPP_create",
                    "--refpkg_name", "PF01280",
                    "--similarity", "0.95",
                    "--bootstraps", str(0),
                    "--screen", "Archaea", "--filter", "Bacteria,Eukaryota",
                    "--min_taxonomic_rank", 'g',
                    "--num_procs", str(self.num_procs),
                    "--stage", "evaluate",
                    "--trim_align", "--cluster", "--fast", "--headless", "--overwrite", "--delete"]
        create(cmd_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__pkl = "./TreeSAPP_create/final_outputs/PF01280_build.pkl"
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
                                 "--trim_align"]
        evaluate(evaluate_command_list)

        # Test rerunning of already processed data
        evaluate(evaluate_command_list)
        self.assertTrue(True)
        return

    def test_info(self):
        from treesapp.commands import info
        with pytest.raises(SystemExit):
            info(["--verbose", "--refpkg_dir", "./"])
        info(["--verbose"])
        self.assertTrue(True)
        return

    def test_layer(self):
        from treesapp.commands import layer
        from .testing_utils import get_test_data
        from treesapp.file_parsers import read_classification_table

        # Layering annotations from multiple reference packages without feature_annotations
        layer_command_list = ["--treesapp_output", self.ts_assign_output]
        layer(layer_command_list)

        # With a different reference package, XmoA, just to be sure
        classification_table = os.path.join(get_test_data("p_amoA_FunGene9.5_isolates_assign/"),
                                            "final_outputs", "classifications.tsv")
        layered_table = os.path.join(get_test_data("p_amoA_FunGene9.5_isolates_assign/"), "final_outputs",
                                     "layered_classifications.tsv")
        pre_lines = read_classification_table(get_test_data(classification_table))
        layer_command_list = ["--treesapp_output", get_test_data("p_amoA_FunGene9.5_isolates_assign/"),
                              "--refpkg_dir", self.refpkg_dir]
        layer(layer_command_list)
        post_lines = read_classification_table(get_test_data(classification_table))

        layered_classifications = read_classification_table(layered_table)
        self.assertEqual(len(pre_lines), len(post_lines))
        for line in layered_classifications:
            self.assertEqual(13, len(line))
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
        from treesapp.clade_annotation import CladeAnnotation
        from .testing_utils import get_test_data

        test_refpkg = ReferencePackage()
        output_dir = "./TreeSAPP_package" + os.sep
        test_refpkg.f__pkl = output_dir + "McrA_build.pkl"
        edit_command_list = ["edit",
                             "refpkg_code", "M0701",
                             "--refpkg_path", self.mcra_pkl,
                             "--output", output_dir]
        package(edit_command_list)

        # Test editing a component file's contents
        edit_command_list = ["edit",
                             "f__msa", self.aa_test_fa,
                             "--refpkg_path", self.mcra_pkl,
                             "--output", output_dir]
        with pytest.raises(SystemExit):
            package(edit_command_list)
        test_refpkg.slurp()
        self.assertTrue(test_refpkg.validate())
        with pytest.raises(SystemExit):
            package(["rename", "-h"])

        # Test adding to a reference package's feature_annotations attribute
        edit_command_list = ["edit",
                             "feature_annotations", "Function",
                             "--taxa_map", get_test_data("Mcr_taxonomy-phenotype_map.tsv"),
                             "--refpkg_path", self.mcra_pkl,
                             "--output", output_dir,
                             "--overwrite"]
        package(edit_command_list)
        test_refpkg.slurp()
        self.assertIsInstance(test_refpkg.feature_annotations, dict)
        self.assertTrue('Function' in test_refpkg.feature_annotations)
        self.assertEqual(5, len(test_refpkg.feature_annotations["Function"]))
        self.assertIsInstance(test_refpkg.feature_annotations["Function"][0], CladeAnnotation)
        for ca in test_refpkg.feature_annotations["Function"]:
            if ca.name == "Short-chain alkane oxidation":
                self.assertEqual(2, len(ca.members))

        # Test view with a complex attribute (feature_annotations)
        view_command_list = ["view",
                             "lineage_ids", "feature_annotations", "tree",
                             "--refpkg_path", output_dir + os.path.basename(self.mcra_pkl),
                             "--output", output_dir]
        package(view_command_list)

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
                              "--profile",
                              "--max_examples", str(max_ex),
                              "--num_proc", str(self.num_procs),
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

    def test_update(self):
        from treesapp.commands import update
        from treesapp.refpkg import ReferencePackage
        from .testing_utils import get_test_data
        # Test with re-training and the resolve workflow
        update_command_list = ["--fastx_input", get_test_data("ENOG4111FIN_PF03967_seed.faa"),
                               "--refpkg_path", self.puha_pkl,
                               "--treesapp_output", get_test_data("assign_SwissProt_PuhA/"),
                               "--output", "./TreeSAPP_update",
                               "--num_proc", str(self.num_procs),
                               "--trim_align", "--cluster", "--fast", "--headless",
                               "--overwrite", "--delete", "--skip_assign", "--resolve"]
        update(update_command_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__pkl = "./TreeSAPP_update/final_outputs/PuhA_build.pkl"
        test_refpkg.slurp()
        self.assertTrue(test_refpkg.validate())
        self.assertEqual(49, test_refpkg.num_seqs)
        self.assertEqual(2, len(test_refpkg.feature_annotations["test"]))

        # Test the workflow when seqs2lineage table is provided and no training
        update_command_list = ["--refpkg_path", self.puha_pkl,
                               "--treesapp_output", get_test_data("assign_SwissProt_PuhA/"),
                               "--seqs2lineage", get_test_data("SwissProt_PuhA_seqs2lineage.txt"),
                               "--output", "./TreeSAPP_update",
                               "--num_procs", str(self.num_procs),
                               "-b", str(0),
                               "--trim_align", "--cluster", "--fast", "--headless",
                               "--overwrite", "--delete", "--skip_assign"]
        update(update_command_list)
        test_refpkg = ReferencePackage()
        test_refpkg.f__pkl = "./TreeSAPP_update/final_outputs/PuhA_build.pkl"
        test_refpkg.slurp()
        self.assertTrue(test_refpkg.validate())
        self.assertEqual(49, test_refpkg.num_seqs)
        self.assertEqual((-2.7691, 7.0), test_refpkg.pfit)
        self.assertTrue(test_refpkg.svc is not None)
        return

    # def test_tmp(self):
    #     from treesapp.commands import create
    #     base_dir = "/home/connor/Bioinformatics/Hallam_projects/RefPkgs/"
    #     cmd = "".format(base_dir)
    #     create(cmd.split())
    #     return


if __name__ == '__main__':
    unittest.main()

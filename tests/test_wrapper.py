import unittest
import os
import shutil
import pytest

from .testing_utils import get_test_data, get_treesapp_root


class ExecutableWrapperTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp import utilities
        self.num_procs = max(int(utilities.available_cpu_count()/2), 2)
        self.ts_dir = get_treesapp_root()
        self.tmp_dir = os.path.join(os.getcwd(), "tests", "wrapper_test_dir") + os.sep
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.mkdir(self.tmp_dir)

        self.test_fasta = get_test_data("create_test.faa")
        return

    def tearDown(self) -> None:
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        return

    def test_run_linclust(self):
        from treesapp.utilities import fetch_executable_path
        from treesapp.wrapper import run_linclust
        from treesapp.fasta import read_fasta_to_dict
        # Fail due to fasta file being a string and not a list (and hence doesn't exist)
        prefix = os.path.join(self.tmp_dir, "linclust_test")
        with pytest.raises(SystemExit):
            run_linclust(mmseqs_exe=fetch_executable_path("mmseqs", self.ts_dir),
                         fa_in=self.test_fasta, output_prefix=prefix, prop_sim=0.99)

        # Fail due to out of range alignment similarity
        with pytest.raises(SystemExit):
            run_linclust(mmseqs_exe=fetch_executable_path("mmseqs", self.ts_dir),
                         fa_in=[self.test_fasta], output_prefix=prefix, prop_sim=50)

        # Run as normal
        run_linclust(mmseqs_exe=fetch_executable_path("mmseqs", self.ts_dir),
                     fa_in=[self.test_fasta], output_prefix=prefix, prop_sim=0.70)

        # Ensure the necessary files exist and are not empty
        rep_fasta = prefix + "_rep_seq.fasta"
        self.assertTrue(os.path.exists(rep_fasta))
        self.assertEqual(91, len(read_fasta_to_dict(rep_fasta)))

        aln_table = prefix + "_cluster_aln.tsv"
        self.assertTrue(os.path.exists(aln_table))

        cluster_table = prefix + "_cluster.tsv"
        self.assertTrue(os.path.exists(cluster_table))
        return

    def test_cluster_sequences(self):
        from treesapp.utilities import fetch_executable_path
        from treesapp.wrapper import cluster_sequences
        # VSEARCH isn't required anymore so test if its installed
        try:
            vsearch = fetch_executable_path("vsearch", self.ts_dir)
        except SystemExit:
            return
        cluster_sequences(software_path=vsearch, similarity=0.90,
                          fasta_input=self.test_fasta, output_prefix=os.path.join(self.tmp_dir, "vsearch_test"))
        self.assertTrue(os.path.isfile(os.path.join(self.tmp_dir, "vsearch_test.uc")))
        return

    def test_support_tree_raxml(self):
        from treesapp.utilities import fetch_executable_path
        from treesapp.wrapper import support_tree_raxml
        from treesapp.entish import load_ete3_tree
        tree = get_test_data("PuhA.raxml.bestTree")
        msa = get_test_data("PuhA.phy")
        f_support = support_tree_raxml(raxml_exe=fetch_executable_path(exe_name="raxml-ng", treesapp_dir=self.ts_dir),
                                       ref_tree=tree, ref_msa=msa,
                                       model="LG", tree_prefix=os.path.join(self.tmp_dir, "PuhA"),
                                       metric="tbe",
                                       n_bootstraps=4, num_threads=self.num_procs, mre=False)
        self.assertTrue(os.path.isfile(f_support))
        bs_tree = load_ete3_tree(f_support)
        self.assertEqual(39, len(bs_tree))
        return

    def test_construct_tree(self):
        from treesapp import utilities
        from treesapp import wrapper
        from treesapp import entish
        raxml_exe = utilities.fetch_executable_path(exe_name="raxml-ng", treesapp_dir=self.ts_dir)
        # Test with an absurd number of threads to ensure RAxML-NG's auto-scaling works
        best_tree = wrapper.construct_tree(tree_builder="RAxML-NG",
                                           executables={"raxml-ng": raxml_exe},
                                           evo_model="WAG+R2",
                                           multiple_alignment_file=get_test_data("PuhA.phy"),
                                           tree_output_dir=self.tmp_dir,
                                           tree_prefix="TMP",
                                           num_trees=1,
                                           num_threads=self.num_procs*2,
                                           verbosity=0)
        self.assertTrue(os.path.isfile(best_tree))
        bs_tree = entish.load_ete3_tree(best_tree)
        self.assertEqual(39, len(bs_tree))
        # self.assertEqual(logging.INFO,
        #                  logging.getLogger().level)
        return

    def test_build_hmm_profile(self):
        from treesapp.wrapper import build_hmm_profile
        from treesapp.utilities import fetch_executable_path
        # Create an output directory with parentheses, which ordinarily trips up HMMER
        os.makedirs(os.path.join(self.tmp_dir, "s__Leucobacter_sp._7(1)"))
        test_hmm = os.path.join(self.tmp_dir, "s__Leucobacter_sp._7(1)", "test.hmm")
        build_hmm_profile(hmmbuild_exe=fetch_executable_path("hmmbuild", self.ts_dir),
                          msa_in=get_test_data("PuhA.mfa"),
                          output_hmm=test_hmm)
        self.assertTrue(os.path.isfile(test_hmm))
        return


if __name__ == '__main__':
    unittest.main()

import unittest
import os
import shutil
import pytest

from .testing_utils import get_test_data, get_treesapp_root


class ExecutableWrapperTester(unittest.TestCase):
    def setUp(self) -> None:
        self.num_procs = 2
        self.ts_dir = get_treesapp_root()
        self.tmp_dir = os.path.join(os.getcwd(), "tests", "wrapper_test_dir")
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
        tree = get_test_data("PuhA.raxml.bestTree")
        msa = get_test_data("PuhA.phy")
        f_support = support_tree_raxml(raxml_exe=fetch_executable_path(exe_name="raxml-ng", treesapp_dir=self.ts_dir),
                                       ref_tree=tree, ref_msa=msa,
                                       model="LG", tree_prefix=os.path.join(self.tmp_dir, "PuhA"),
                                       n_bootstraps=4, num_threads=self.num_procs, mre=False)
        self.assertTrue(os.path.isfile(f_support))
        return

    def test_run_graftm_graft(self):
        from treesapp.utilities import fetch_executable_path
        from treesapp.wrapper import run_graftm_graft
        try:
            graftm_exe = fetch_executable_path("graftM", self.ts_dir)
        except SystemExit:
            return
        classification_tbl = os.path.join(self.tmp_dir, "EggNOG_McrA", "EggNOG_McrA_read_tax.tsv")

        # Run graftM graft
        run_graftm_graft(graftm_exe,
                         input_path=get_test_data("EggNOG_McrA.faa"),
                         output_dir=self.tmp_dir,
                         gpkg_path=get_test_data(os.path.join("refpkgs", "7.27_mcrA.gpkg")),
                         classifier="graftm",
                         num_threads=self.num_procs)

        self.assertTrue(os.path.isfile(classification_tbl))
        with open(classification_tbl) as tbl_handler:
            tbl_lines = tbl_handler.readlines()
        self.assertEqual(54, len(tbl_lines))
        return


if __name__ == '__main__':
    unittest.main()

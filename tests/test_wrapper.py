import unittest
import os
import shutil
import pytest

from .testing_utils import get_test_data, get_treesapp_root


class ExecutableWrapperTester(unittest.TestCase):
    def setUp(self) -> None:
        self.ts_dir = get_treesapp_root()
        self.tmp_dir = os.path.join(os.getcwd(), "wrapper_test_dir")
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
        # Fail due to temporary directory not being provided
        with pytest.raises(SystemExit):
            run_linclust(mmseqs_exe=fetch_executable_path("mmseqs", self.ts_dir),
                         fa_in=self.test_fasta,
                         output_prefix=os.path.join(self.tmp_dir, "linclust_test"),
                         prop_sim=0.99,
                         tmp_dir="")

        # Fail due to out of range alignment similarity
        with pytest.raises(SystemExit):
            run_linclust(mmseqs_exe=fetch_executable_path("mmseqs", self.ts_dir),
                         fa_in=self.test_fasta,
                         output_prefix=os.path.join(self.tmp_dir, "linclust_test"),
                         prop_sim=50,
                         tmp_dir=self.tmp_dir)

        # Run as normal
        run_linclust(mmseqs_exe=fetch_executable_path("mmseqs", self.ts_dir),
                     fa_in=self.test_fasta,
                     output_prefix=os.path.join(self.tmp_dir, "linclust_test"),
                     prop_sim=0.80,
                     tmp_dir=self.tmp_dir)
        # Ensure the necessary files exist and are not empty
        rep_fasta = os.path.join(self.tmp_dir, "linclust_test_rep_seq.fasta")
        cluster_table = os.path.join(self.tmp_dir, "linclust_test_cluster.tsv")
        self.assertTrue(os.path.exists(rep_fasta))
        self.assertTrue(os.path.exists(cluster_table))
        self.assertEqual(97, len(read_fasta_to_dict(rep_fasta)))
        return

    def test_cluster_sequences(self):
        from treesapp.utilities import fetch_executable_path
        from treesapp.wrapper import cluster_sequences
        cluster_sequences(software_path=fetch_executable_path("vsearch", self.ts_dir),
                          fasta_input=self.test_fasta, output_prefix=os.path.join(self.tmp_dir, "vsearch_test"))
        self.assertTrue(os.path.isfile(os.path.join(self.tmp_dir, "vsearch_test.uc")))
        return


if __name__ == '__main__':
    unittest.main()

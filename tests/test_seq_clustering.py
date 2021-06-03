import os
import shutil
import unittest
import pytest

from .testing_utils import get_test_data


class ClusteringTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_dir = os.path.join(os.getcwd(), "tests", "phyloclust_test_dir")
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.mkdir(self.tmp_dir)

    def tearDown(self) -> None:
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        return

    def test_read_linclust_clusters(self):
        from treesapp.seq_clustering import read_linclust_clusters
        with pytest.raises(SystemExit):
            read_linclust_clusters(get_test_data("linclust_test_cluster_aln.tsv"))
        cluster_dict = read_linclust_clusters(get_test_data("linclust_test_cluster.tsv"))
        self.assertEqual(86, len(cluster_dict))
        self.assertEqual(85, max([int(x) for x in cluster_dict.keys()]))  # Ensure the cluster names are zero-based
        return

    def test_create_mmseqs_clusters(self):
        from treesapp.seq_clustering import create_mmseqs_clusters
        clusters_dict = create_mmseqs_clusters(clusters_tbl=get_test_data("linclust_test_cluster.tsv"),
                                               aln_tbl=get_test_data("linclust_test_cluster_aln.tsv"))
        self.assertEqual(86, len(clusters_dict))
        return

    def test_read_uc(self):
        from treesapp.seq_clustering import read_uc
        # Fail if input doesn't exist
        with pytest.raises(SystemExit):
            read_uc("tests/tmp.uc.txt")

        # Run as normal
        cluster_dict = read_uc(get_test_data("vsearch_test.uc"))
        self.assertEqual(98, len(cluster_dict))
        return

    def test_dereplicate_by_clustering(self):
        from treesapp.seq_clustering import dereplicate_by_clustering
        from treesapp import fasta
        from treesapp.utilities import fetch_executable_path
        test_fa = fasta.FASTA(get_test_data("PuhA.fa"))
        test_fa.load_fasta()
        cluster_map = dereplicate_by_clustering(fasta_inst=test_fa,
                                                prop_similarity=0.799,
                                                mmseqs_exe=fetch_executable_path("mmseqs", treesapp_dir="./"),
                                                tmp_dir=self.tmp_dir)
        self.assertEqual(len(cluster_map), test_fa.n_seqs())
        return


if __name__ == '__main__':
    unittest.main()

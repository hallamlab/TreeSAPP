import os
import unittest
import pytest

from .testing_utils import get_test_data, get_treesapp_root


class TreesappTester(unittest.TestCase):
    def test_gather_ref_packages(self):
        from treesapp.file_parsers import gather_ref_packages
        from . import testing_utils as utils
        refpkg_dir = utils.get_test_data(os.path.join("refpkgs"))
        #
        refpkg_dict = gather_ref_packages(refpkg_data_dir=refpkg_dir, targets=["McrA"])
        self.assertEqual(1, len(refpkg_dict))

        # Test to ensure all core reference packages are present in the main refpkg pkl directory
        core = {"McrA", "McrB", "McrG", "DsrAB", "XmoA",
                "NapA", "NxrA", "NifD", "NirK", "NirS", "NxrB", "NorB", "NosZ",
                "HydA", "PilA"}
        refpkg_dict = gather_ref_packages(refpkg_data_dir=os.path.join(get_treesapp_root(), "data"))
        self.assertEqual(0, len(core.difference(set(refpkg_dict.keys()))))
        return

    def test_read_linclust_clusters(self):
        from treesapp.file_parsers import read_linclust_clusters
        with pytest.raises(SystemExit):
            read_linclust_clusters(get_test_data("linclust_test_cluster_aln.tsv"))
        cluster_dict = read_linclust_clusters(get_test_data("linclust_test_cluster.tsv"))
        self.assertEqual(86, len(cluster_dict))
        self.assertEqual(85, max([int(x) for x in cluster_dict.keys()]))  # Ensure the cluster names are zero-based
        return

    def test_create_mmseqs_clusters(self):
        from treesapp.file_parsers import create_mmseqs_clusters
        clusters_dict = create_mmseqs_clusters(clusters_tbl=get_test_data("linclust_test_cluster.tsv"),
                                               aln_tbl=get_test_data("linclust_test_cluster_aln.tsv"))
        self.assertEqual(86, len(clusters_dict))
        return

    def test_read_uc(self):
        from treesapp.file_parsers import read_uc
        # Fail if input doesn't exist
        with pytest.raises(SystemExit):
            read_uc("tests/tmp.uc.txt")

        # Run as normal
        cluster_dict = read_uc(get_test_data("vsearch_test.uc"))
        self.assertEqual(98, len(cluster_dict))
        return

    def test_read_graftm_classifications(self):
        from treesapp.file_parsers import read_graftm_classifications
        return

    def test_grab_graftm_taxa(self):
        from treesapp.file_parsers import grab_graftm_taxa
        return

    def test_best_discrete_matches(self):
        from treesapp.file_parsers import best_discrete_matches
        from treesapp.hmmer_tbl_parser import HmmMatch
        # Create test HmmMatch instances
        m1 = HmmMatch()
        m2 = HmmMatch()
        m3 = HmmMatch()
        m1.start, m1.end, m1.target_hmm = 80, 200, "NorB"
        m2.start, m2.end, m2.target_hmm = 15, 220, "NorC"
        m3.start, m3.end, m3.target_hmm = 240, 391, "NorB"

        test_matches = [m1, m2, m3]
        best_matches = best_discrete_matches(test_matches)
        self.assertEqual(2, len(best_matches))
        self.assertEqual(3, len(test_matches))
        self.assertTrue(m2 in best_matches and m3 in best_matches)
        return

    def test_read_phylip_to_dict(self):
        from treesapp.file_parsers import read_phylip_to_dict
        from .testing_utils import get_test_data
        with pytest.raises(SystemExit):
            read_phylip_to_dict("test.phy")

        phy_dict = read_phylip_to_dict(get_test_data("test.phy"))
        self.assertEqual(5, len(phy_dict))
        self.assertEqual(76, len(phy_dict["cox2_leita"]))
        return

    def test_read_stockholm_to_dict(self):
        from treesapp.file_parsers import read_stockholm_to_dict
        from .testing_utils import get_test_data
        with pytest.raises(SystemExit):
            read_stockholm_to_dict("test.sto")
        sto_dict = read_stockholm_to_dict(get_test_data("test.sto"))
        self.assertEqual(236, len(sto_dict))
        self.assertEqual(680, len(sto_dict["KKH90701"]))
        return

    def test_load_classifified_sequences_from_assign_output(self):
        from treesapp.file_parsers import load_classified_sequences_from_assign_output
        refpkg_pquery_map = load_classified_sequences_from_assign_output(get_test_data("marker_test_results"))
        self.assertEqual(["McrA", "McrB"], list(refpkg_pquery_map.keys()))
        self.assertEqual(13, sum([len(pqueries) for pqueries in refpkg_pquery_map.values()]))
        for rp, pqueries in refpkg_pquery_map.items():
            for pq in pqueries:
                self.assertTrue(len(pq.seq) > 0)
                self.assertEqual(rp, pq.ref_name)

        # Test refpkg filtering
        refpkg_pquery_map = load_classified_sequences_from_assign_output(get_test_data("test_output_TarA"),
                                                                         refpkg_name="DsrAB")
        self.assertEqual(70, len(refpkg_pquery_map["DsrAB"]))
        for rp, pqueries in refpkg_pquery_map.items():
            for pq in pqueries:
                self.assertEqual("DsrAB", pq.ref_name)

        return

    def test_write_classification_table(self):
        from treesapp.file_parsers import write_classification_table, load_classified_sequences_from_assign_output
        refpkg_pquery_map = load_classified_sequences_from_assign_output(get_test_data("marker_test_results"))
        test_output = "tests/tmp_classifications.tsv"
        write_classification_table(tree_saps=refpkg_pquery_map, output_file=test_output, sample_name="Testing")
        with open(test_output) as t_out:
            self.assertEqual(14, len(t_out.readlines()))

        # Try appending
        write_classification_table(tree_saps=refpkg_pquery_map,
                                   sample_name="Testing_testing",
                                   output_file=test_output,
                                   append=True)
        with open(test_output) as t_out:
            self.assertEqual(27, len(t_out.readlines()))

        # Clean up output
        if os.path.isfile(test_output):
            os.remove(test_output)
        return


if __name__ == '__main__':
    unittest.main()

import unittest
import pytest

from .testing_utils import get_test_data, get_treesapp_root


class TreesappTester(unittest.TestCase):
    def test_gather_ref_packages(self):
        from treesapp.file_parsers import gather_ref_packages
        import os
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


if __name__ == '__main__':
    unittest.main()

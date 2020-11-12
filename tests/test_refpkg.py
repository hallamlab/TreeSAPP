import unittest
import pytest
import os
from shutil import rmtree


@pytest.fixture(scope="class")
def refpkg_class(request):
    from treesapp import refpkg
    from . import testing_utils as utils
    request.cls.db = refpkg.ReferencePackage("McrA")
    request.cls.db.f__json = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
    request.cls.db.slurp()
    return


@pytest.mark.usefixtures("refpkg_class")
class RefPkgTester(unittest.TestCase):
    def setUp(self) -> None:
        self.new_pkl_path = "./test_write_json" + self.db.refpkg_suffix
        self.disband_path = "_".join([self.db.prefix, self.db.refpkg_code, self.db.date])
        if os.path.isdir(self.disband_path):
            rmtree(self.disband_path)

        return

    def tearDown(self) -> None:
        if os.path.isdir(self.disband_path):
            rmtree(self.disband_path)
        if os.path.isfile(self.new_pkl_path):
            os.remove(self.new_pkl_path)
        return

    def test_band(self):
        self.db.band()
        return

    def test_slurp(self):
        from . import testing_utils as utils
        self.db.f__json = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.db.slurp()
        self.assertEqual("McrA", self.db.prefix)
        return

    def test_disband(self):
        self.db.disband("./")
        self.assertTrue(os.path.isfile(os.path.join(self.disband_path, "McrA.fa")))

    def test_remove_taxon_from_lineage_ids(self):
        # Ensure the initial state is as expected
        self.assertEqual(246, self.db.num_seqs)
        self.assertEqual(246, len(self.db.lineage_ids))

        # Test
        self.db.remove_taxon_from_lineage_ids("d__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales;"
                                              " f__Methanobacteriaceae; g__Methanosphaera")
        self.assertEqual(193, self.db.num_seqs)
        self.db.remove_taxon_from_lineage_ids("d__Archaea; p__Euryarchaeota; c__Methanobacteria")
        self.assertEqual(146, len(self.db.lineage_ids))

        self.db.remove_taxon_from_lineage_ids("d__Archaea")
        self.assertEqual(3, len(self.db.lineage_ids))
        self.assertEqual(3, self.db.num_seqs)

    def test_get_fasta(self):
        refpkg_fa = self.db.get_fasta()
        self.assertEqual(246, len(refpkg_fa.fasta_dict))

    def test_get_ete_tree(self):
        from treesapp.refpkg import ReferencePackage
        rt = self.db.get_ete_tree()
        self.assertEqual(246, len(rt))
        blank = ReferencePackage()
        with pytest.raises(SystemExit):
            blank.get_ete_tree()
        return

    def test_hmm_length(self):
        from treesapp.refpkg import ReferencePackage
        blank = ReferencePackage()
        self.assertEqual(blank.profile, [])
        with pytest.raises(ValueError):
            blank.hmm_length()
        return

    def test_bail(self):
        self.db.bail()

    def test_get_internal_node_leaf_map(self):
        from . import testing_utils as utils
        self.db.f__json = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.db.slurp()
        node_map = self.db.get_internal_node_leaf_map()
        self.assertEqual(2*self.db.num_seqs, len(node_map))
        self.assertEqual(self.db.num_seqs, len(node_map[max(node_map.keys())]))
        return

    def test_pickle_package(self):
        self.db.f__json = self.new_pkl_path
        self.db.pickle_package()
        self.db.slurp()
        self.assertTrue("McrA" == self.db.prefix)

    def test_taxonomically_label_tree(self):
        from treesapp.taxonomic_hierarchy import Taxon
        labelled_rt = self.db.taxonomically_label_tree()
        self.assertEqual(246, len(labelled_rt))
        self.assertTrue('taxon' in labelled_rt.features)
        self.assertEqual('Root', labelled_rt.taxon.name)
        self.assertEqual('root', labelled_rt.taxon.rank)
        self.assertTrue(isinstance(labelled_rt.taxon, Taxon))
        return

    def test_enumerate_taxonomic_lineages(self):
        from treesapp.refpkg import ReferencePackage
        mock_rp = ReferencePackage()
        self.assertEqual(0, len(mock_rp.enumerate_taxonomic_lineages()))
        taxa_counts = self.db.enumerate_taxonomic_lineages()
        self.assertIsInstance(taxa_counts, dict)
        self.assertTrue("r__Root" in taxa_counts)
        return


if __name__ == '__main__':
    unittest.main()

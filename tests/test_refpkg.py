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
        self.db.remove_taxon_from_lineage_ids("d__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales;"
                                              " f__Methanobacteriaceae; g__Methanosphaera")
        self.assertEqual(190, self.db.num_seqs)
        self.db.remove_taxon_from_lineage_ids("d__Archaea; p__Euryarchaeota; c__Methanobacteria")
        self.assertEqual(143, len(self.db.lineage_ids))

        self.db.remove_taxon_from_lineage_ids("d__Archaea")
        self.assertEqual(0, len(self.db.lineage_ids))
        self.assertEqual(0, self.db.num_seqs)

    def test_get_fasta(self):
        refpkg_fa = self.db.get_fasta()
        self.assertEqual(246, len(refpkg_fa.fasta_dict))

    def test_bail(self):
        self.db.bail()

    def test_get_internal_node_leaf_map(self):
        from . import testing_utils as utils
        self.db.f__json = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.db.slurp()
        node_map = self.db.get_internal_node_leaf_map()
        self.assertEqual((2*self.db.num_seqs)-1, len(node_map))
        self.assertEqual(self.db.num_seqs, len(node_map[max(node_map.keys())]))

    def test_pickle_package(self):
        self.db.f__json = self.new_pkl_path
        self.db.pickle_package()
        self.db.slurp()
        self.assertTrue("McrA" == self.db.prefix)


if __name__ == '__main__':
    unittest.main()

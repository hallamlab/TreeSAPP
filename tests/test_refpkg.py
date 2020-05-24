import unittest
import pytest


@pytest.fixture(scope="class")
def refpkg_class(request):
    from treesapp import refpkg
    from . import testing_utils as utils
    request.cls.db = refpkg.ReferencePackage("McrA")
    request.cls.db.f__json = utils.get_test_data("band_test.json")
    request.cls.db.slurp()
    return


@pytest.mark.usefixtures("refpkg_class")
class RefPkgTester(unittest.TestCase):
    def test_band(self):
        self.db.band()
        return

    def test_slurp(self):
        self.db.slurp()
        self.assertEqual("McrA", self.db.prefix)
        return

    def test_disband(self):
        import os
        self.db.disband("./")
        self.assertTrue(os.path.isfile("./McrA__/McrA.fa"))

    def test_remove_taxon_from_lineage_ids(self):
        self.db.remove_taxon_from_lineage_ids("d__Archaea; p__Euryarchaeota; c__Methanobacteria")
        self.assertEqual(119, len(self.db.lineage_ids))
        self.assertEqual(119, self.db.num_seqs)

        self.db.remove_taxon_from_lineage_ids("d__Archaea")
        self.assertEqual(0, len(self.db.lineage_ids))
        self.assertEqual(0, self.db.num_seqs)

if __name__ == '__main__':
    unittest.main()
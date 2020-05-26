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
        from . import testing_utils as utils
        self.db.f__json = utils.get_test_data("band_test.json")
        self.db.slurp()
        self.assertEqual("McrA", self.db.prefix)
        return

    def test_disband(self):
        import os
        self.db.disband("./")
        self.assertTrue(os.path.isfile("./McrA_M0701_/McrA.fa"))

    def test_remove_taxon_from_lineage_ids(self):
        self.db.remove_taxon_from_lineage_ids("d__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales;"
                                              " f__Methanobacteriaceae; g__Methanosphaera")
        self.assertEqual(168, self.db.num_seqs)
        self.db.remove_taxon_from_lineage_ids("d__Archaea; p__Euryarchaeota; c__Methanobacteria")
        self.assertEqual(121, len(self.db.lineage_ids))

        self.db.remove_taxon_from_lineage_ids("d__Archaea")
        self.assertEqual(0, len(self.db.lineage_ids))
        self.assertEqual(0, self.db.num_seqs)

    def test_get_fasta(self):
        refpkg_fa = self.db.get_fasta()
        self.assertEqual(246, len(refpkg_fa.fasta_dict))

    def test_bail(self):
        self.db.bail()


if __name__ == '__main__':
    unittest.main()

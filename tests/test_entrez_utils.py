import unittest
from Bio import Entrez


class MyTestCase(unittest.TestCase):
    def test_repair_lineages(self):
        self.assertEqual(True, False)

    def test_fetch_entrez_lineages(self):
        from . import testing_utils as utils
        from treesapp.classy import Creator
        from treesapp.fasta import FASTA
        test_fa = FASTA(utils.get_test_data("create_test.faa"))
        test_fa.load_fasta()
        test_ts = Creator()
        test_ts.fetch_entrez_lineages(test_fa, 'prot')

    def test_prep_for_entrez_query(self):
        from treesapp.entrez_utils import prep_for_entrez_query
        handle = prep_for_entrez_query()
        for record in Entrez.read(handle):
            self.assertTrue(record["TaxId"] == "158330")
            self.assertEqual("Cypripedioideae", record["ScientificName"])
        return


if __name__ == '__main__':
    unittest.main()

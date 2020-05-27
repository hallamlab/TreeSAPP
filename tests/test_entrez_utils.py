import unittest


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


if __name__ == '__main__':
    unittest.main()

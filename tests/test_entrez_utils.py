import os
import unittest

from Bio import Entrez

from . import testing_utils as utils


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.fasta import FASTA
        from treesapp.classy import Creator

        self.test_fa = FASTA(utils.get_test_data("create_test.faa"))
        self.test_fa.load_fasta()

        self.create_inst = Creator()
        self.create_inst.current_stage = self.create_inst.stages[1]
        self.create_inst.acc_to_lin = "./test_create_acc_to_lin.tsv"
        return

    def tearDown(self) -> None:
        if os.path.isfile(self.create_inst.acc_to_lin):
            os.remove(self.create_inst.acc_to_lin)

    def test_repair_lineages(self):
        # TODO: provide EntrezRecord instances in ref_seq_dict
        from treesapp.entrez_utils import repair_lineages
        repair_lineages(ref_seq_dict={}, t_hierarchy=self.create_inst.ref_pkg.taxa_trie)
        return

    def test_fetch_entrez_lineages(self):
        self.create_inst.fetch_entrez_lineages(self.test_fa, 'prot')
        return

    def test_prep_for_entrez_query(self):
        from treesapp.entrez_utils import prep_for_entrez_query
        handle = prep_for_entrez_query()
        for record in Entrez.read(handle):
            self.assertTrue(record["TaxId"] == "158330")
            self.assertEqual("Cypripedioideae", record["ScientificName"])
        return

    def test_lineage_format_check(self):
        from treesapp.entrez_utils import Lineage
        tlin = Lineage()
        tlin.Lineage = ""
        self.assertFalse(tlin.lineage_format_check())
        tlin.Lineage = "r__Root"
        self.assertFalse(tlin.lineage_format_check())
        tlin.Lineage = "f__Methanomassiliicoccaceae;g__Candidatus Methanoplasma;s__Candidatus Methanoplasma termitum"
        self.assertTrue(tlin.lineage_format_check())
        self.assertEqual(3, len(tlin.Lineage.split(tlin.lin_sep)))
        # Retry now that tlin.Lineage has been reformatted
        self.assertFalse(tlin.lineage_format_check())
        return

    def test_verify_rank_occupancy(self):
        from treesapp.entrez_utils import Lineage
        tlin = Lineage()
        incomplete_lin = "d__Archaea; p__Euryarchaeota; c__Methanomicrobia; o__; f__; g__Candidatus Syntrophoarchaeum"
        complete_lin = "d__Domain; p__Phylum; c__Class; o__Order; f__Family"
        no_prefix_lin = "Archaea; Euryarchaeota; Methanomicrobia"

        # Test one: Incomplete lineage should be truncated
        tlin.Lineage = incomplete_lin
        tlin.verify_rank_occupancy()
        self.assertEqual("d__Archaea; p__Euryarchaeota; c__Methanomicrobia", tlin.Lineage)
        # Test two: complete lineage should be the same before and after
        tlin.Lineage = complete_lin
        self.assertEqual(complete_lin, tlin.Lineage)
        # Test three: no prefix
        tlin.Lineage = no_prefix_lin
        self.assertEqual(no_prefix_lin, tlin.Lineage)

        return


if __name__ == '__main__':
    unittest.main()

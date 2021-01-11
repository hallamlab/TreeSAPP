import os
import unittest

from .testing_utils import get_test_data


class UpdaterTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp import refpkg as ts_rp
        self.ref_pkg = ts_rp.ReferencePackage()
        self.ref_pkg.f__json = get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.ref_pkg.slurp()
        return

    def test_decide_length_filter(self):
        from treesapp import update_refpkg as ts_update_mod
        min_length = ts_update_mod.decide_length_filter(self.ref_pkg, proposed_min_length=800)
        self.assertEqual(800, min_length)

        # Ensure the incorrectly-provided percentage is converted to a proportion
        min_length = ts_update_mod.decide_length_filter(self.ref_pkg, min_hmm_proportion=80)
        self.assertEqual(443, min_length)
        return

    def test_simulate_entrez_records(self):
        from treesapp.update_refpkg import simulate_entrez_records
        from treesapp.fasta import FASTA, register_headers
        # Set up test data
        fa = FASTA("dummy.faa")
        fa.fasta_dict = {"AB-750_E14_NODE_11_length_29119_cov_8.74474_ID_15_7 # 6278 # 6853 # -1|HgcB|4_117": "MA",
                         "AB-750_E14_NODE_11_length_29119_cov_8.74474_ID_15_7 # 6278 # 6853 # -1|HgcB|119_185": "MA",
                         "2509738965_Desulfomonile_tiedjei_DCB|HgcB|1_200": "NA"}
        fa.header_registry = register_headers(list(fa.fasta_dict.keys()))
        lin_map = {"AB-750_E14_NODE_11_length_29119_cov_8.74474_ID_15_7": "d__Bacteria;p__Proteobacteria",
                   "2509738965_Desulfomonile_tiedjei_DCB|HgcB|1_200": "d__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; o__Syntrophobacterales; f__Syntrophaceae"}
        fa.change_dict_keys("num")

        # Test to ensure sequences with same prefix have different EntrezRecord instances
        entrez_records = simulate_entrez_records(fasta_records=fa, seq_lineage_map=lin_map)
        self.assertEqual(3, len(entrez_records))
        self.assertEqual(3, len(set(entrez_records.values())))
        self.assertEqual("# 6278 # 6853 # -1|HgcB|4_117", entrez_records['1'].description)
        tmp = [er.rebuild_header() for (num_id, er) in entrez_records.items()]
        self.assertTrue("AB-750_E14_NODE_11_length_29119_cov_8.74474_ID_15_7 # 6278 # 6853 # -1|HgcB|119_185" in tmp)
        return


if __name__ == '__main__':
    unittest.main()

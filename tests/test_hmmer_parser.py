import unittest
from collections import namedtuple

args_tuple = namedtuple("args_tuple", ["max_e", "max_ie", "min_acc", "min_score", "perc_aligned", "query_aligned"])


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from .testing_utils import get_test_data
        from treesapp.hmmer_tbl_parser import prep_args_for_parsing
        args = args_tuple(max_e=1E-5, max_ie=1E-3, min_acc=0.7, min_score=20, perc_aligned=60, query_aligned=80)
        self.thresholds = prep_args_for_parsing(args)
        self.nxra_domtbl = get_test_data("NxrA_search_to_ORFs_domtbl.txt")
        self.norc_domtbl = get_test_data("NorC_search_to_ORFs_domtbl.txt")
        return

    def test_parse_domain_tables(self):
        from treesapp import file_parsers
        hmm_matches = file_parsers.parse_domain_tables(thresholds=self.thresholds,
                                                       hmm_domtbl_files={("NxrA", "ORFs"): self.nxra_domtbl,
                                                                         ("NorC", "ORFs"): self.norc_domtbl})
        self.assertTrue('NorC' in hmm_matches.keys())
        orf_names = {match.orf for match in hmm_matches["NxrA"]}
        self.assertEqual(5, len(hmm_matches["NorC"]))
        self.assertEqual({'455656', '1982434', '612303', '86147', '1717364', '30974', '1799218'}, orf_names)

        return


if __name__ == '__main__':
    unittest.main()

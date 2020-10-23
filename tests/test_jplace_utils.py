import unittest
from .testing_utils import get_test_data
from treesapp.phylo_seq import PhyloPlace


class JPlaceTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.test_jplace = get_test_data("epa_result.jplace")
        self.placement_dict = {'p': [[0, -56570.7911467148, 1.0, 0.859477616, 1.2276349525]],
                               'n': ['scaffold_52752_c1_1 partial=11;start_type=Edge;rbs_motif=None|McrA|1_204']}
        return

    def test_jplace_parser(self):
        from treesapp.jplace_utils import jplace_parser
        jplace_dat = jplace_parser(self.test_jplace)
        self.assertEqual(3, len(jplace_dat.placements))
        self.assertIsInstance(jplace_dat.placements, list)
        return

    def test_demultiplex_pqueries(self):
        from treesapp.jplace_utils import demultiplex_pqueries, jplace_parser
        jplace_dat = jplace_parser(self.test_jplace)
        pqueries = demultiplex_pqueries(jplace_dat)
        for pquery in pqueries:
            self.assertIsInstance(pquery.placements, list)
        self.assertEqual(['edge_num', 'likelihood', 'like_weight_ratio', 'distal_length', 'pendant_length'],
                         pqueries[0].fields)
        return


if __name__ == '__main__':
    unittest.main()

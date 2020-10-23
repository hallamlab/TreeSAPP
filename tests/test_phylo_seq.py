import unittest
import pytest


class PhyloSeqtests(unittest.TestCase):
    def setUp(self) -> None:
        self.placement_dict = {'p': [[0, -56570.791, 1.0, 0.859, 1.227]],
                               'n': ['scaffold_52752_c1_1 partial=11;start_type=Edge;rbs_motif=None|McrA|1_204']}
        return

    def test_phylo_place(self):
        from treesapp.phylo_seq import PhyloPlace
        bad_dict = {}
        for k, v in self.placement_dict.items():
            bad_dict[k] = v.copy()
        bad_dict['n'].append("second seq name")
        with pytest.raises(SystemExit):
            PhyloPlace(bad_dict)
        pplace = PhyloPlace(self.placement_dict)
        self.assertEqual(pplace.like_weight_ratio, 1.0)
        self.assertEqual(pplace.distal_length, 0.859)
        return


if __name__ == '__main__':
    unittest.main()

import unittest


class PlacementTrainerTestCase(unittest.TestCase):
    def test_reduce_examples(self):
        from treesapp.placement_trainer import reduce_examples
        test_size = 10  # Maximum is 23 (number of seq names in the dict)
        candidates_dict = {'p': {'bacp1': ['seq1', 'seq2', 'seq3', 'seq4'],
                                 'bacp2': ['seq5', 'seq6', 'seq7'],
                                 'arcp1': ['seq8', 'seq9', 'seq10']},
                           'o': {'bacc1': ['seq1', 'seq2'],
                                 'bacc2': ['seq3', 'seq4'],
                                 'bacc3': ['seq5', 'seq6', 'seq7'],
                                 'arcc1': ['seq8', 'seq9'],
                                 'arcc2': ['seq10']},
                           's': {'bacs1': ['seq1'],
                                 'bacs2': ['seq3'],
                                 'arcs1': ['seq9']}}
        trimmed_candidates = reduce_examples(candidate_seqs=candidates_dict, max_examples=test_size)
        i = 0
        for r, t in trimmed_candidates.items():
            i += sum([len(n) for n in t.values()])
        self.assertEqual(test_size, i)


if __name__ == '__main__':
    unittest.main()

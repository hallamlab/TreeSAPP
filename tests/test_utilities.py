import unittest


class ClassifierTester(unittest.TestCase):
    def test_median(self):
        from treesapp.utilities import median
        l1 = [1, 2, 3, 4, 5]
        self.assertEqual(3, median(l1))
        l2 = [5.2]
        self.assertEqual(5.2, median(l2))
        self.assertEqual(None, median([]))
        return

    def test_mean(self):
        from treesapp.utilities import mean
        self.assertEqual(4.5, mean([4, 5]))
        self.assertEqual(4, mean([4]))
        return


if __name__ == '__main__':
    unittest.main()

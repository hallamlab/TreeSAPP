import unittest
import pytest
import re
import os

from .testing_utils import get_test_data

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

    def test_get_file_lines(self):
        from .testing_utils import get_test_data
        from treesapp.utilities import get_file_lines
        self.assertEqual(12, len(get_file_lines(file_path=get_test_data("colours_file.txt"))))
        self.assertEqual(110, len(get_file_lines(file_path=get_test_data("create_test.faa"),
                                                 re_pattern=re.compile("^>.*"))))
        return

    def test_concatenate_files(self):
        from treesapp import utilities
        # Fail if output can't be opened
        with pytest.raises(SystemExit):
            utilities.concatenate_files([get_test_data("create_test.faa"),
                                         get_test_data("fasta_parser_test.fasta")],
                                        "bad/path/cat_test.fasta")

        utilities.concatenate_files([get_test_data("create_test.faa"),
                                     get_test_data("fasta_parser_test.fasta")],
                                    "cat_test.fasta")
        self.assertTrue(os.path.isfile("cat_test.fasta"))
        # Ensure the number of lines is expected
        c = 0
        with open("cat_test.fasta", 'r') as out:
            for _ in out:
                c += 1
        self.assertEqual(266, c)
        os.remove("cat_test.fasta")
        return


if __name__ == '__main__':
    unittest.main()

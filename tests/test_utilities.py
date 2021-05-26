import unittest
import pytest
import shutil
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

    def test_find_msa_type(self):
        from treesapp.utilities import find_msa_type
        # Test failure with an unrecognized file extension
        with pytest.raises(SystemExit):
            find_msa_type({"PuhA": ["mock_msa.stk"]})
        # Test failure with multiple different file extensions
        with pytest.raises(SystemExit):
            find_msa_type({"PuhA": ["f1.mfa"],
                           "TyrA": ["f2.sto"]})

        # Test success
        self.assertEqual("Phylip", find_msa_type({"PuhA": ["f1.phy"]}))
        return

    def test_match_file(self):
        from treesapp.utilities import match_file
        test_dir = "./tests/test_data/s__[Mycobacterium]_[s]tephanolepidis/"
        test_file_one = test_dir + "tree_data_evaluate.raxml.bestModel"

        os.makedirs(test_dir, exist_ok=True)
        tmp_fh = open(test_file_one, 'w')
        tmp_fh.close()
        tmp_path = match_file(test_file_one)
        self.assertTrue(os.path.exists(tmp_path))

        if os.path.isdir(test_dir):
            shutil.rmtree(test_dir)
        return

    def test_complement_nucs(self):
        from treesapp.utilities import complement_nucs
        self.assertEqual("AAAA", complement_nucs("TTTT"))
        self.assertEqual("ANNN", complement_nucs("UXBY"))

        return

    def test_reverse_complement(self):
        from treesapp.utilities import reverse_complement
        self.assertEqual("AAAA", reverse_complement("TTTT"))
        self.assertEqual("TCGA", reverse_complement("TCGA"))
        self.assertEqual("GTCC", reverse_complement("GGAC"))
        return

    def test_elegant_pair(self):
        from treesapp.utilities import elegant_pair
        self.assertEqual(elegant_pair(138, 139, True),
                         elegant_pair(139, 138, True))
        self.assertEqual(17, elegant_pair(1, 4))
        self.assertEqual(21, elegant_pair(4, 1))
        return


if __name__ == '__main__':
    unittest.main()

import unittest
import pytest
import shutil
import os

from .testing_utils import get_test_data


class MccTester(unittest.TestCase):
    def setUp(self) -> None:
        self.num_procs = 4
        self.refpkg_dir = get_test_data(os.path.join("refpkgs"))
        self.ts_test_dir = "./MCC_TreeSAPP"
        self.gm_test_dir = "./MCC_GraftM"

        if not os.path.isdir(self.gm_test_dir):
            os.mkdir(self.gm_test_dir)
        return

    def tearDown(self) -> None:
        for dir_name in [self.ts_test_dir, self.gm_test_dir]:
            if os.path.isdir(dir_name):
                shutil.rmtree(dir_name)
        return

    def test_mcc_calculator(self):
        from treesapp import mcc_calculator
        from treesapp.fasta import read_fasta_to_dict
        from .testing_utils import get_test_data
        test_fa = get_test_data("EggNOG_McrA.faa")
        annotations_map = get_test_data("EggNOG_McrA_annot_map.tsv")

        # Ensure the number of sequences in test fasta and annotations map are as expected
        self.assertEqual(54, len(read_fasta_to_dict(test_fa)))
        with open(annotations_map) as annotations:
            self.assertEqual(52, len(annotations.readlines()))

        cmd = ["--fastx_input", test_fa,
               "--annot_map", annotations_map,
               "--output", self.ts_test_dir,
               "--refpkg_dir", self.refpkg_dir,
               "--targets", "McrA",
               "--molecule", "prot",
               "--tool", "treesapp",
               "--num_procs", str(self.num_procs),
               "--delete", "--svm", "--overwrite"]
        mcc_calculator.mcc_calculator(cmd)

        # Compare the number of MCC classifications to number of TP and FN
        with open(os.path.join(self.ts_test_dir, "EggNOG_McrA_MCC_table.tsv")) as mcc_tbl:
            tbl_lines = mcc_tbl.readlines()
            tbl_lines.pop(0)
            cls_counts = [line.split()[2:] for line in tbl_lines]
            for dist in cls_counts:
                tp, tn, fp, fn = [int(x) for x in dist]
                self.assertTrue(fp >= 2)
                self.assertTrue(fn == 0)
                self.assertEqual(54, sum([tp, tn, fp]))

        # Read MCC_table.tsv to ensure the correct number of sequences were classified
        n_lines = 0
        tp = 0
        fn = 0
        with open(os.path.join(self.ts_test_dir, "EggNOG_McrA_classifications.tsv")) as mcc_classified:
            mcc_classified.readline()
            line = mcc_classified.readline()
            while line:
                if line.split()[3] == "True":
                    tp += 1
                else:
                    fn += 1
                line = mcc_classified.readline()
                n_lines += 1
        self.assertEqual(52, n_lines)
        self.assertEqual(52, tp)
        self.assertEqual(0, fn)

        ##
        # Test MCC calculator with GraftM
        ##
        cmd = ["--fastx_input", test_fa,
               "--annot_map", annotations_map,
               "--tool", "graftm",
               "--output", self.gm_test_dir,
               "--targets", "McrA",
               "--molecule", "prot",
               "--num_procs", str(self.num_procs),
               "--overwrite"]
        with pytest.raises(SystemExit):
            mcc_calculator.mcc_calculator(cmd)
        return


if __name__ == '__main__':
    unittest.main()

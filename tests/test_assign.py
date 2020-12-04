import os
import shutil
import unittest
import pytest

from collections import namedtuple
from .testing_utils import get_test_data


class AssignerTester(unittest.TestCase):
    def setUp(self) -> None:
        self.output_dir = "./tests/assigner_test/"
        self.mock_args = namedtuple('args', ['reads', 'reverse', 'paired', 'rpkm',
                                             'molecule', 'overwrite', 'stage', 'output'])
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)

    def tearDown(self) -> None:
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

    def test_decide_stage(self):
        from treesapp import assign
        ts_assigner = assign.Assigner()
        ts_assigner.output_dir = self.output_dir
        ts_assigner.var_output_dir = self.output_dir
        ts_assigner.final_output_dir = self.output_dir
        # Test normal condition for protein
        ts_assigner.decide_stage(self.mock_args(molecule="prot", overwrite=False, stage="continue", output=self.output_dir,
                                                reads="", reverse="", rpkm=False, paired="pe"))
        self.assertTrue(ts_assigner.stage_status("clean"))
        self.assertFalse(ts_assigner.stage_status("orf-call"))
        self.assertFalse(ts_assigner.stage_status("abundance"))

        # Test a nucleotide input but with no rpkm
        ts_assigner.change_stage_status("orf-call", True)
        ts_assigner.change_stage_status("abundance", True)
        ts_assigner.decide_stage(self.mock_args(molecule="dna", overwrite=False, stage="continue", output=self.output_dir,
                                                reads="", reverse="", rpkm=False, paired="pe"))
        self.assertTrue(ts_assigner.stage_status("orf-call"))
        self.assertFalse(ts_assigner.stage_status("abundance"))

        # Test an RPKM run with bad fastq files
        with pytest.raises(SystemExit):
            ts_assigner.decide_stage(
                self.mock_args(molecule="dna", overwrite=False, stage="continue", output=self.output_dir,
                               reads="test_TarA.1.fq", reverse="", rpkm=True, paired="pe"))

        # Test a successful RPKM run
        ts_assigner.decide_stage(self.mock_args(molecule="dna", overwrite=False,
                                                stage="continue", output=self.output_dir,
                                                reads=get_test_data("test_TarA.1.fq"), reverse=get_test_data("test_TarA.2.fq"),
                                                rpkm=True, paired="pe"))
        self.assertTrue(ts_assigner.stage_status("abundance"))
        return


if __name__ == '__main__':
    unittest.main()

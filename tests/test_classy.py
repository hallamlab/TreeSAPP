import unittest
import pytest
import os
from collections import namedtuple


@pytest.fixture(scope="class")
def treesapp_instance(request):
    from treesapp.classy import TreeSAPP
    request.cls.db = TreeSAPP("")
    request.cls.db.output_dir = "./temp_tmp"
    return


@pytest.mark.usefixtures("treesapp_instance")
class TreeSAPPClassTester(unittest.TestCase):
    def test_furnish_with_arguments(self):
        from .testing_utils import get_test_data
        args = namedtuple("args", ["molecule", "output", "input", "executables"])
        args.input = get_test_data("marker_test.faa")
        args.output = self.db.output_dir
        args.molecule = "prot"
        args.executables = {'prodigal': '/home/connor/bin/prodigal', 'BMGE.jar': '/usr/local/bin/BMGE.jar',
                            'hmmbuild': '/usr/local/bin/hmmbuild',
                            'hmmalign': '/usr/local/bin/hmmalign',
                            'hmmsearch': '/usr/local/bin/hmmsearch',
                            'epa-ng': '/usr/local/bin/epa-ng', 'raxml-ng': '/usr/local/bin/raxml-ng'}
        self.db.furnish_with_arguments(args)
        self.assertEqual(len(args.executables), len(self.db.executables))
        self.assertEqual(get_test_data("marker_test.faa"), self.db.input_sequences)
        self.assertEqual(args.molecule, self.db.molecule_type)
        return

    def test_check_previous_output(self):
        self.db.final_output_dir = self.db.output_dir + os.sep + "final_outputs"
        self.db.var_output_dir = self.db.output_dir + os.sep + "intermediates"
        self.db.check_previous_output(overwrite=True)
        self.assertTrue(os.path.isdir(self.db.output_dir))
        return

    def test_validate_new_dir(self):
        from treesapp.utilities import validate_new_dir
        ts_inst = self.db
        ts_inst.output_dir = validate_new_dir(ts_inst.output_dir)
        self.assertEqual(os.sep, ts_inst.output_dir[-1])
        with pytest.raises(SystemExit):
            validate_new_dir(os.path.join(ts_inst.output_dir, "nope", "nope"))
        return


if __name__ == '__main__':
    unittest.main()

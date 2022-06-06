import os
import unittest

from .testing_utils import get_test_data


class MyTestCase(unittest.TestCase):
    def test_trim_multiple_alignments(self):
        from treesapp import multiple_alignment
        from treesapp import refpkg
        test_fa = get_test_data('PuhA.mfa')
        output_file = os.path.join("tests", "test_data", "PuhA.trim.mfa")
        test_rp = refpkg.ReferencePackage(refpkg_name="PuhA")
        test_rp.f__pkl = get_test_data(filename=os.path.join("refpkgs", "PuhA_build.pkl"))
        test_rp.slurp()

        result = multiple_alignment.trim_multiple_alignment_farmer({"PuhA": [test_fa]},
                                                                   min_seq_length=10,
                                                                   n_proc=1,
                                                                   ref_pkgs={"PuhA": test_rp})
        self.assertTrue(os.path.isfile(output_file))
        self.assertIsInstance(result, dict)
        self.assertTrue("PuhA" in result.keys())
        self.assertEqual(os.path.basename(output_file),
                         os.path.basename(result["PuhA"].pop()))

        if os.path.isfile(output_file):
            os.remove(output_file)
        return


if __name__ == '__main__':
    unittest.main()

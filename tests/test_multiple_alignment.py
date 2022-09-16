import os
import unittest

from .testing_utils import get_test_data


class MyTestCase(unittest.TestCase):
    def test_trim_multiple_alignments(self):
        from treesapp import multiple_alignment
        from treesapp import refpkg
        test_fa = get_test_data('PuhA.mfa')
        trim_file = os.path.join("tests", "test_data", "PuhA.trim.mfa")
        qc_file = os.path.join("tests", "test_data", "PuhA.trim.qc.mfa")
        test_rp = refpkg.ReferencePackage(refpkg_name="PuhA")
        test_rp.f__pkl = get_test_data(filename=os.path.join("refpkgs", "PuhA_build.pkl"))
        test_rp.slurp()

        result = multiple_alignment.trim_multiple_alignment_farmer([{"qry_ref_mfa": test_fa,
                                                                     "refpkg_name": "PuhA",
                                                                     "gap_tuned": True,
                                                                     "avg_id": 88}],
                                                                   min_seq_length=10,
                                                                   n_proc=1,
                                                                   ref_pkgs={"PuhA": test_rp},
                                                                   for_placement=False)
        self.assertTrue(os.path.isfile(trim_file))
        self.assertIsInstance(result, dict)
        self.assertTrue("PuhA" in result.keys())
        self.assertEqual(os.path.basename(qc_file),
                         os.path.basename(result["PuhA"].pop()))

        for f_path in [trim_file, qc_file]:
            if os.path.isfile(f_path):
                os.remove(f_path)
        return


if __name__ == '__main__':
    unittest.main()

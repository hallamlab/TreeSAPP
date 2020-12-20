import os
import unittest

from .testing_utils import get_test_data


class UpdaterTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp import refpkg as ts_rp
        self.ref_pkg = ts_rp.ReferencePackage()
        self.ref_pkg.f__json = get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.ref_pkg.slurp()
        return

    def test_decide_length_filter(self):
        from treesapp import update_refpkg as ts_update_mod
        min_length = ts_update_mod.decide_length_filter(self.ref_pkg, proposed_min_length=800)
        self.assertEqual(800, min_length)

        # Ensure the incorrectly-provided percentage is converted to a proportion
        min_length = ts_update_mod.decide_length_filter(self.ref_pkg, min_hmm_proportion=80)
        self.assertEqual(443, min_length)
        return


if __name__ == '__main__':
    unittest.main()

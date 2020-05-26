import unittest


class TreesappTester(unittest.TestCase):
    def test_gather_ref_packages(self):
        from treesapp.file_parsers import gather_ref_packages
        import os
        from . import testing_utils as utils
        refpkg_dir = utils.get_test_data(os.path.join("refpkgs"))
        refpkg_dict = gather_ref_packages(refpkg_data_dir=refpkg_dir, targets=["McrA"])
        self.assertEqual(1, len(refpkg_dict))
        return


if __name__ == '__main__':
    unittest.main()

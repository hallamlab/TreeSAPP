import unittest
import pytest


@pytest.fixture(scope="class")
def refpkg_class(request):
    import os
    from treesapp import refpkg
    from . import testing_utils as utils
    request.cls.db = refpkg.ReferencePackage("McrA")
    request.cls.db.f__json = "./band_test.json"
    request.cls.db.gather_package_files(pkg_path=os.path.join(utils.get_treesapp_path(), "treesapp", "data"))
    return


@pytest.mark.usefixtures("refpkg_class")
class RefPkgTester(unittest.TestCase):
    def test_band(self):
        self.db.band()
        return

    def test_slurp(self):
        self.db.slurp()
        self.assertEqual("McrA", self.db.prefix)
        return

    def test_disband(self):
        import os
        self.db.disband("./")
        self.assertTrue(os.path.isfile("./McrA__/McrA.fa"))


if __name__ == '__main__':
    unittest.main()

import unittest
import os
import shutil

from .testing_utils import get_test_data, get_treesapp_root


class CladeExclusionTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.refpkg import ReferencePackage
        from treesapp.utilities import fetch_executable_path

        # Create and slurp the reference package that will be used for clade exclusions
        self.ref_pkg = ReferencePackage()
        self.puha_pkl = get_test_data(os.path.join("refpkgs", "PuhA_build.pkl"))
        self.ref_pkg.f__json = get_test_data(self.puha_pkl)
        self.ref_pkg.slurp()

        # Create the directory for temporary files
        self.tmp_dir = "./ce_intermediates"
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.mkdir(self.tmp_dir)

        # Find the executables
        self.exe_map = {}
        self.treesapp_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
        for dep in ["hmmbuild", "hmmalign", "raxml-ng", "mafft", "BMGE.jar"]:
            self.exe_map[dep] = fetch_executable_path(dep, get_treesapp_root())

        self.num_procs = 2
        return

    def tearDown(self) -> None:
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        return

    def test_prep_graftm_ref_files(self):
        from treesapp.clade_exclusion_evaluator import prep_graftm_ref_files
        taxon_str = 'd__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales'
        ce_dir = os.path.join(self.tmp_dir, taxon_str.split("; ")[-1])
        prep_graftm_ref_files(tmp_dir=ce_dir, target_clade=taxon_str,
                              refpkg=self.ref_pkg, executables=self.exe_map)
        self.assertTrue(os.path.exists(os.path.join(ce_dir, "PuhA_lineage_ids.txt")))
        self.assertTrue(os.path.exists(os.path.join(ce_dir, "PuhA.fa")))
        self.assertTrue(os.path.exists(os.path.join(ce_dir, "PuhA.mfa")))
        return

    # def test_build_graftm_package(self):
    #     from treesapp.clade_exclusion_evaluator import build_graftm_package
    #     pkg_path = os.path.join(self.tmp_dir, "PuhA.gpkg")
    #     build_graftm_package(gpkg_path=pkg_path, threads=self.num_procs,
    #                          tax_file=get_test_data("PuhA_lineage_ids.txt"),
    #                          fa_file=get_test_data("PuhA.fa"),
    #                          mfa_file=get_test_data("PuhA.mfa"))
    #     self.assertTrue(os.path.isdir(pkg_path))
    #     return


if __name__ == '__main__':
    unittest.main()

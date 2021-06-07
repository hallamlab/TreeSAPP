import os
import shutil
import unittest

from .testing_utils import get_test_data, get_treesapp_root


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp import refpkg
        from treesapp import utilities
        self.num_procs = 4
        self.ts_dir = get_treesapp_root()
        self.tmp_dir = os.path.join(os.getcwd(), "tests", "graftm_test_dir") + os.sep
        if not os.path.isdir(self.tmp_dir):
            os.mkdir(self.tmp_dir)

        self.ref_pkg = refpkg.ReferencePackage()
        self.ref_pkg.f__pkl = get_test_data(get_test_data(os.path.join("refpkgs", "PuhA_build.pkl")))
        self.ref_pkg.slurp()
        try:
            self.graftm_exe = utilities.fetch_executable_path("graftM", self.ts_dir)
        except SystemExit:
            self.graftm_exe = None
        return

    def tearDown(self) -> None:
        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        return

    def test_read_graftm_classifications(self):
        from treesapp import graftm_utils
        from treesapp import phylo_seq
        assignments = graftm_utils.read_graftm_classifications(get_test_data("graftm_raw_read_tax.tsv"))
        self.assertIsInstance(assignments, list)
        self.assertEqual(9, len(assignments))
        self.assertIsInstance(assignments[0], phylo_seq.PQuery)
        return

    def test_grab_graftm_taxa(self):
        from treesapp import graftm_utils
        from pygtrie import StringTrie
        tax_trie = graftm_utils.grab_graftm_taxa(tax_ids_file=get_test_data("graftm_taxonomy.csv"))
        self.assertIsInstance(tax_trie, StringTrie)
        self.assertEqual(211, len(tax_trie))
        return

    def test_run_graftm_graft(self):
        from treesapp import graftm_utils as g_utils
        if not self.graftm_exe:
            return
        classification_tbl = os.path.join(self.tmp_dir, "EggNOG_McrA", "EggNOG_McrA_read_tax.tsv")

        # Run graftM graft
        g_utils.run_graftm_graft(graftm_exe=self.graftm_exe,
                                 input_path=get_test_data("EggNOG_McrA.faa"),
                                 output_dir=self.tmp_dir,
                                 gpkg_path=get_test_data(os.path.join("refpkgs", "McrA.gpkg")),
                                 classifier="graftm",
                                 num_threads=self.num_procs)

        self.assertTrue(os.path.isfile(classification_tbl))
        with open(classification_tbl) as tbl_handler:
            tbl_lines = tbl_handler.readlines()
        self.assertEqual(54, len(tbl_lines))
        return

    def test_prep_graftm_ref_files(self):
        from treesapp import graftm_utils
        from treesapp import utilities
        # Find the executables
        exe_map = {}
        for dep in ["hmmbuild", "hmmalign", "raxml-ng", "mafft", "BMGE.jar"]:
            exe_map[dep] = utilities.fetch_executable_path(dep, self.ts_dir)

        taxon_str = 'd__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales'
        ce_dir = os.path.join(self.tmp_dir, taxon_str.split("; ")[-1])
        graftm_utils.prep_graftm_ref_files(tmp_dir=ce_dir, target_clade=taxon_str,
                                           ref_pkg=self.ref_pkg, executables=exe_map)
        self.assertTrue(os.path.exists(os.path.join(ce_dir, "PuhA_lineage_ids.txt")))
        self.assertTrue(os.path.exists(os.path.join(ce_dir, "PuhA.fa")))
        self.assertTrue(os.path.exists(os.path.join(ce_dir, "PuhA.mfa")))
        return

    def test_build_graftm_package(self):
        from treesapp import graftm_utils
        if not self.graftm_exe:
            return
        pkg_path = os.path.join(self.tmp_dir, "PuhA.gpkg")
        graftm_utils.build_graftm_package(gpkg_path=pkg_path, threads=self.num_procs,
                                          tax_file=get_test_data("PuhA_lineage_ids.txt"),
                                          fa_file=get_test_data("PuhA.fa"),
                                          mfa_file=get_test_data("PuhA.mfa"))
        self.assertTrue(os.path.isdir(pkg_path))
        return


if __name__ == '__main__':
    unittest.main()

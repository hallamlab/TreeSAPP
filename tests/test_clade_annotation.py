import os
import unittest

from .testing_utils import get_test_data


class TestCladeAnnotation(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.clade_annotation import CladeAnnotation
        from treesapp.refpkg import ReferencePackage

        self.ref_pkg = ReferencePackage()
        self.ref_pkg.f__pkl = get_test_data(os.path.join("refpkgs", "XmoA_build.pkl"))
        self.ref_pkg.slurp()

        self.test_clade_annot = CladeAnnotation(name="test_feature", key="Paralog")
        self.test_clade_annot.members = {'5_XmoA', '2_XmoA', '8_XmoA', '3_XmoA', '10_XmoA', '7_XmoA', '1_XmoA', '4_XmoA', '9_XmoA', '6_XmoA'}
        return

    def test_get_internal_nodes(self):
        internal_node_map = self.ref_pkg.get_internal_node_leaf_map()
        internal_nodes = self.test_clade_annot.get_internal_nodes(internal_node_map)
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 185], internal_nodes)
        return

    def test_summarise(self):
        summary = self.test_clade_annot.summarise()
        self.assertIsInstance(summary, str)
        self.assertTrue(len(summary) > 0)
        return


if __name__ == '__main__':
    unittest.main()

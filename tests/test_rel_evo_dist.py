import os
import unittest

from .testing_utils import get_test_data


class RelDistTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.refpkg import ReferencePackage
        from ete3 import Tree
        self.test_refpkg = ReferencePackage()
        self.test_refpkg.f__json = get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.test_refpkg.slurp()

        self.label_tree = self.test_refpkg.taxonomically_label_tree()

        simple_nwk = "((D:0.7,F:0.5)E:0.06,(B:0.2,H:0.8)B:0.8);"
        self.simple_tree = Tree(simple_nwk, format=1)
        return

    def test__avg_descendant_rate(self):
        from treesapp.rel_evo_dist import RedTree
        red = RedTree()
        red._avg_descendant_rate(self.simple_tree)
        mrca = self.simple_tree.get_tree_root()
        self.assertEqual(4, mrca.num_taxa)
        self.assertEqual(1.53, mrca.mean_dist)
        for l in mrca.get_leaves():
            self.assertEqual(0, l.mean_dist)
        return

    def test_decorate_rel_dist(self):
        from treesapp.rel_evo_dist import RedTree
        red = RedTree()
        red.decorate_rel_dist(self.simple_tree)
        self.assertEqual(0.0, self.simple_tree.rel_dist)
        return

    def test_rel_dist_to_named_clades(self):
        from treesapp.rel_evo_dist import RedTree
        red = RedTree()
        rel_dists = red.rel_dist_to_named_clades(self.label_tree)

        self.assertTrue(rel_dists['root']['Root'] > rel_dists['domain']['Archaea'])
        self.assertEqual(1, len(rel_dists['domain']))
        self.assertEqual(0,
                         len({'root', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'}.difference(
                             set(rel_dists.keys())))
                         )
        return


if __name__ == '__main__':
    unittest.main()

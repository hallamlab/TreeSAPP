import unittest
import pytest
import os

from ete3 import Tree, TreeNode

from .testing_utils import get_test_data


@pytest.fixture(scope="class")
def refpkg_class(request):
    from treesapp import refpkg
    from . import testing_utils as utils
    request.cls.db = refpkg.ReferencePackage("McrA")
    request.cls.db.f__json = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
    request.cls.db.slurp()
    return


@pytest.mark.usefixtures("refpkg_class")
class EntishTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.entish import load_ete3_tree, label_internal_nodes_ete
        with open(get_test_data("jplace.tree")) as test_tree:
            self.placement_tree = test_tree.readline()
        self.test_ete_tree = load_ete3_tree(self.placement_tree)
        self.mock_tree = Tree("(A:1,(B:1,(E:1,D:1):0.5):0.5);")
        self.multifurcating_tree_str = "(B,(A, C, D), E);"
        label_internal_nodes_ete(self.test_ete_tree)
        self.refpkg_tree = self.db.tree
        return

    def test_get_ete_edge(self):
        from treesapp.entish import get_ete_edge
        # Test an interesting edge
        p, c = get_ete_edge(self.test_ete_tree, 490)  # type: (TreeNode, TreeNode)
        self.assertEqual('490', c.name)
        self.assertIsNone(p)
        self.assertEqual(self.test_ete_tree.get_leaves(), c.get_leaves())

        # Test an edge that doesn't exist
        self.assertEqual(None, get_ete_edge(self.test_ete_tree, 492))

        # Test a leaf node
        p, c = get_ete_edge(self.mock_tree, 1)
        self.assertEqual(0, len(c.children))

        # Test root edge
        p, c = get_ete_edge(self.mock_tree, 6)
        self.assertEqual(4, len(c.get_leaves()))
        self.assertEqual(None, p)
        return

    def test_label_internal_nodes_ete(self):
        from treesapp.entish import label_internal_nodes_ete
        label_internal_nodes_ete(self.mock_tree)
        multi_tree = Tree(self.multifurcating_tree_str)
        for n in self.mock_tree:
            self.assertIsInstance(n.name, str)
        self.assertEqual(7, len(multi_tree.get_edges()))
        label_internal_nodes_ete(multi_tree)
        self.assertEqual(9, len(multi_tree.get_edges()))
        self.assertEqual('8', multi_tree.get_tree_root().name)
        return

    def test_edge_from_node_name(self):
        from treesapp.entish import edge_from_node_name
        edge_name = edge_from_node_name(self.mock_tree, 'D')
        self.assertEqual(3, edge_name)
        return

    def test_map_internal_nodes_leaves(self):
        from treesapp.entish import map_internal_nodes_leaves
        node_map = map_internal_nodes_leaves(self.placement_tree)
        self.assertEqual(489, len(node_map))
        return

    def test_verify_bifurcations(self):
        from treesapp.entish import verify_bifurcations
        pre_edges = len(Tree(self.multifurcating_tree_str).get_edges())
        temp_tree = verify_bifurcations(self.multifurcating_tree_str)
        post_edges = len(Tree(temp_tree).get_edges())
        self.assertTrue(post_edges > pre_edges)
        self.assertEqual(len(Tree(self.multifurcating_tree_str).get_leaves()),
                         len(Tree(temp_tree).get_leaves()))

        # Test for different formats with internal node identifiers
        rp_tree = verify_bifurcations(self.refpkg_tree)
        self.assertEqual(246, len(Tree(rp_tree).get_leaves()))
        self.assertTrue(len(Tree(rp_tree).get_edges()) >= len(Tree(self.refpkg_tree).get_edges()))

        rf = Tree(rp_tree).robinson_foulds(Tree(self.refpkg_tree), unrooted_trees=True)
        self.assertEqual(0, rf[0])
        return


if __name__ == "__main__":
    unittest.main()

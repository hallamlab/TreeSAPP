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
    request.cls.db.f__pkl = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
    request.cls.db.slurp()
    return


@pytest.mark.usefixtures("refpkg_class")
class EntishTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.entish import load_ete3_tree, label_internal_nodes_ete
        with open(get_test_data("jplace.tree")) as test_tree:
            self.placement_tree = test_tree.readline()
        self.test_ete_tree = load_ete3_tree(self.placement_tree)
        self.mock_tree_str = "(A:1,(B:1,(E:1,D:1):0.5):0.5);"
        self.mock_tree = Tree(self.mock_tree_str)
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
        label_internal_nodes_ete(self.mock_tree, attr_type=str)
        multi_tree = Tree(self.multifurcating_tree_str)
        for n in self.mock_tree:
            self.assertIsInstance(n.name, str)
        self.assertEqual(7, len(multi_tree.get_edges()))
        label_internal_nodes_ete(multi_tree)
        self.assertEqual(9, len(multi_tree.get_edges()))
        self.assertEqual('8', multi_tree.get_tree_root().name)

        label_tree = Tree(self.mock_tree_str)
        label_internal_nodes_ete(ete_tree=label_tree, attr_type=int, attr="i_node")
        for n in label_tree:
            self.assertIsInstance(n.i_node, int)
        return

    def test_edge_from_node_name(self):
        from treesapp import entish
        test_tree = self.mock_tree.copy()
        entish.label_internal_nodes_ete(test_tree)
        edge_name = entish.edge_from_node_name(test_tree, 'D')
        self.assertEqual(3, edge_name)
        edge_name = entish.edge_from_node_name(test_tree, 4)
        self.assertEqual(4, edge_name)
        return

    def test_map_internal_nodes_leaves(self):
        from treesapp.entish import map_internal_nodes_leaves
        node_map = map_internal_nodes_leaves(self.placement_tree)
        self.assertEqual(491, len(node_map))
        # Ensure the root node is in the node map by ensuring a node maps to all leaves
        self.assertEqual(len(self.test_ete_tree.get_leaves()),
                         max([len(node_map[x]) for x in node_map]))
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
        self.assertEqual(251, len(Tree(rp_tree).get_leaves()))
        self.assertTrue(len(Tree(rp_tree).get_edges()) >= len(Tree(self.refpkg_tree).get_edges()))

        rf = Tree(rp_tree).robinson_foulds(Tree(self.refpkg_tree), unrooted_trees=True)
        self.assertEqual(0, rf[0])
        return

    def test_match_leaves_to_internal_nodes(self):
        from treesapp.entish import match_leaves_to_internal_nodes
        test_leaves = ["1_R", "2_R", "3_R", "8_R"]
        test_inodes = {0: ["1_R"], 1: ["2_R"], 3: ["1_R", "2_R"],
                       4: ["3_R"], 5: ["4_R"], 6: ["3_R", "4_R"], 7: ["1_R", "2_R"] + ["3_R", "4_R"],
                       8: ["8_R"], 9: ["1_R", "2_R"] + ["3_R", "4_R"] + ["8_R"]}
        i_nodes = match_leaves_to_internal_nodes(leaf_names=test_leaves, internal_node_leaf_map=test_inodes)
        self.assertEqual([3, 4, 8], sorted(i_nodes, key=int))
        self.assertEqual(4, len(test_leaves))
        self.assertEqual(9, len(test_inodes))
        return

    def test_collapse_ete_tree(self):
        from treesapp.entish import collapse_ete_tree
        mock_tree = Tree("(A:1,(B:1,(E:0.1,D:0.4):0.1):0.5);")
        collapse_ete_tree(mock_tree, min_branch_length=0.6)
        self.assertEqual(['A', 'B'], mock_tree.get_leaf_names())
        self.assertEqual(3, len(mock_tree.get_edges()))

        rp_tree = Tree('(((((((((217_McrA:1.40284)1:0.072783)1:0.076191)1:0.060539)1:0.048257,488:0.039767)1:0.07175)'
                       '1:0.109454)1:0.036157)1:0.0282665,(129:0.064633,(((67:0.055106)1:0.086543)1:0.054828,'
                       '((201_McrA:1.01787,(203_McrA:0.361599)1:0.213719)1:0.443498,(((200_McrA:0.412511)1:0.277927)'
                       '1:0.119218,197_McrA:0.955209)1:0.276102)1:0.121564)1:0.035701)1:0.0282665)1:0.107436;')
        collapse_ete_tree(rp_tree, min_branch_length=0.67)
        self.assertEqual(3, len(rp_tree))

        rp_tree = Tree('((49_McrA:0.035,((37_McrA:0.029,9_McrA:0.025)1:0.014,'
                       '(7_McrA:0.031,194_McrA:0.031)1:0.019)1:0.02)1:0.015,'
                       '(6_McrA:0.037,((191_McrA:0.031,192_McrA:0.051)1:0.018,27_McrA:0.079)1:0.022)1:0.015);')
        collapse_ete_tree(rp_tree, min_branch_length=0.0375)
        self.assertEqual(4, len(rp_tree))
        self.assertEqual(7, len(rp_tree.get_edges()))
        return


if __name__ == "__main__":
    unittest.main()

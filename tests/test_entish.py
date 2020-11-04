import unittest

from ete3 import Tree, TreeNode

from .testing_utils import get_test_data


class EntishTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.entish import load_ete3_tree, label_internal_nodes_ete
        with open(get_test_data("jplace.tree")) as test_tree:
            self.placement_tree = test_tree.readline()
        self.test_ete_tree = load_ete3_tree(self.placement_tree)
        self.mock_tree = Tree("(A:1,(B:1,(E:1,D:1):0.5):0.5);")
        self.multifurcating_tree = Tree("(B,(A, C, D), E);")
        label_internal_nodes_ete(self.test_ete_tree)
        return

    def test_get_ete_edge(self):
        from treesapp.entish import get_ete_edge
        # Test an interesting edge
        p, c = get_ete_edge(self.test_ete_tree, 31)  # type: (TreeNode, TreeNode)
        self.assertEqual('490', p.name)
        self.assertEqual(self.test_ete_tree.get_leaves(), p.get_leaves())

        # Test an edge that doesn't exist
        self.assertEqual(None, get_ete_edge(self.test_ete_tree, 492))

        # Test a leaf node
        p, c = get_ete_edge(self.mock_tree, 1)
        self.assertEqual(0, len(c.children))

        # Test root edge
        p, c = get_ete_edge(self.mock_tree, 7)
        self.assertEqual(4, len(c.get_leaves()))
        self.assertEqual(None, p)
        return

    def test_label_internal_nodes_ete(self):
        from treesapp.entish import label_internal_nodes_ete
        label_internal_nodes_ete(self.mock_tree)
        for n in self.mock_tree:
            self.assertIsInstance(n.name, str)
        self.assertEqual(7, len(self.multifurcating_tree.get_edges()))
        label_internal_nodes_ete(self.multifurcating_tree)
        self.assertEqual(9, len(self.multifurcating_tree.get_edges()))
        self.assertEqual('8', self.multifurcating_tree.get_tree_root().name)
        return

    def test_edge_from_node_name(self):
        from treesapp.entish import edge_from_node_name
        edge_name = edge_from_node_name(self.mock_tree, 'D')
        self.assertEqual(4, edge_name)
        return

    def test_map_internal_nodes_leaves(self):
        from treesapp.entish import map_internal_nodes_leaves
        node_map = map_internal_nodes_leaves(self.placement_tree)
        self.assertEqual(489, len(node_map))
        return


if __name__ == "__main__":
    unittest.main()

import unittest


class EntishTester(unittest.TestCase):
    def test_map_internal_nodes_leaves(self):
        from . import testing_utils as utils
        from treesapp.entish import map_internal_nodes_leaves
        placement_tree = utils.get_test_data("jplace.tree")
        with open(placement_tree) as tree:
            tree_text = tree.readline()
            node_map = map_internal_nodes_leaves(tree_text)
        self.assertEqual(489, len(node_map))
        return


if __name__ == "__main__":
    unittest.main()

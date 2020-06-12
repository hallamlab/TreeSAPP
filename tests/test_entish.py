import pytest


def test_map_internal_nodes_leaves():
    from . import testing_utils as utils
    from treesapp.entish import map_internal_nodes_leaves
    placement_tree = utils.get_test_data("jplace.tree")
    with open(placement_tree) as tree:
        tree_text = tree.readline()
        node_map = map_internal_nodes_leaves(tree_text)
    assert len(node_map) == 489

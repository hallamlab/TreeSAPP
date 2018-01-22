__author__ = 'Connor Morgan-Lang'


from utilities import Autovivify


def get_node(tree, pos):
    node = ""
    pos += 1
    c = tree[pos]
    while c != '}':
        node += c
        pos += 1
        c = tree[pos]
    return int(node), pos


def subtrees_to_dictionary(subtrees_string, tree_info):
    subtree_list = subtrees_string.split(';')
    for subtree in subtree_list:
        node = subtree.split(')')[-1]
        tree_info['subtree_of_node'][node] = subtree
    return tree_info


def create_tree_info_hash():
    tree_info = Autovivify()
    return tree_info

"""Modified from sourmash's sourmash/tests/sourmash_tst_utils.py"""

import os
from pkg_resources import Requirement, resource_filename, ResolutionError

import ete3


def random_ete_tree(leaf_names: list, branch_len_dist=None) -> ete3.Tree:
    if not branch_len_dist:
        branch_len_dist = (0, 1)
    rand_tree = ete3.Tree()
    rand_tree.populate(size=len(leaf_names),
                       names_library=leaf_names,
                       branch_range=branch_len_dist,
                       random_branches=True)
    return rand_tree


def get_treesapp_root():
    return resource_filename(Requirement.parse("treesapp"), 'treesapp')


def get_test_data(filename):
    filepath = None
    try:
        filepath = resource_filename(Requirement.parse("treesapp"), "tests/test_data/" + filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), 'test_data', filename)
    return filepath


def get_treesapp_file(filename):
    return resource_filename(Requirement.parse("treesapp"), filename)


def get_treesapp_path():
    return resource_filename(Requirement.parse("treesapp"), "")

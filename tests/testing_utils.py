"""Modified from sourmash's sourmash/tests/sourmash_tst_utils.py"""

import os
from pkg_resources import Requirement, resource_filename, ResolutionError


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

"""Modified from sourmash's sourmash/tests/sourmash_tst_utils.py"""

import os
from pkg_resources import Requirement, resource_filename, ResolutionError


def get_test_data(filename):
    filepath = None
    try:
        filepath = resource_filename(Requirement.parse("treesapp"), "tests/test_data/" + filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), 'tests/test_data', filename)
    return filepath
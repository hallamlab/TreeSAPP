import unittest
import pytest
import os


@pytest.fixture(scope="class")
def mcra_refpkg_class(request):
    from treesapp import refpkg
    from . import testing_utils as utils
    request.cls.db = refpkg.ReferencePackage("McrA")
    request.cls.db.f__json = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
    request.cls.db.slurp()
    return


@pytest.fixture(scope="class")
def xmoa_refpkg_class(request):
    from treesapp import refpkg
    from . import testing_utils as utils
    request.cls.db = refpkg.ReferencePackage("XmoA")
    request.cls.db.f__json = utils.get_test_data(os.path.join("refpkgs", "XmoA_build.pkl"))
    request.cls.db.slurp()
    return


@pytest.mark.usefixtures()
class MyTestCase(unittest.TestCase):
    @pytest.mark.usefixtures("mcra_refpkg_class")
    def test_names_for_nodes_range(self):
        from .testing_utils import get_test_data
        from treesapp.annotate_extra import names_for_nodes
        from treesapp.file_parsers import read_colours_file
        range_annot_file = get_test_data("McrA_Metabolism.txt")
        groups_dict = read_colours_file(range_annot_file, "McrA")
        clusters = names_for_nodes(groups_dict, self.db.get_internal_node_leaf_map(),
                                   self.db.generate_tree_leaf_references_from_refpkg())
        self.assertEqual(6, len(clusters))
        return

    @pytest.mark.usefixtures("xmoa_refpkg_class")
    def test_convert_outer_to_inner_nodes(self):
        from treesapp.utilities import convert_outer_to_inner_nodes
        group_dict = {'AmoA_AOB': [('40', '40'), ('96', '96')],
                      'PxmA': [('109', '109')],
                      'BmoA': [('148', '148')],
                      'AmoA_AOA': [('141', '141')],
                      'EmoA': [('121', '121')],
                      'PmoA': [('176', '176'), ('73', '73'), ('59', '59'), ('94_XmoA', '94_XmoA')]}
        clusters = convert_outer_to_inner_nodes(internal_node_map=self.db.get_internal_node_leaf_map(),
                                                clusters=group_dict)
        annot_i_nodes = []
        for i_nodes in clusters.values():
            annot_i_nodes += i_nodes
        self.assertEqual(10, len(annot_i_nodes))
        # With variable types representing nodes
        group_dict = {'EmoA': [('69', '69')],
                      'PxmA': [('57', '57')],
                      'AmoA_AOA': [('9', '9'), (0, 0), (10, 10), ('13', '13'), (14, 14)],
                      'AmoA_AOB': [('44', '44'), ('108', '108')],
                      'PmoA': [('175', '175'), ('103', '103'), ('83', '83'), (84, 84)],
                      'BmoA': [('21', '21')]}
        clusters = convert_outer_to_inner_nodes(internal_node_map=self.db.get_internal_node_leaf_map(),
                                                clusters=group_dict)
        annot_i_nodes = []
        for i_nodes in clusters.values():
            annot_i_nodes += i_nodes
        self.assertEqual(14, len(annot_i_nodes))
        return

    @pytest.mark.usefixtures("xmoa_refpkg_class")
    def test_names_for_nodes_inodes(self):
        from .testing_utils import get_test_data
        from treesapp.annotate_extra import names_for_nodes
        from treesapp.file_parsers import read_colours_file
        inode_annot_file = get_test_data("XmoA_Function.txt")
        groups_dict = read_colours_file(inode_annot_file, "XmoA")
        clusters = names_for_nodes(groups_dict, self.db.get_internal_node_leaf_map(),
                                   self.db.generate_tree_leaf_references_from_refpkg())
        self.assertEqual(6, len(clusters))

    def test_read_colours_file(self):
        from treesapp.file_parsers import read_colours_file
        from .testing_utils import get_test_data
        cols_dict = read_colours_file(annotation_file=get_test_data("colours_file.txt"), refpkg_name="McrA")
        self.assertEqual(6, len(cols_dict))

    # def test_annotate_internal_nodes(self):
    #     from annotate_extra import annotate_internal_nodes
    #     return


if __name__ == '__main__':
    unittest.main()

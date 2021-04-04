import unittest
import pytest
import os

from . import testing_utils as utils


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp import refpkg
        from treesapp import annotate_extra
        self.mcra_refpkg = refpkg.ReferencePackage("McrA")
        self.mcra_refpkg.f__pkl = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.mcra_refpkg.slurp()

        self.xmoa_refpkg = refpkg.ReferencePackage("XmoA")
        self.xmoa_refpkg.f__pkl = utils.get_test_data(os.path.join("refpkgs", "XmoA_build.pkl"))
        self.xmoa_refpkg.slurp()

        self.leaf_node_map = self.xmoa_refpkg.get_internal_node_leaf_map()

        # Example classified sequences
        self.cs1 = annotate_extra.ClassifiedSequence("McrA")
        self.cs1.query_name, self.cs1.i_node = 'scaffold_3993_c2_1', '289'
        self.cs2 = annotate_extra.ClassifiedSequence("McrA")
        self.cs2.query_name, self.cs2.i_node = 'scaffold_52752_c1_1', '224'
        self.cs3 = annotate_extra.ClassifiedSequence("McrA")
        self.cs3.query_name, self.cs3.i_node = 'scaffold_52147_c1_2 # 260', '2'

        return

    def test_convert_outer_to_inner_nodes(self):
        from treesapp.entish import convert_outer_to_inner_nodes
        group_dict = {'AmoA_AOB': [('40', '40'), ('96', '96')],
                      'PxmA': [('109', '109')],
                      'BmoA': [('148', '148')],
                      'AmoA_AOA': [('141', '141')],
                      'EmoA': [('121', '121')],
                      'PmoA': [('176', '176'), ('73', '73'), ('59', '59'), ('94_XmoA', '94_XmoA')]}
        clusters = convert_outer_to_inner_nodes(internal_node_map=self.xmoa_refpkg.get_internal_node_leaf_map(),
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
        clusters = convert_outer_to_inner_nodes(internal_node_map=self.xmoa_refpkg.get_internal_node_leaf_map(),
                                                clusters=group_dict)
        annot_i_nodes = []
        for i_nodes in clusters.values():
            annot_i_nodes += i_nodes
        self.assertEqual(14, len(annot_i_nodes))
        return

    def test_read_colours_file(self):
        from treesapp.file_parsers import read_colours_file
        cols_dict = read_colours_file(annotation_file=utils.get_test_data("colours_file.txt"), refpkg_name="McrA")
        self.assertEqual(6, len(cols_dict))

    def test_annotate_internal_nodes(self):
        from treesapp.annotate_extra import annotate_internal_nodes
        from treesapp.clade_annotation import CladeAnnotation
        xmoa_clade_annotations = self.xmoa_refpkg.feature_annotations["Paralog"]
        marker_tree_info, leaves_in_clusters = annotate_internal_nodes(internal_node_map=self.leaf_node_map,
                                                                       clade_annotations=xmoa_clade_annotations)
        self.assertEqual(4, len(marker_tree_info))
        self.assertTrue('PmoA' in marker_tree_info)
        self.assertEqual(16, len(marker_tree_info["PmoA"]))
        self.assertEqual(46, len(leaves_in_clusters))

        # Test an internal node that doesn't exist in the internal node-to-leaf node map
        bad_ca = CladeAnnotation(name="", key="")
        bad_ca.members = {'1', '5', '9', '241'}
        with pytest.raises(SystemExit):
            annotate_internal_nodes(internal_node_map=self.leaf_node_map, clade_annotations=[bad_ca])

        with pytest.raises(AssertionError):
            annotate_internal_nodes(internal_node_map=self.leaf_node_map, clade_annotations=[])
        return

    def test_map_queries_to_annotations(self):
        from treesapp.annotate_extra import map_queries_to_annotations
        # Functional test
        map_queries_to_annotations(master_dat={"McrA": [self.cs1, self.cs2, self.cs3]},
                                   marker_tree_info={"Function": {"McrA": {"SCO": {290, 289},
                                                                           "Aceticlastic": {340, 2}}}})
        self.assertEqual({"Function": "SCO"}, self.cs1.layers)
        self.assertEqual({"Function": "Unknown"}, self.cs2.layers)
        self.assertEqual({"Function": "Aceticlastic"}, self.cs3.layers)

        # Fails because the value in the marker_tree_info is a set, not a dictionary
        with pytest.raises(TypeError):
            map_queries_to_annotations(master_dat={"McrA": [self.cs1, self.cs2, self.cs3]},
                                       marker_tree_info={"Function": {"McrA": {290}}})

        return


if __name__ == '__main__':
    unittest.main()

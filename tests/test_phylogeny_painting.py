import os
import unittest
import pytest

from .testing_utils import get_test_data


class PhyPainterTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.refpkg import ReferencePackage
        self.mcra_refpkg = ReferencePackage()
        self.mcra_refpkg.f__pkl = get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.mcra_refpkg.slurp()
        return

    def test_filter_unwanted_taxa(self):
        from treesapp.phylogeny_painting import PhyPainter
        taxon_leaf_map, unique_taxa = self.mcra_refpkg.map_rank_representatives_to_leaves(rank_name='genus')
        filtered_taxa = PhyPainter.filter_unwanted_taxa(taxon_leaf_map, unique_taxa,
                                                        taxa_filter="Crenarchaeota")
        self.assertEqual(0, len(filtered_taxa))

        filtered_taxa = PhyPainter.filter_unwanted_taxa(taxon_leaf_map, unique_taxa,
                                                        taxa_filter="Methanomassiliicoccales")
        self.assertEqual(4, len(filtered_taxa))
        return

    def test_primer(self):
        from treesapp.treesapp_args import TreeSAPPArgumentParser, add_colour_arguments
        from treesapp.phylogeny_painting import PhyPainter
        # Set up the argument parser
        arg_parser = TreeSAPPArgumentParser()
        add_colour_arguments(arg_parser)

        painter = PhyPainter()
        painter.primer(args=arg_parser.parse_args(["-r", self.mcra_refpkg.f__pkl,
                                                   "-n", "Pathway",
                                                   "--unknown_colour", "white"]))
        self.assertEqual(1, len(painter.refpkg_dict))
        self.assertEqual("Pathway", painter.feature_name)
        self.assertEqual('#FFFFFF', painter.unknown_col)

        return

    def test_set_unknown_colour_from_mpl(self):
        from treesapp.phylogeny_painting import PhyPainter
        painter = PhyPainter()
        with pytest.raises(SystemExit):
            painter.set_unknown_colour_from_mpl("rainbow")

        painter.set_unknown_colour_from_mpl("teal")
        self.assertEqual('#008080', painter.unknown_col)
        return

    def test_find_mono_clades(self):
        from treesapp.phylogeny_painting import PhyPainter
        painter = PhyPainter()
        node_map = painter.find_mono_clades(taxon_leaf_map={'f__Methanomethyliaceae': ['67_McrA', '68_McrA', '69_McrA'],
                                                            'f__Methanocorpusculaceae': ['159_McrA', '160_McrA', '161_McrA'],
                                                            'f__Methanocellaceae': ['154_McrA', '155_McrA', '156_McrA']},
                                            ref_pkg=self.mcra_refpkg)
        leaves_included = 0
        for _taxon, mono_clade_list in node_map.items():
            for mono_clade in mono_clade_list:
                leaves_included += len(mono_clade)
        self.assertEqual(9, leaves_included)
        return


class ColourTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.refpkg import ReferencePackage
        self.mcra_metabolism = get_test_data("Mcr_taxonomy-phenotype_map.tsv")

        self.mcra_refpkg = ReferencePackage()
        self.mcra_refpkg.f__pkl = get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.mcra_refpkg.slurp()
        return

    def test_convert_taxa_to_phenotypes(self):
        from treesapp.file_parsers import read_phenotypes
        from treesapp.phylogeny_painting import convert_taxa_to_phenotypes

        phenotypes = read_phenotypes(self.mcra_metabolism)
        self.assertEqual(20, len(phenotypes))
        taxon_leaf_map = self.mcra_refpkg.map_taxa_to_leaf_nodes(list(phenotypes.keys()))
        phenotype_leaf_map = convert_taxa_to_phenotypes(phenotypes_map=phenotypes,
                                                        taxon_leaf_map=taxon_leaf_map)
        self.assertEqual(6, len(phenotype_leaf_map))
        self.assertEqual(19, len(phenotype_leaf_map["Aceticlastic"]))
        return


if __name__ == '__main__':
    unittest.main()

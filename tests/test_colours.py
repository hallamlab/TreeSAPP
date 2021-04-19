import os
import unittest


class ColourTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from .testing_utils import get_test_data
        from treesapp.refpkg import ReferencePackage
        self.mcra_pkl = get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.mcra_metabolism = get_test_data("Mcr_taxonomy-phenotype_map.tsv")

        self.mcra_refpkg = ReferencePackage()
        self.mcra_refpkg.f__pkl = self.mcra_pkl
        self.mcra_refpkg.slurp()
        return

    def test_convert_taxa_to_phenotypes(self):
        from treesapp.file_parsers import read_phenotypes
        from treesapp.phylogeny_painting import convert_taxa_to_phenotypes

        phenotypes = read_phenotypes(self.mcra_metabolism)
        taxon_leaf_map = self.mcra_refpkg.map_taxa_to_leaf_nodes(list(phenotypes.keys()))
        phenotype_leaf_map = convert_taxa_to_phenotypes(phenotypes_map=phenotypes,
                                                        taxon_leaf_map=taxon_leaf_map)
        self.assertEqual(5, len(phenotype_leaf_map))
        self.assertEqual(19, len(phenotype_leaf_map["Aceticlastic"]))
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


if __name__ == '__main__':
    unittest.main()

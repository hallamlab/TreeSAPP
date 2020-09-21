import os
import unittest


class ColourTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from .testing_utils import get_test_data
        from treesapp.refpkg import ReferencePackage
        self.mcra_pkl = get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        self.mcra_metabolism = get_test_data("Mcr_taxonomy-phenotype_map.tsv")

        self.mcra_refpkg = ReferencePackage()
        self.mcra_refpkg.f__json = self.mcra_pkl
        self.mcra_refpkg.slurp()
        return

    def test_map_taxa_to_leaf_nodes(self):
        from treesapp.phylogeny_painting import map_taxa_to_leaf_nodes

        leaf_map = map_taxa_to_leaf_nodes(leaf_names=["g__Methanothrix", "g__Methanosarcina"], refpkg=self.mcra_refpkg)
        self.assertEqual(4, len(leaf_map["g__Methanothrix"]))
        self.assertEqual(11, len(leaf_map["g__Methanosarcina"]))
        return

    def test_read_phenotypes_map(self):
        from treesapp.phylogeny_painting import read_phenotypes

        phenotypes = read_phenotypes(self.mcra_metabolism)
        self.assertIsInstance(cls=dict, obj=phenotypes)
        self.assertEqual(15, len(phenotypes))

    def test_convert_taxa_to_phenotypes(self):
        from treesapp.phylogeny_painting import convert_taxa_to_phenotypes, read_phenotypes, map_taxa_to_leaf_nodes

        phenotypes = read_phenotypes(self.mcra_metabolism)
        phenotype_leaf_map = convert_taxa_to_phenotypes(phenotypes_map=phenotypes,
                                                        taxon_leaf_map=map_taxa_to_leaf_nodes(list(phenotypes.keys()),
                                                                                              self.mcra_refpkg))
        self.assertEqual(5, len(phenotype_leaf_map))
        self.assertEqual(15, len(phenotype_leaf_map["Aceticlastic"]))
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

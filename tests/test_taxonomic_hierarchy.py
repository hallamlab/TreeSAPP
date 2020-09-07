import pytest
import unittest

from treesapp.phylo_seq import TreeLeafReference

_dup_lineage = "cellular organisms; Bacteria; Terrabacteria group; Actinobacteria; Actinobacteria;" \
               " Actinomycetales; Actinomycetaceae; Actinomyces; Actinomyces nasicola"
_dup_lineage_ex = [{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'},
                   {'TaxId': '2', 'ScientificName': 'Bacteria', 'Rank': 'superkingdom'},
                   {'TaxId': '1783272', 'ScientificName': 'Terrabacteria group', 'Rank': 'no rank'},
                   {'TaxId': '201174', 'ScientificName': 'Actinobacteria', 'Rank': 'phylum'},
                   {'TaxId': '1760', 'ScientificName': 'Actinobacteria', 'Rank': 'class'},
                   {'TaxId': '2037', 'ScientificName': 'Actinomycetales', 'Rank': 'order'},
                   {'TaxId': '2049', 'ScientificName': 'Actinomycetaceae', 'Rank': 'family'},
                   {'TaxId': '1654', 'ScientificName': 'Actinomyces', 'Rank': 'genus'},
                   {'ScientificName': "Actinomyces nasicola", 'Rank': "species"}]
_trunc_lineage = "cellular organisms; Bacteria; Terrabacteria group; Actinobacteria; Actinobacteria"
_trunc_lineage_ex = [{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'},
                     {'TaxId': '2', 'ScientificName': 'Bacteria', 'Rank': 'superkingdom'},
                     {'TaxId': '1783272', 'ScientificName': 'Terrabacteria group', 'Rank': 'no rank'},
                     {'TaxId': '201174', 'ScientificName': 'Actinobacteria', 'Rank': 'phylum'},
                     {'TaxId': '1760', 'ScientificName': 'Actinobacteria', 'Rank': 'class'}]
_diff_lineage = "cellular organisms; Bacteria; Terrabacteria group; Actinobacteria; Actinobacteria;" \
                " Bifidobacteriales; Bifidobacteriaceae; Bifidobacterium"
_diff_lineage_ex = [{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'},
                    {'TaxId': '2', 'ScientificName': 'Bacteria', 'Rank': 'superkingdom'},
                    {'TaxId': '1783272', 'ScientificName': 'Terrabacteria group', 'Rank': 'no rank'},
                    {'TaxId': '201174', 'ScientificName': 'Actinobacteria', 'Rank': 'phylum'},
                    {'TaxId': '1760', 'ScientificName': 'Actinobacteria', 'Rank': 'class'},
                    {'TaxId': '85004', 'ScientificName': 'Bifidobacteriales', 'Rank': 'order'},
                    {'TaxId': '31953', 'ScientificName': 'Bifidobacteriaceae', 'Rank': 'family'},
                    {'TaxId': '1678', 'ScientificName': 'Bifidobacterium', 'Rank': 'genus'}]
_conflict_node_one = TreeLeafReference('1', 'clashy')
_conflict_node_one.lineage = "d__Bacteria; p__Proteobacteria; n__environmental samples"
_conflict_node_two = TreeLeafReference('2', 'clashier')
_conflict_node_two.lineage = "d__Bacteria; p__Firmicutes; c__Fake; n__environmental samples"
_conflict_node_three = TreeLeafReference('3', 'clashiest')
_conflict_node_three.lineage = "d__Bacteria; p__Firmicutes; c__Fake; n__environmental samples"
_conflict_node_four = TreeLeafReference('4', 'redirect')
_conflict_node_four.lineage = "d__Bacteria; p__Firmicutes; n__RemoveMe; c__Fake"


@pytest.fixture(scope="class")
def th_class(request):
    from treesapp import taxonomic_hierarchy
    request.cls.db = taxonomic_hierarchy.TaxonomicHierarchy()
    request.cls.db.feed(lineage=_dup_lineage, lineage_ex=_dup_lineage_ex)
    request.cls.db.feed(lineage=_trunc_lineage, lineage_ex=_trunc_lineage_ex)
    request.cls.db.feed(lineage=_diff_lineage, lineage_ex=_diff_lineage_ex)
    return


@pytest.mark.usefixtures("th_class")
class TaxonomicHierarchyTester(unittest.TestCase):
    def test_fixture(self):
        self.assertEqual(True, hasattr(self, "db"))

    def test_feed(self):
        """
        Since lineages are already fed into self.db (via TaxonomicHierarchy.feed()), other functions are used for
        robustly testing whether feed() is operating as expected.
        This function currently just ensures lineages with incorrect LineageEx dictionaries will not be entered.
        :return: None
        """
        with pytest.raises(SystemExit):
            self.db.feed("Archaea; Crenarchaeota", [{'ScientificName': 'Archaea', 'Rank': 'superkingdom'}])
        with pytest.raises(SystemExit):
            self.db.feed("Archaea; Crenarchaeota", [{'ScientificName': 'Archaea', 'Rank': 'superkingdom'},
                                                    {'ScientificName': 'Euryarchaeota', 'Rank': 'phylum'}])
        return

    def test_feed_leaf_nodes(self):
        from treesapp.taxonomic_hierarchy import TaxonomicHierarchy
        test_th = TaxonomicHierarchy()
        test_th.feed_leaf_nodes([_conflict_node_one, _conflict_node_two, _conflict_node_three])
        post_feed = test_th.lineages_fed
        self.assertEqual(3, post_feed)
        self.assertTrue("n__environmental samples_1" in test_th.hierarchy)
        self.assertTrue("n__environmental samples_2" not in test_th.hierarchy)
        self.assertEqual(0, len(test_th.conflicts))
        return

    def test_emit(self):
        actino_line = "Bacteria; Terrabacteria group; Actinobacteria;" \
                      " Actinobacteria; Actinomycetales; Actinomycetaceae; Actinomyces"
        bifido_line = "Bacteria; Terrabacteria group; Actinobacteria;" \
                      " Actinobacteria; Bifidobacteriales; Bifidobacteriaceae; Bifidobacterium"
        self.assertEqual(actino_line, self.db.emit("g__Actinomyces"))
        self.assertEqual(bifido_line, self.db.emit("g__Bifidobacterium"))
        return

    def test_check_lineage(self):
        t1_lineage = "n__cellular organisms; d__Bacteria; n__Terrabacteria group; p__Actinobacteria; c__Actinobacteria"
        t1_organism = "Actinobacteria"
        t2_lineage = "d__Bacteria; p__Actinobacteria; c__Actinobacteria; " \
                     "o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces"
        t2_organism = "s__Actinomyces nasicola"
        self.assertEqual(3,
                         len(self.db.check_lineage(lineage=t1_lineage,
                                                   organism=t1_organism).split(self.db.lin_sep)))
        self.assertEqual(7,
                         len(self.db.check_lineage(lineage=t2_lineage,
                                                   organism=t2_organism).split(self.db.lin_sep)))
        return

    def test_check_lineage_nonexistant(self):
        with pytest.raises(RuntimeError):
            self.db.check_lineage(lineage="d__Archaea", organism="Archaea")
        return

    def test_check_lineage_rankless(self):
        lin = "d__Bacteria; p__Cyanobacteria; o__Synechococcales; f__Prochloraceae; g__Prochlorococcus"
        org = "s__Prochlorococcus marinus"
        self.db.feed("Bacteria; Cyanobacteria; Synechococcales; Prochloraceae; Prochlorococcus",
                     [{'ScientificName': 'Bacteria', 'Rank': 'domain'},
                      {'ScientificName': 'Cyanobacteria', 'Rank': 'phylum'},
                      {'ScientificName': 'Synechococcales', 'Rank': 'order'},
                      {'ScientificName': 'Prochloraceae', 'Rank': 'family'},
                      {'ScientificName': 'Prochlorococcus', 'Rank': 'genus'}])
        self.assertEqual("d__Bacteria; p__Cyanobacteria", self.db.check_lineage(lineage=lin, organism=org))
        return

    def test_clean_lineage_string(self):
        l1 = "n__cellular organisms; d__Bacteria; n__Terrabacteria group"
        l2 = "d__Bacteria"
        l3 = "d__Bacteria; p__Candidatus Omnitrophica; s__Candidatus Omnitrophica bacterium CG02__42_8"
        self.assertEqual(self.db.clean_lineage_string(l1, with_prefix=True),
                         self.db.clean_lineage_string(l2, with_prefix=True))
        self.assertEqual("Bacteria", self.db.clean_lineage_string(l1, with_prefix=False))
        self.assertEqual("d__Bacteria; p__Candidatus Omnitrophica; s__Candidatus Omnitrophica bacterium CG02_42_8",
                         self.db.clean_lineage_string(l3, with_prefix=True))
        return

    def test_remove_leaf_nodes(self):
        self.db.feed("Archaea; Euryarchaeota", [{'ScientificName': 'Archaea', 'Rank': 'superkingdom'},
                                                {'ScientificName': 'Euryarchaeota', 'Rank': 'phylum'}])
        before = self.db.lineages_fed
        self.db.remove_leaf_nodes("p__Euryarchaeota")  # To test whether conversion of single string to set works
        self.assertEqual(before-1, self.db.lineages_fed)
        self.assertRaises(KeyError, lambda: self.db.hierarchy["d__Archaea"])

    def test_trim_lineages_to_rank(self):
        dict_in = {'1': 'd__Archaea; p__Candidatus Bathyarchaeota',
                   '2': 'd__Archaea; p__Candidatus Micrarchaeota',
                   '3': 'd__Archaea; p__Euryarchaeota',
                   '4': 'd__Bacteria',
                   '7': 'd__Bacteria; p__Acidobacteria',
                   '9': 'd__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales;'
                        ' f__Actinomycetaceae; g__Actinomyces; s__Actinomyces nasicola'}
        out = self.db.trim_lineages_to_rank(dict_in, "class")
        self.assertEqual({'9': 'd__Bacteria; p__Actinobacteria; c__Actinobacteria'}, out)

    def test_rank_representatives(self):
        self.assertEqual({'Bacteria'}, self.db.rank_representatives("domain"))
        self.assertEqual({'Actinobacteria'}, self.db.rank_representatives("class"))
        self.assertEqual({'Actinomyces nasicola'}, self.db.rank_representatives("species"))

    def test_unprefix_rank_existance(self):
        self.db.clean_trie = True
        self.db.build_multifurcating_trie(key_prefix=False, value_prefix=True)
        self.assertTrue("Bacteria" in self.db.trie)
        self.assertEqual("c__Actinobacteria", self.db.trie["Bacteria; Actinobacteria; Actinobacteria"])

    def test_rm_bad_taxa(self):
        self.assertEqual(["Archaea"], self.db.rm_bad_taxa_from_lineage(["cellular organisms", "Archaea"]))

    def test_get_prefixed_lineage_from_bare(self):
        self.db.clean_trie = True
        self.db.build_multifurcating_trie(key_prefix=False, value_prefix=True)
        self.assertEqual("d__Bacteria",
                         self.db.get_prefixed_lineage_from_bare("cellular organisms; Bacteria; "
                                                                "Terrabacteria group; Actinobacteria"))
        self.assertEqual("d__Bacteria; p__Actinobacteria; c__Actinobacteria",
                         self.db.get_prefixed_lineage_from_bare("Bacteria; Actinobacteria; Actinobacteria"))

    def test_build_multifurcating_trie(self):
        self.db.trie_key_prefix = True
        self.db.build_multifurcating_trie(key_prefix=False, value_prefix=True)
        self.assertFalse(self.db.trie_key_prefix)
        self.assertTrue(self.db.trie_value_prefix)

    def test_strip_taxon(self):
        self.assertEqual("Bacteria; Proteobacteria", self.db.strip_rank_prefix("Bacteria; Proteobacteria"))
        self.assertEqual("Bacteria; Mock", self.db.strip_rank_prefix("d__Bacteria; n__Mock"))

    def test_resolve_conflicts(self):
        from treesapp.taxonomic_hierarchy import TaxonomicHierarchy
        test_th = TaxonomicHierarchy()
        test_leaf_nodes = [_conflict_node_two, _conflict_node_three, _conflict_node_four]
        test_th.feed_leaf_nodes(test_leaf_nodes)
        self.assertEqual(len(test_leaf_nodes), test_th.lineages_fed)
        test_th.resolve_conflicts()
        self.assertEqual(4, len(test_th.hierarchy))
        coverage_values = [test_th.hierarchy[t].coverage for t in test_th.hierarchy]
        self.assertTrue(len(test_leaf_nodes) == max(coverage_values))
        return

    # def test_rm_absent_taxa_from_lineage(self):
    #     return
    #
    # def test_evaluate_hierarchy_clash(self):
    #     from taxonomic_hierarchy import TaxonomicHierarchy
    #     return


class TaxonTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.taxonomic_hierarchy import Taxon
        self.lin_sep = "; "
        self.taxon_sep = "__"

        # Lineages used in tests:
        # "d__Archaea; p__Crenarchaeota; c__Thermoprotei; n__environmental samples"
        # "d__Archaea; p__Euryarchaeota; c__Methanomicrobia; o__Methanomicrobiales; n__environmental samples"

        # Instantiate some Taxon class objects
        self.t_arc = Taxon("Archaea", "domain")
        self.t_cren = Taxon("Crenarchaeota", "phylum")
        self.t_therm = Taxon("Thermoprotei", "class")

        self.t_eury = Taxon("Euryarchaeota", "phylum")
        self.t_methia = Taxon("Methanomicrobia", "class")
        self.t_methales = Taxon("Methanomicrobiales", "order")
        self.t_envsam = Taxon("environmental samples", "no rank")

        self.t_bac = Taxon("Bacteria", "domain")
        self.t_bac.parent = None

        # Populate their prefix and parent attributes
        previous = None
        for taxon in [self.t_arc, self.t_cren, self.t_therm]:  # type: Taxon
            taxon.prefix = taxon.rank[0]
            taxon.parent = previous
            previous = taxon
        previous = None
        for taxon in [self.t_arc, self.t_eury, self.t_methia, self.t_methales, self.t_envsam]:  # type: Taxon
            taxon.prefix = taxon.rank[0]
            taxon.parent = previous
            previous = taxon

        # Off to the races
        return

    def test_lca(self):
        from treesapp.taxonomic_hierarchy import Taxon
        self.assertEqual("Archaea", Taxon.lca(self.t_eury, self.t_cren).name)
        self.assertEqual("Archaea", Taxon.lca(self.t_envsam, self.t_therm).name)
        self.assertEqual("Methanomicrobia", Taxon.lca(self.t_methia, self.t_methales).name)
        self.assertEqual(None, Taxon.lca(self.t_bac, self.t_arc))
        return

    def test_tax_dist(self):
        self.assertEqual(0, self.t_arc.tax_dist(self.t_arc))
        self.assertEqual(0, self.t_cren.tax_dist(self.t_cren))
        self.assertEqual(1, self.t_bac.tax_dist(self.t_arc))
        self.assertEqual(2, self.t_arc.tax_dist(self.t_methia))
        self.assertEqual(2, self.t_methia.tax_dist(self.t_arc))
        return

    def test_lineage_slice(self):
        from treesapp.taxonomic_hierarchy import Taxon
        self.assertEqual(2, len(Taxon.lineage_slice(self.t_therm, self.t_arc)))
        return


if __name__ == '__main__':
    unittest.main()

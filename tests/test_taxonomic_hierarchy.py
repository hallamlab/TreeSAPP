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

    def test_check_lineage(self):
        t1_lineage = "n__cellular organisms; d__Bacteria; n__Terrabacteria group; p__Actinobacteria; c__Actinobacteria"
        t1_organism = "Actinobacteria"
        t2_lineage = "d__Bacteria; p__Actinobacteria; c__Actinobacteria; " \
                     "o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces"
        t2_organism = "s__Actinomyces nasicola"
        t3_lineage = "r__Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales"
        t4_lineage = "d__Bacteria; p__Cyanobacteria; o__Synechococcales; f__Prochloraceae; g__Prochlorococcus"
        t4_organism = "s__Prochlorococcus marinus"

        # Remove the taxa that have no taxonomic rank
        self.assertEqual(4,
                         len(self.db.check_lineage(lineage=t1_lineage, organism=t1_organism,
                                                   verbosity=1).split(self.db.lin_sep)))
        # Add the organism as the new species
        self.assertEqual(8,
                         len(self.db.check_lineage(lineage=t2_lineage, organism=t2_organism,
                                                   verbosity=1).split(self.db.lin_sep)))
        # Nothing to change
        self.assertEqual(t3_lineage,
                         self.db.check_lineage(lineage=t3_lineage, organism=t3_lineage.split(self.db.lin_sep)[-1],
                                               verbosity=1))
        # Test when a taxon is out of order
        self.db.feed("Bacteria; Cyanobacteria; Synechococcales; Prochloraceae; Prochlorococcus",
                     [{'ScientificName': 'Bacteria', 'Rank': 'domain'},
                      {'ScientificName': 'Cyanobacteria', 'Rank': 'phylum'},
                      {'ScientificName': 'Synechococcales', 'Rank': 'order'},
                      {'ScientificName': 'Prochloraceae', 'Rank': 'family'},
                      {'ScientificName': 'Prochlorococcus', 'Rank': 'genus'}])
        self.assertEqual("r__Root; d__Bacteria; p__Cyanobacteria",
                         self.db.check_lineage(lineage=t4_lineage, organism=t4_organism))

        # Test a lineage that doesn't exist in the hierarchy
        with pytest.raises(RuntimeError):
            self.db.check_lineage(lineage="d__Archaea", organism="Archaea")
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

    def test_emit(self):
        actino_line = "Root; Bacteria; Terrabacteria group; Actinobacteria;" \
                      " Actinobacteria; Actinomycetales; Actinomycetaceae; Actinomyces"
        bifido_line = "r__Root; d__Bacteria; n__Terrabacteria group; p__Actinobacteria;" \
                      " c__Actinobacteria; o__Bifidobacteriales; f__Bifidobacteriaceae; g__Bifidobacterium"
        self.assertEqual(actino_line, self.db.emit("g__Actinomyces"))
        self.assertEqual(bifido_line, self.db.emit("g__Bifidobacterium", with_prefix=True))
        return

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

        # Test loading a TreeLeafReference instance with unclassified taxa
        unclassified_leaf = TreeLeafReference('5', '2667527398_667992784')
        unclassified_leaf.lineage = "d__Archaea; p__Euryarchaeota; c__Methanomicrobia; o__Methanomicrobiales; f__unclassified; g__unclassified"
        test_th.feed_leaf_nodes([unclassified_leaf])
        self.assertEqual(4, test_th.lineages_fed)
        self.assertEqual(1, test_th.get_taxon("d__Archaea").coverage)
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

        dict_in.update({"10": 'd__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales;'
                        ' f__Actinomycetaceae; s__Actinomyces nasicola'})
        self.assertEqual(1, len(self.db.trim_lineages_to_rank(dict_in, "genus")))
        return

    def test_rank_representatives(self):
        self.assertEqual({'Bacteria'}, self.db.rank_representatives("domain"))
        self.assertEqual({'Actinobacteria'}, self.db.rank_representatives("class"))
        self.assertEqual({'Actinomyces nasicola'}, self.db.rank_representatives("species"))

    def test_unprefix_rank_existance(self):
        self.db.clean_trie = True
        self.db.build_multifurcating_trie(key_prefix=False, value_prefix=True)
        self.assertTrue("Bacteria" in self.db.trie)
        self.assertEqual("c__Actinobacteria", self.db.trie["Bacteria; Actinobacteria; Actinobacteria"])
        return

    def test_rm_bad_taxa_from_lineage(self):
        self.assertEqual(["Archaea"], self.db.rm_bad_taxa_from_lineage(["cellular organisms", "Archaea"]))
        return

    def test_rm_absent_taxa_from_lineage(self):
        ll_1 = ["d__Bacteria", "n__Absent", "p__Actinobacteria"]
        ll_2 = ["Bacteria", "Actinomyces", "Actinobacteria"]
        ll_3 = "d__Bacteria; p__Actinobacteria; c__Actinobacteria".split(self.db.lin_sep)
        # Test a lineage containing a taxon that is not present in the taxonomic hierarchy
        self.db.rm_absent_taxa_from_lineage(ll_1, True)
        self.assertEqual(2, len(ll_1))
        # Test a lineage where the taxa are out of order
        self.db.rm_absent_taxa_from_lineage(ll_2)
        self.assertEqual(3, len(ll_2))
        # Test a lineage where all taxa are in the hierarchy -> lineage should be unchanged
        self.db.rm_absent_taxa_from_lineage(ll_3, True)
        self.assertEqual(3, len(ll_3))
        return

    def test_get_prefixed_lineage_from_bare(self):
        self.db.clean_trie = True
        self.db.build_multifurcating_trie(key_prefix=False, value_prefix=True)
        # Determine whether a prefix for the root taxon is needed for the expected values
        if self.db.rooted:
            root_name = self.db.root_taxon + self.db.lin_sep
        else:
            root_name = ""
        # Test a lineage with taxa with no rank
        self.assertEqual(root_name + "d__Bacteria",
                         self.db.get_prefixed_lineage_from_bare("cellular organisms; Bacteria; "
                                                                "Terrabacteria group; Actinobacteria"))
        # Test a lineage where two ranks share the same name
        self.assertEqual(root_name + "d__Bacteria; p__Actinobacteria; c__Actinobacteria",
                         self.db.get_prefixed_lineage_from_bare("Bacteria; Actinobacteria; Actinobacteria"))
        # Test a lineage that isn't found in the taxonomic hierarchy
        with pytest.raises(SystemExit):
            self.db.get_prefixed_lineage_from_bare("Archaea; Lokiarchaeota")
        return

    def test_get_bare_taxon(self):
        with pytest.raises(RuntimeError):
            self.db.get_bare_taxon("d__Bacteria")
        with pytest.raises(SystemExit):
            self.db.get_bare_taxon("E. coli")
        taxon = self.db.get_bare_taxon("Actinobacteria")
        self.assertIsNone(taxon)
        taxon = self.db.get_bare_taxon("Bacteria")
        self.assertEqual("d__Bacteria", taxon.prefix_taxon())
        return

    def test_reroot_lineage(self):
        # Test without rank-prefixes on the taxa - should exit
        with pytest.raises(RuntimeError):
            self.db.reroot_lineage("Bacteria; Actinobacteria; Actinobacteria")
        # Test with normal behaviour - adds the root taxon
        self.assertEqual("r__Root; d__Bacteria", self.db.reroot_lineage("d__Bacteria"))
        # Test when the lineage is already rooted - nothing should be changed
        self.assertEqual("r__Root; d__Bacteria", self.db.reroot_lineage("r__Root; d__Bacteria"))
        return

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

    def test_root_domains(self):
        from treesapp.taxonomic_hierarchy import Taxon
        with pytest.raises(TypeError):
            self.db.root_domains("Root")
        root_taxon = Taxon(name="Root", rank="root")
        root_taxon = self.db.root_domains(root_taxon)
        self.assertEqual(root_taxon, self.db.hierarchy["d__Bacteria"].parent)
        # Remove the Taxon 'r__Root' so other tests are not affected
        self.db.redirect_hierarchy_paths(root_taxon)
        return

    def test_match_organism(self):
        # Prepare the lineages for testing, depending on whether the domains are rooted or not
        ll_1 = ["d__Bacteria"]
        ll_2 = ["d__Eukaryota"]
        ll_3 = "d__Bacteria; p__Actinobacteria; c__Actinobacteria".split('; ')
        if self.db.rooted:
            ll_1 = [self.db.root_taxon] + ll_1
            ll_2 = [self.db.root_taxon] + ll_2
            ll_3 = [self.db.root_taxon] + ll_3

        # Test when organism isn't in the taxonomic hierarchy
        organism = self.db.match_organism(organism="doesn't exist",
                                          lineage=ll_1)
        self.assertEqual("", organism)
        # Test when the lineage isn't in the taxonomic hierarchy
        with pytest.raises(SystemExit):
            self.db.match_organism(organism="Homo sapiens", lineage=ll_2)

        # Test when the organism isn't a descendent of the lineage
        organism = self.db.match_organism(organism="Methanosarcina barkeri",
                                          lineage=ll_3)
        self.assertEqual("", organism)
        # Test when the lineage is the wrong type (string)
        with pytest.raises(TypeError):
            self.db.match_organism(organism="", lineage="")

        # Behaviour when inputs are correct
        organism = self.db.match_organism(organism="Bifidobacterium",
                                          lineage=ll_3)
        self.assertEqual("g__Bifidobacterium", organism)
        return

    def test_redirect_hierarchy_paths(self):
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

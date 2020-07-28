import pytest
import unittest

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

    def test_emit(self):
        actino_line = "cellular organisms; Bacteria; Terrabacteria group; Actinobacteria;" \
                      " Actinobacteria; Actinomycetales; Actinomycetaceae; Actinomyces"
        bifido_line = "cellular organisms; Bacteria; Terrabacteria group; Actinobacteria;" \
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
        with pytest.raises(SystemExit):
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
        self.assertEqual(["Archaea"], self.db.rm_bad_taxa(["cellular organisms", "Archaea"]))

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


if __name__ == '__main__':
    unittest.main()

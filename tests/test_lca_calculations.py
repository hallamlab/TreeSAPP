import unittest


class LCATester(unittest.TestCase):
    def setUp(self) -> None:
        from pygtrie import StringTrie
        self.lineage_trie = StringTrie(separator="; ")
        lineages = ["r__Root; d__Bacteria; p__Proteobacteria; c__Deltaproteobacteria",
                    "r__Root; d__Bacteria; p__Verrucomicrobia; c__Verrucomicrobiae; o__Verrucomicrobiales",
                    "r__Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Bifidobacteriales"]
        for lin in lineages:
            taxa = lin.split("; ")
            while taxa:
                self.lineage_trie["; ".join(taxa)] = True
                taxa.pop(-1)

        return

    def test_identify_excluded_clade(self):
        from treesapp import lca_calculations
        # Test query lineages without 'r__Root'
        red = lca_calculations.identify_excluded_clade(assignment_dict={"r__Root; d__Bacteria; p__Proteobacteria":
                                                                            ["d__Bacteria; p__Verrucomicrobia; c__Verrucomicrobiae"]},
                                                       trie=self.lineage_trie)
        self.assertEqual(1, len(red["order"]))

        red = lca_calculations.identify_excluded_clade(assignment_dict={"r__Root; d__Bacteria; p__Proteobacteria":
                                                                            ["r__Root; d__Bacteria; p__Verrucomicrobia; c__Verrucomicrobiae; o__Verrucomicrobiales"]},
                                                       trie=self.lineage_trie)
        self.assertEqual(1, len(red["family"]))
        self.assertEqual(9, len(red))

        return

    def test_taxonomic_distinctness(self):
        from treesapp import lca_calculations
        from treesapp import taxonomic_hierarchy
        hierarchy = taxonomic_hierarchy.TaxonomicHierarchy()
        hierarchy.validate_rank_prefixes()
        hierarchy.rank_prefix_map.update({'d': "domain", 'p': "phylum", 'c': "class", 'o': "order"})
        lineages = ["r__Root; d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria",
                    "r__Root; d__Archaea; p__Helarchaeota",
                    "r__Root; d__Bacteria; p__Verrucomicrobia; n__fake; c__Verrucomicrobiae; o__Verrucomicrobiales"]
        for lin in lineages:
            prev = None
            for t in lin.split(hierarchy.lin_sep):
                prefix, taxon = t.split(hierarchy.taxon_sep)
                tt = taxonomic_hierarchy.Taxon(name=taxon, rank=hierarchy.rank_prefix_map[prefix])
                tt.parent = prev
                hierarchy.hierarchy[t] = tt
                prev = tt
        queries = {"q1": hierarchy.get_taxon("o__Verrucomicrobiales"),
                   "q2": hierarchy.get_taxon("o__Verrucomicrobiales"),
                   "q3": hierarchy.get_taxon("p__Helarchaeota"),
                   "q4": hierarchy.get_taxon("c__Alphaproteobacteria")}
        td = lca_calculations.taxonomic_distinctness(query_taxa=queries,
                                                     rank="class",
                                                     rank_depths=hierarchy.accepted_ranks_depths)
        self.assertEqual(1.25, td)
        return


if __name__ == '__main__':
    unittest.main()

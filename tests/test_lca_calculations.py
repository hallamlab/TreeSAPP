import unittest


class MyTestCase(unittest.TestCase):
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
        from treesapp.lca_calculations import identify_excluded_clade
        # Test query lineages without 'r__Root'
        red = identify_excluded_clade(assignment_dict={"r__Root; d__Bacteria; p__Proteobacteria":
                                                           ["d__Bacteria; p__Verrucomicrobia; c__Verrucomicrobiae"]},
                                      trie=self.lineage_trie)
        self.assertEqual(1, len(red["order"]))

        red = identify_excluded_clade(assignment_dict={"r__Root; d__Bacteria; p__Proteobacteria":
                                                           ["r__Root; d__Bacteria; p__Verrucomicrobia; c__Verrucomicrobiae; o__Verrucomicrobiales"]},
                                      trie=self.lineage_trie)
        self.assertEqual(1, len(red["family"]))
        self.assertEqual(9, len(red))

        return


if __name__ == '__main__':
    unittest.main()

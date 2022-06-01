import os
import pytest
import shutil
import unittest

from Bio import Entrez

from . import testing_utils as utils


class EntrezUtilitiesTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.fasta import FASTA
        from treesapp.classy import Creator
        from treesapp.entrez_utils import EntrezRecord

        self.test_fa = FASTA(utils.get_test_data("create_test.faa"))
        self.test_fa.load_fasta()
        # Subset to complete quicker
        self.test_fa.keep_only(self.test_fa.get_seq_names()[0:18])

        self.create_inst = Creator()
        self.create_inst.output_dir = "./tests/entrez_utils_test/"
        self.create_inst.final_output_dir = self.create_inst.output_dir + "final_outputs" + os.sep
        self.create_inst.var_output_dir = self.create_inst.output_dir + "intermediates" + os.sep
        self.create_inst.check_previous_output()
        self.create_inst.current_stage = self.create_inst.stages[1]
        self.create_inst.acc_to_lin = os.path.join(self.create_inst.output_dir, "test_create_acc_to_lin.tsv")

        self.accession2taxid = utils.get_test_data("create_test.accession2taxid")
        self.test_entrez_records = []
        test_accs = ["BAJ94456", "CUW39146.1", "NP_001076868", "WP_056230317.1", "XP_005707548.1", "WP_000101794"]
        for acc in test_accs:
            self.test_entrez_records.append(EntrezRecord(acc=acc, ver=""))
        return

    def tearDown(self) -> None:
        if os.path.isdir(self.create_inst.output_dir):
            shutil.rmtree(self.create_inst.output_dir)
        return

    def test_validate_target_db(self):
        from treesapp.entrez_utils import validate_target_db
        self.assertEqual("Taxonomy", validate_target_db("tax"))
        self.assertEqual("nucleotide", validate_target_db("dna"))
        with pytest.raises(SystemExit):
            validate_target_db("nuc")
        return

    def test_fetch_entrez_lineages(self):
        entrez_record_dict = self.create_inst.fetch_entrez_lineages(self.test_fa, 'prot')
        self.assertEqual(self.test_fa.n_seqs(), len(entrez_record_dict))
        passed = 0
        failed = 0
        for er in entrez_record_dict.values():
            if er.taxon_rank != "root":
                self.assertEqual(7, er.bitflag)
                passed += 1
            else:
                self.assertEqual(5, er.bitflag)
                failed += 1
        self.assertEqual(17, passed)
        self.assertEqual(1, failed)
        return

    def test_prep_for_entrez_query(self):
        from treesapp.entrez_utils import prep_for_entrez_query
        handle = prep_for_entrez_query()
        for record in Entrez.read(handle):
            self.assertTrue(record["TaxId"] == "158330")
            self.assertEqual("Cypripedioideae", record["ScientificName"])
        return

    def test_map_accession2taxid(self):
        from treesapp import entrez_utils as e_utils
        # Test failure if accession2taxid file doesn't exist
        with pytest.raises(SystemExit):
            e_utils.map_accession2taxid(self.test_entrez_records, "test_data/fake.accession2taxid")

        # Test normal operating conditions
        er_acc_dict = e_utils.map_accession2taxid(query_accessions=self.test_entrez_records,
                                                  accession2taxid_list=self.accession2taxid)
        self.assertEqual(6, len(er_acc_dict))
        self.assertEqual('112509', er_acc_dict["BAJ94456"].pop().ncbi_tax)
        self.assertEqual('', er_acc_dict["XP_005707548.1"].pop().ncbi_tax)

        # Clear for other tests
        for q in self.test_entrez_records:  # type: e_utils.EntrezRecord
            q.lineage = ""
        return

    def test_map_accessions_to_lineages(self):
        from treesapp import entrez_utils as e_utils
        from treesapp.taxonomic_hierarchy import TaxonomicHierarchy
        taxa_hrcy = TaxonomicHierarchy()

        # Test normal operating conditions
        e_utils.map_accessions_to_lineages(query_accession_list=self.test_entrez_records, t_hierarchy=taxa_hrcy,
                                           accession_to_taxid=self.accession2taxid, molecule="prot")
        self.assertEqual(75, len(taxa_hrcy.hierarchy))
        self.assertEqual(6, taxa_hrcy.lineages_fed)
        self.assertTrue("d__Bacteria" in taxa_hrcy.hierarchy)

        # Clear for other tests
        for q in self.test_entrez_records:  # type: e_utils.EntrezRecord
            q.lineage = ""
        return

    def test_map_orf_lineages(self):
        from treesapp.entrez_utils import map_orf_lineages
        from treesapp.fasta import Header
        h1 = Header('Q3J170')
        h1.first_split = h1.original
        h2 = Header('P06008_orf_2')
        h2.first_split = h2.original
        h3 = Header('P060081_orf_3')
        h3.first_split = h3.original
        h4 = Header('AB-746_I02_AB-902_NODE_1_length_133761_cov_569.254_ID_1_26 # 25226 # 25834 # 1 # ID=10766_26;partial=00;start_type=ATG;rbs_motif=AGGAG(G)/GGAGG;rbs_spacer=13-15bp;gc_cont=0.335')
        h4.first_split = h4.original.split()[0]
        h4.find_accession()
        h5 = Header("1147.D082_03820")
        h5.first_split = h5.original
        h5.find_accession()
        h6 = Header("SI072_135m_bin.14_k147_20226")
        h6.first_split = h6.original
        h6.find_accession()
        self.assertEqual("SI072_135m_bin.14_k147_20226", h6.version)
        self.assertEqual("AB-746_I02_AB-902_NODE_1_length_133761_cov_569.254_ID_1_26", h4.version)

        with pytest.raises(SystemExit):
            map_orf_lineages(seq_lineage_tbl=utils.get_test_data("SwissProt_PuhA_seqs2lineage.txt"), header_registry={})

        # Test when there are names that are prefixes of others and ensure they are mapped correctly
        orf_lin_map, found = map_orf_lineages(seq_lineage_tbl=utils.get_test_data("SwissProt_PuhA_seqs2lineage.txt"),
                                              header_registry={'1': h1, '2': h2, '3': h3})
        self.assertEqual(2, len(found))
        self.assertTrue(h2.original in orf_lin_map)
        self.assertEqual("d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales;"
                         " f__Rhodobacteraceae; g__Rhodobacter; s__Rhodobacter sphaeroides",
                         orf_lin_map[h1.first_split])

        # Test a complicated header
        orf_lin_map, found = map_orf_lineages(seq_lineage_tbl=utils.get_test_data("test_seqs2lineage.txt"),
                                              header_registry={'1': h4, '2': h5, '3': h6})
        self.assertEqual(3, len(found))
        return

    def test_ref_seq_lineage_scanner(self):
        from treesapp import entrez_utils
        er_one = entrez_utils.EntrezRecord(".", ".")
        er_one.lineage = 'r__Root; d__Bacteria; p__Tenericutes; c__Mollicutes; o__Mycoplasmatales;' \
                         ' f__Mycoplasmataceae; g__Mycoplasmopsis; s__Mycoplasma cynos; t__Mycoplasma cynos C142'
        matches = entrez_utils.ref_seq_lineage_scanner(ref_seq_dict={'1': er_one}, taxon_name="g__Mycoplasmopsis")
        self.assertEqual([er_one], matches)
        return

    def test_repair_lineages(self):
        from treesapp import entrez_utils
        from treesapp.refpkg import ReferencePackage
        # Load reference package to compare with
        mcra_refpkg = ReferencePackage()
        mcra_refpkg.f__pkl = utils.get_test_data(os.path.join("refpkgs", "McrA_build.pkl"))
        mcra_refpkg.slurp()
        # Set up mock test data
        l1 = "Archaea; Euryarchaeota; environmental samples"
        l2 = "Archaea; Euryarchaeota; Stenosarchaea group; Methanomicrobia; Methanosarcinales; " \
             "Candidatus Methanoperedenaceae; Candidatus Methanoperedens; environmental samples"
        l3 = "r__Root; d__Archaea; p__Candidatus Thermoplasmatota; c__Thermoplasmata_1; " \
             "o__Methanomassiliicoccales_1; f__Candidatus Methanomethylophilaceae; " \
             "g__Candidatus Methanomethylophilus; s__Candidatus Methanomethylophilus sp. 1R26"
        mock_one_er = entrez_utils.EntrezRecord("AAW80308", "AAW80308.1")
        mock_two_er = entrez_utils.EntrezRecord("QGT32507", "QGT32507.1")
        mock_three_er = entrez_utils.EntrezRecord("seq_3", "seq_3.1")
        mock_one_er.lineage = l1
        mock_two_er.lineage = l2
        mock_three_er.lineage = l3 + "; t__Methanomethylophilus sp. 1R26 strain mock"

        # Test
        entrez_utils.repair_lineages(ref_seq_dict={'1': mock_one_er, '2': mock_two_er, '3': mock_three_er},
                                     t_hierarchy=mcra_refpkg.taxa_trie)
        self.assertEqual("r__Root; d__Archaea; p__Euryarchaeota", mock_one_er.lineage)
        self.assertEqual("r__Root; d__Archaea; p__Euryarchaeota", mock_two_er.lineage)
        self.assertEqual(l3, mock_three_er.lineage)
        return

    def test_repair_conflict_lineages(self):
        from treesapp.entrez_utils import repair_conflict_lineages, EntrezRecord
        from treesapp.taxonomic_hierarchy import TaxonomicHierarchy
        from treesapp.phylo_seq import TreeLeafReference
        # Test an empty TaxonomicHierarchy with zero conflicts
        repair_conflict_lineages(t_hierarchy=self.create_inst.ref_pkg.taxa_trie, ref_seq_dict={})

        triple_conflict = ["p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium; s__Clostridium autoethanogenum",
                           "p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Alkaliphilus; s__Alkaliphilus metalliredigens",
                           "p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Selenomonas; s__Selenomonas flueggei",
                           'p__Firmicutes; c__Negativicutes; o__Selenomonadales; f__Veillonellaceae; g__Mitsuokella; s__Mitsuokella multacida',
                           'p__Firmicutes; c__Negativicutes; o__Selenomonadales; f__Selenomonadaceae; g__Mitsuokella; s__Mitsuokella multacida',
                           'p__Firmicutes; c__Negativicutes; o__Veillonellales; f__Veillonellaceae; g__Megasphaera; s__Megasphaera elsdenii']
        # Test a TaxonomicHierarchy with several conflicting taxa
        er_dict = {}
        leaf_nodes = []
        index = 1
        for lineage in triple_conflict:
            leaf = TreeLeafReference(str(index), 'test')
            er = EntrezRecord(acc=str(index), ver=str(index) + ".1")
            leaf.lineage = lineage
            er.lineage = lineage
            er_dict[str(index)] = er
            leaf_nodes.append(leaf)
            index += 1

        taxonomy = TaxonomicHierarchy()
        taxonomy.feed_leaf_nodes(leaf_nodes)

        # Retry with domain included
        for leaf in leaf_nodes:
            leaf.lineage = "d__Bacteria; " + leaf.lineage
        taxonomy.feed_leaf_nodes(leaf_nodes)

        # Orders: 'Veillonellales' -> 'Selenomonadales', 'Selenomonadales' -> 'Clostridiales'
        # Families: 'Selenomonadaceae' -> 'Veillonellaceae'
        # Change the coverage of the different lineages to control the representative
        for t in taxonomy.get_taxon("s__Clostridium autoethanogenum").lineage():
            t.coverage += 3
        taxonomy.conflicts = {(taxonomy.get_taxon("o__Veillonellales"), taxonomy.get_taxon("o__Selenomonadales")),
                              (taxonomy.get_taxon("o__Selenomonadales"), taxonomy.get_taxon("o__Clostridiales")),
                              (taxonomy.get_taxon("f__Veillonellaceae"), taxonomy.get_taxon("f__Selenomonadaceae"))}

        repair_conflict_lineages(taxonomy, er_dict)
        self.assertEqual("r__Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales;"
                         " f__Veillonellaceae; g__Megasphaera; s__Megasphaera elsdenii",
                         "; ".join([t.prefix_taxon() for t in taxonomy.get_taxon("s__Megasphaera elsdenii").lineage()]))
        self.assertEqual("r__Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales;"
                         " f__Veillonellaceae; g__Mitsuokella; s__Mitsuokella multacida",
                         "; ".join([t.prefix_taxon() for t in taxonomy.get_taxon("s__Mitsuokella multacida").lineage()]))
        self.assertFalse("o__Veillonellales" in taxonomy.hierarchy)
        self.assertFalse("o__Selenomonadales" in taxonomy.hierarchy)
        self.assertFalse("c__Negativicutes" in taxonomy.hierarchy)
        return


class LineageTester(unittest.TestCase):
    def test_verify_rank_occupancy(self):
        from treesapp.entrez_utils import Lineage
        tlin = Lineage()
        incomplete_lin = "d__Archaea; p__Euryarchaeota; c__Methanomicrobia; o__; f__; g__Candidatus Syntrophoarchaeum"
        complete_lin = "d__Domain; p__Phylum; c__Class; o__Order; f__Family"
        no_prefix_lin = "Archaea; Euryarchaeota; Methanomicrobia"

        # Test one: Incomplete lineage should be truncated
        tlin.Lineage = incomplete_lin
        tlin.verify_rank_occupancy()
        self.assertEqual("d__Archaea; p__Euryarchaeota; c__Methanomicrobia", tlin.Lineage)
        # Test two: complete lineage should be the same before and after
        tlin.Lineage = complete_lin
        self.assertEqual(complete_lin, tlin.Lineage)
        # Test three: no prefix
        tlin.Lineage = no_prefix_lin
        self.assertEqual(no_prefix_lin, tlin.Lineage)

    def test_lineage_format_check(self):
        from treesapp.entrez_utils import Lineage
        tlin = Lineage()
        tlin.Lineage = ""
        self.assertFalse(tlin.lineage_format_check())
        tlin.Lineage = "r__Root"
        self.assertFalse(tlin.lineage_format_check())
        tlin.Lineage = "f__Methanomassiliicoccaceae;g__Candidatus Methanoplasma;s__Candidatus Methanoplasma termitum"
        self.assertTrue(tlin.lineage_format_check())
        self.assertEqual(3, len(tlin.Lineage.split(tlin.lin_sep)))
        # Retry now that tlin.Lineage has been reformatted
        self.assertFalse(tlin.lineage_format_check())
        return

    def test_build_lineage(self):
        from treesapp.entrez_utils import Lineage
        test_lin = Lineage()
        lineage = "d__Archaea; p__Euryarchaeota; c__Methanomicrobia; o__Methanophagales"
        test_lin.Domain, test_lin.Phylum, test_lin.Class, test_lin.Order = lineage.split("; ")
        test_lin.build_lineage()
        self.assertEqual(lineage, test_lin.Lineage)
        return


if __name__ == '__main__':
    unittest.main()

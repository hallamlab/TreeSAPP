import os
import unittest
import pytest

from .testing_utils import get_test_data


class TreesappTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp import fasta
        self.test_fasta_data = {"Bad_header_name": "ACGT"}
        self.test_data_file = os.path.join("tests", "file_parser_test.fasta")
        fasta.write_new_fasta(fasta_dict=self.test_fasta_data,
                              fasta_name=self.test_data_file)
        return

    def tearDown(self) -> None:
        if os.path.isfile(self.test_data_file):
            os.remove(self.test_data_file)
        return

    def test_best_discrete_matches(self):
        from treesapp.file_parsers import best_discrete_matches
        from treesapp.hmmer_tbl_parser import HmmMatch
        # Create test HmmMatch instances
        m1 = HmmMatch()
        m2 = HmmMatch()
        m3 = HmmMatch()
        m1.start, m1.end, m1.target_hmm = 80, 200, "NorB"
        m2.start, m2.end, m2.target_hmm = 15, 220, "NorC"
        m3.start, m3.end, m3.target_hmm = 240, 391, "NorB"

        test_matches = [m1, m2, m3]
        best_matches = best_discrete_matches(test_matches)
        self.assertEqual(2, len(best_matches))
        self.assertEqual(3, len(test_matches))
        self.assertTrue(m2 in best_matches and m3 in best_matches)
        return

    def test_read_phylip_to_dict(self):
        from treesapp.file_parsers import read_phylip_to_dict
        from .testing_utils import get_test_data
        with pytest.raises(SystemExit):
            read_phylip_to_dict("test.phy")

        phy_dict = read_phylip_to_dict(get_test_data("test.phy"))
        self.assertEqual(5, len(phy_dict))
        self.assertEqual(76, len(phy_dict["cox2_leita"]))
        return

    def test_read_stockholm_to_dict(self):
        from treesapp.file_parsers import read_stockholm_to_dict
        from .testing_utils import get_test_data
        with pytest.raises(SystemExit):
            read_stockholm_to_dict("test.sto")
        sto_dict = read_stockholm_to_dict(get_test_data("test.sto"))
        self.assertEqual(236, len(sto_dict))
        self.assertEqual(680, len(sto_dict["KKH90701"]))
        return

    def test_load_classifified_sequences_from_assign_output(self):
        from treesapp.file_parsers import load_classified_sequences_from_assign_output
        from treesapp.classy import TreeSAPP
        assigner_instance = TreeSAPP("phylotu")
        refpkg_pquery_map = load_classified_sequences_from_assign_output(get_test_data("marker_test_results"),
                                                                         assigner_instance)
        self.assertEqual(["McrA", "McrB"], list(refpkg_pquery_map.keys()))
        self.assertEqual(13, sum([len(pqueries) for pqueries in refpkg_pquery_map.values()]))
        for rp, pqueries in refpkg_pquery_map.items():
            for pq in pqueries:
                self.assertTrue(len(pq.seq) > 0)
                self.assertEqual(rp, pq.ref_name)

        # Test refpkg filtering
        refpkg_pquery_map = load_classified_sequences_from_assign_output(get_test_data("test_output_TarA"),
                                                                         assigner_instance,
                                                                         refpkg_name="DsrAB")
        self.assertEqual(70, len(refpkg_pquery_map["DsrAB"]))
        for rp, pqueries in refpkg_pquery_map.items():
            for pq in pqueries:
                self.assertEqual("DsrAB", pq.ref_name)

        return

    def test_write_classification_table(self):
        from treesapp.file_parsers import write_classification_table, load_classified_sequences_from_assign_output
        from treesapp.classy import TreeSAPP
        assigner_instance = TreeSAPP("phylotu")
        refpkg_pquery_map = load_classified_sequences_from_assign_output(get_test_data("marker_test_results"),
                                                                         assigner_instance)
        test_output = "tests/tmp_classifications.tsv"
        write_classification_table(tree_saps=refpkg_pquery_map, output_file=test_output, sample_name="Testing")
        with open(test_output) as t_out:
            self.assertEqual(14, len(t_out.readlines()))

        # Try appending
        write_classification_table(tree_saps=refpkg_pquery_map,
                                   sample_name="Testing_testing",
                                   output_file=test_output,
                                   append=True)
        with open(test_output) as t_out:
            self.assertEqual(27, len(t_out.readlines()))

        # Clean up output
        if os.path.isfile(test_output):
            os.remove(test_output)
        return

    def test_read_phenotypes_map(self):
        from treesapp.file_parsers import read_phenotypes
        mcra_metabolism = get_test_data("Mcr_taxonomy-phenotype_map.tsv")

        phenotypes = read_phenotypes(mcra_metabolism)
        self.assertIsInstance(cls=dict, obj=phenotypes)
        self.assertEqual(20, len(phenotypes))
        return

    def test_read_lineage_map(self):
        from treesapp.file_parsers import read_lineage_map
        mock_lineage_map = get_test_data("McrA_lineage_ids - GTDB_map.tsv")
        lineages = read_lineage_map(mock_lineage_map)
        self.assertEqual(149, len(lineages))
        return

    def test_read_lineage_ids(self):
        from treesapp.file_parsers import read_lineage_ids
        mock_lineages = get_test_data("McrA_lineage_ids.tsv")
        lineages = read_lineage_ids(mock_lineages)
        self.assertIsInstance(lineages, dict)
        self.assertEqual(50, len(lineages))
        self.assertEqual([str(x) for x in range(1, 51)], list(lineages.keys()))

        with pytest.raises(SystemExit):
            read_lineage_ids(get_test_data("McrA_lineage_ids - GTDB_map.tsv"))
        return

    def test_check_seq_name_integer_compatibility(self):
        from treesapp import file_parsers
        msa, n_refs = file_parsers.check_seq_name_integer_compatibility(seq_dict=self.test_fasta_data)
        self.assertEqual(-1, n_refs)
        self.assertEqual("Bad_header_name", msa.popitem()[0])
        return

    def test_validate_alignment_trimming(self):
        from treesapp import file_parsers
        from treesapp import fasta
        test_msa = get_test_data("PuhA.mfa")
        headers = set([str(i) + "_PuhA" for i in range(1, 48)])
        tmp_fasta = fasta.FASTA(file_name=test_msa)
        tmp_fasta.load_fasta()

        # Test bad file extension
        with pytest.raises(SystemExit):
            file_parsers.validate_alignment_trimming(["PuhA.stk"], headers)

        # Fail due to a bad sequence name in a fasta
        with pytest.raises(SystemExit):
            file_parsers.validate_alignment_trimming([self.test_data_file], set(self.test_fasta_data))

        # Ensure success
        success, fail, msg = file_parsers.validate_alignment_trimming(msa_files=[test_msa],
                                                                      unique_ref_headers=set(tmp_fasta.get_seq_names()))
        self.assertEqual(32, len(success[test_msa]))
        self.assertEqual([], fail)
        self.assertIsInstance(msg, str)

        # MSA fails due to more sequence names in unique_ref_headers than MSA
        success, fail, msg = file_parsers.validate_alignment_trimming(msa_files=[test_msa],
                                                                      unique_ref_headers=headers)
        self.assertEqual({}, success)
        self.assertEqual([test_msa], fail)
        self.assertIsInstance(msg, str)

        return


if __name__ == '__main__':
    unittest.main()

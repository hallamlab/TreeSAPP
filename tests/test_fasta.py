import unittest


class MyTestCase(unittest.TestCase):
    def test_get_header_format(self):
        from fasta import load_fasta_header_regexes, get_header_format
        header_regexes = load_fasta_header_regexes()
        _, header_db, _ = get_header_format("650716056_50856936 Metig_1237", header_regexes)
        self.assertEqual("unformatted", header_db)
        seq_name = ">1000565.METUNv1_03969"
        _, header_db, _ = get_header_format(seq_name, header_regexes)
        self.assertEqual("eggnog", header_db)
        return

    def test_sequence_info_groups(self):
        from fasta import load_fasta_header_regexes, get_header_format, sequence_info_groups
        header_regexes = load_fasta_header_regexes()
        # Test one: unformatted
        seq_name = "650716056_50856936 Metig_1237"
        header_format_re, header_db, header_molecule = get_header_format(seq_name, header_regexes)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual("650716056_50856936", seq_info_tuple.accession)
        # Test two: NCBI protein
        seq_name = ">AKB49151.1 Methyl coenzyme M reductase alpha subunit [Methanosarcina sp. Kolksee]"
        header_format_re, header_db, header_molecule = get_header_format(seq_name, header_regexes)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual("AKB49151", seq_info_tuple.accession)
        # Test three: EggNOG
        seq_name = ">1000565.METUNv1_03969"
        header_format_re, header_db, header_molecule = get_header_format(seq_name, header_regexes)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual("1000565.METUNv1_03969", seq_info_tuple.accession)
        self.assertEqual("1000565", seq_info_tuple.taxid)
        # Test four: SwissProt, representing accessions with database prefix (e.g. 'sp') with accession enclosed by '|'
        seq_name = 'sp|P06008|RCEH_BLAVI Reaction center protein H chain OS=Blastochloris viridis GN=puhA'
        header_format_re, header_db, header_molecule = get_header_format(seq_name, header_regexes)
        self.assertEqual("sp", header_db)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual("P06008", seq_info_tuple.accession)
        return

    def test_sequence_info_groups_assigned(self):
        from fasta import load_fasta_header_regexes, get_header_format, sequence_info_groups
        header_regexes = load_fasta_header_regexes("PuhA")
        seq_name = 'sp|P06008|RCEH_BLAVI Reaction center protein H chain OS=Blastochloris viridis GN=puhA|PuhA|1_258'
        header_format_re, header_db, header_molecule = get_header_format(seq_name, header_regexes)
        self.assertEqual("ts_assign", header_db)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name, header_regexes)
        self.assertEqual("P06008", seq_info_tuple.accession)
        # self.assertEqual("sp|P06008|RCEH_BLAVI", seq_info_tuple.version)
        return


if __name__ == '__main__':
    unittest.main()

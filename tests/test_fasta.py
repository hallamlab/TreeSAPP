import os
import unittest


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.fasta import load_fasta_header_regexes
        from .testing_utils import get_test_data
        self.header_regexes = load_fasta_header_regexes(code_name="PuhA")

        self.test_fa = get_test_data("McrA_eval.faa")
        self.test_fq = get_test_data("test_TarA.1.fq")
        self.troublesome_fa = get_test_data("fasta_parser_test.fasta")
        return

    def tearDown(self) -> None:
        from glob import glob
        test_fa_prefix, _ = os.path.splitext(os.path.basename(self.test_fa))
        test_fq_prefix, _ = os.path.splitext(os.path.basename(self.test_fq))
        for file_name in sorted(glob(test_fa_prefix + "*") + glob(test_fq_prefix + "*")):
            if os.path.isfile(file_name):
                os.remove(file_name)
        return

    def test_get_header_format(self):
        from treesapp.fasta import get_header_format
        _, header_db, _ = get_header_format("650716056_50856936 Metig_1237", self.header_regexes)
        self.assertEqual("unformatted", header_db)
        seq_name = ">1000565.METUNv1_03969"
        _, header_db, _ = get_header_format(seq_name, self.header_regexes)
        self.assertEqual("eggnog", header_db)
        return

    def test_sequence_info_groups(self):
        from treesapp.fasta import get_header_format
        from treesapp.fasta import sequence_info_groups
        # Test one: unformatted
        seq_name = "650716056_50856936 Metig_1237"
        header_format_re, header_db, _ = get_header_format(seq_name, self.header_regexes)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual("650716056_50856936", seq_info_tuple.accession)
        # Test two: NCBI protein
        seq_name = ">AKB49151.1 Methyl coenzyme M reductase alpha subunit [Methanosarcina sp. Kolksee]"
        header_format_re, header_db, _ = get_header_format(seq_name, self.header_regexes)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual("AKB49151", seq_info_tuple.accession)
        # Test three: EggNOG
        seq_name = ">1000565.METUNv1_03969"
        header_format_re, header_db, _ = get_header_format(seq_name, self.header_regexes)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual("1000565.METUNv1_03969", seq_info_tuple.accession)
        self.assertEqual("1000565", seq_info_tuple.taxid)
        # Test four: SwissProt, representing accessions with database prefix (e.g. 'sp') with accession enclosed by '|'
        seq_name = 'sp|P06008|RCEH_BLAVI Reaction center protein H chain OS=Blastochloris viridis GN=puhA'
        header_format_re, header_db, _ = get_header_format(seq_name, self.header_regexes)
        self.assertEqual("sp", header_db)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual("P06008", seq_info_tuple.accession)
        return

    def test_sequence_info_groups_assigned(self):
        from treesapp.fasta import get_header_format, sequence_info_groups
        seq_name = 'sp|P06008|RCEH_BLAVI Reaction center protein H chain OS=Blastochloris viridis|PuhA|1_258'
        header_format_re, header_db, _ = get_header_format(seq_name, self.header_regexes)
        self.assertEqual("ts_assign", header_db)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name, self.header_regexes)
        self.assertEqual("P06008", seq_info_tuple.accession)
        self.assertEqual("sp|P06008|RCEH_BLAVI", seq_info_tuple.version)
        self.assertEqual("sp|P06008|RCEH_BLAVI Reaction center protein H chain OS=Blastochloris viridis|PuhA|1_258",
                         seq_info_tuple.description)
        return

    def test_read_fasta_to_dict(self):
        from treesapp.fasta import read_fasta_to_dict

        # Test normal functionality with a fasta file
        fa_dict = read_fasta_to_dict(self.test_fa)
        self.assertEqual(236, len(fa_dict))

        # Test reading a file that doesn't exist
        self.assertEqual({}, read_fasta_to_dict("./fake_file.fa"))
        return

    def test_split_fa_size(self):
        from treesapp.fasta import split_fa, read_fasta_to_dict
        split_files = split_fa(fastx=self.test_fa, outdir="./", file_num=4)
        num_split_seqs = sum([len(read_fasta_to_dict(f)) for f in split_files])
        self.assertEqual(len(read_fasta_to_dict(self.test_fa)), num_split_seqs)
        self.assertEqual(4, len(split_files))
        return

    def test_split_fa_count(self):
        from treesapp.fasta import split_fa, read_fasta_to_dict
        split_files = split_fa(fastx=self.test_fa, outdir="./", file_num=1, max_seq_count=236)
        self.assertEqual(1, len(split_files))
        split_files = split_fa(fastx=self.test_fa, outdir="./", file_num=1, max_seq_count=24)
        num_split_seqs = sum([len(read_fasta_to_dict(f)) for f in split_files])
        self.assertEqual(len(read_fasta_to_dict(self.test_fa)), num_split_seqs)
        self.assertEqual(10, len(split_files))
        return

    def test_fq2fa(self):
        from treesapp.fasta import fq2fa
        split_files = fq2fa(fastx=self.test_fq, outdir="./", max_seq_count=6)
        self.assertEqual(2, len(split_files))

    def test_format_read_fasta(self):
        from treesapp.fasta import format_read_fasta
        fasta_dict = format_read_fasta(fasta_input=self.troublesome_fa, molecule="prot", min_seq_length=50)
        self.assertEqual(5, len(fasta_dict))
        return


if __name__ == '__main__':
    unittest.main()

import os
import pytest
import unittest

from .testing_utils import get_test_data


class FastaTester(unittest.TestCase):
    def setUp(self) -> None:
        from treesapp.fasta import load_fasta_header_regexes
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
        seq_name = ">WP_048080940.1"
        header_format_re, header_db, _ = get_header_format(seq_name, self.header_regexes)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual("WP_048080940", seq_info_tuple.accession)
        self.assertEqual("WP_048080940.1", seq_info_tuple.version)

        # Test three: EggNOG
        seq_name = ">1000565.METUNv1_03969"
        header_format_re, header_db, _ = get_header_format(seq_name, self.header_regexes)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual("1000565.METUNv1_03969", seq_info_tuple.accession)
        self.assertEqual("1000565", seq_info_tuple.taxid)
        seq_name = "1148.1006613"
        header_format_re, header_db, _ = get_header_format(seq_name, self.header_regexes)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual(seq_name, seq_info_tuple.version)
        self.assertEqual("1148", seq_info_tuple.taxid)
        seq_name = "4155.Migut.B00751.1.p"
        header_format_re, header_db, _ = get_header_format(seq_name, self.header_regexes)
        self.assertEqual("eggnog", header_db)
        seq_info_tuple = sequence_info_groups(header_format_re.match(seq_name), header_db, seq_name)
        self.assertEqual(seq_name, seq_info_tuple.version)
        self.assertEqual("4155", seq_info_tuple.taxid)

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

    def test_split_fa(self):
        from treesapp.fasta import split_fa, read_fasta_to_dict
        # Test splitting by number of files
        split_files = split_fa(fastx=self.test_fa, outdir="./", file_num=4)
        num_split_seqs = sum([len(read_fasta_to_dict(f)) for f in split_files])
        self.assertEqual(len(read_fasta_to_dict(self.test_fa)), num_split_seqs)
        self.assertEqual(4, len(split_files))

        # Test splitting by count
        split_files = split_fa(fastx=self.test_fa, outdir="./", file_num=1, max_seq_count=236)
        self.assertEqual(1, len(split_files))
        split_files = split_fa(fastx=self.test_fa, outdir="./", file_num=1, max_seq_count=24)
        num_split_seqs = sum([len(read_fasta_to_dict(f)) for f in split_files])
        self.assertEqual(len(read_fasta_to_dict(self.test_fa)), num_split_seqs)
        self.assertEqual(10, len(split_files))
        return

    def test_parameterize_sub_file_size(self):
        from treesapp.fasta import parameterize_sub_file_size
        fsize, nseqs = parameterize_sub_file_size(file_name=self.test_fq, num_files=4, max_seqs=0, record_len=4)
        self.assertEqual(1032, fsize)
        self.assertEqual(0, nseqs)
        fsize, nseqs = parameterize_sub_file_size(file_name=self.test_fa, num_files=0, max_seqs=10)
        self.assertEqual(10, nseqs)
        self.assertEqual(0, fsize)
        return

    def test_fq2fa(self):
        from treesapp.fasta import fq2fa
        split_files = fq2fa(fastx=self.test_fq, outdir="./", file_num=4)
        line_num = 0
        for f in split_files:
            with open(f, 'r') as f_handler:
                line_num += len(f_handler.readlines())
        self.assertEqual(24, line_num)
        self.assertEqual(4, len(split_files))
        split_files = fq2fa(fastx=self.test_fq, outdir="./", max_seq_count=6)
        for f_path in split_files:
            self.assertTrue(os.path.isfile(f_path))
        return

    def test_format_read_fasta(self):
        from treesapp.fasta import format_read_fasta
        fasta_dict = format_read_fasta(fasta_input=self.troublesome_fa, molecule="prot", min_seq_length=50)
        self.assertEqual(5, len(fasta_dict))
        return

    def test_trim_to_length(self):
        from treesapp.fasta import FASTA, register_headers
        _test_length = 10
        mock_fa = FASTA("mock.fa")
        mock_fa.fasta_dict = {"seq_1": "ADJKHJHFAJKHDFKJH", "seq_2": "MCNMVBASUIQRUIYA", "seq_3": "AAAAAAA"}
        mock_fa.header_registry = register_headers(list(mock_fa.fasta_dict.keys()))

        mock_fa.trim_to_length(_test_length)
        seq_lengths = [len(seq) for seq in mock_fa.fasta_dict.values()]
        # Check sequence lengths are correct
        self.assertTrue(min(seq_lengths) == max(seq_lengths) == _test_length)
        # Check FASTA attributes
        self.assertEqual(["seq_1", "seq_2"], mock_fa.get_seq_names())
        self.assertEqual(2, mock_fa.n_seqs())
        self.assertEqual(2, len(mock_fa.header_registry))
        return

    def test_remove_shorter_than(self):
        from treesapp.fasta import FASTA
        test_fa = FASTA(self.test_fa)
        test_fa.load_fasta()
        self.assertEqual(0, test_fa.remove_shorter_than(min_len=30))
        self.assertEqual(1, test_fa.remove_shorter_than(min_len=1000))
        return

    def test_clear(self):
        from treesapp.fasta import FASTA
        test_fa = FASTA(self.test_fa)
        test_fa.load_fasta()
        test_fa.clear()
        self.assertEqual({}, test_fa.fasta_dict)
        return


class HeaderTester(unittest.TestCase):
    def setUp(self) -> None:
        self.eggnog_one = '1148.1006613'
        self.eggnog_two = '28072.Nos7524_3177'
        self.eggnog_three = '59919.PMM0223'

    def test_find_accession(self):
        from treesapp.fasta import Header
        h = Header(self.eggnog_one)
        h.find_accession()
        self.assertEqual('1148.1006613', h.accession)
        self.assertEqual('1148.1006613', h.version)

        h.original = self.eggnog_two
        h.find_accession()
        self.assertEqual(self.eggnog_two, h.version)

        h.original = self.eggnog_three
        h.find_accession()
        self.assertEqual(self.eggnog_three, h.version)

        bin_ex = "SI072_165m_bin.14_k147_157247_1"
        h = Header(bin_ex)
        h.find_accession()
        self.assertEqual("SI072_165m_bin", h.accession)
        self.assertEqual(bin_ex, h.version)
        return

    def test_get_info(self):
        from treesapp.fasta import Header
        h = Header(self.eggnog_one)
        info = h.get_info()
        self.assertEqual(3, len(info.split('\n')))
        return


class FastaUtilitiesTester(unittest.TestCase):
    def setUp(self) -> None:
        self.test_aa_fa = get_test_data("McrA_eval.faa")
        self.test_nuc_fa = get_test_data("marker_test_suite.fna")
        self.troublesome_fa = get_test_data("fasta_parser_test.fasta")
        return

    def test_read_fasta_to_dict(self):
        from treesapp.fasta import read_fasta_to_dict

        # Test normal functionality with a fasta file
        fa_dict = read_fasta_to_dict(self.test_aa_fa)
        self.assertEqual(236, len(fa_dict))

        fa_dict = read_fasta_to_dict(self.test_nuc_fa, num_records=4)
        self.assertEqual(4, len(fa_dict))

        # Test reading a file that doesn't exist
        self.assertEqual({}, read_fasta_to_dict("./fake_file.fa"))
        return

    def test_read_fastq_to_dict(self):
        from treesapp.fasta import read_fastq_to_dict
        # Test a FASTQ file
        fasta_dict = read_fastq_to_dict(fastq_file=get_test_data("SRR3669912_1.fastq"),
                                        num_records=80)
        self.assertEqual(80, len(fasta_dict))
        return

    def test_read_fastx_to_dict(self):
        from treesapp.fasta import read_fastx_to_dict
        # Test file doesn't exist
        self.assertEqual({}, read_fastx_to_dict("test_data" + os.sep + "fake.fa"))
        # Test a file that is neither a FASTA or FASTQ
        with pytest.raises(SystemExit):
            read_fastx_to_dict(fastx=get_test_data("test.sto"))
        # Test a FASTQ file
        fasta_dict = read_fastx_to_dict(fastx=get_test_data("SRR3669912_1.fastq"))
        self.assertEqual(160, len(fasta_dict))
        # Test a FASTA file
        fasta_dict = read_fastx_to_dict(self.test_aa_fa)
        self.assertEqual(236, len(fasta_dict))
        return

    def test_guess_sequence_type(self):
        from treesapp.fasta import guess_sequence_type, read_fasta_to_dict
        # Test bad keyword argument
        with pytest.raises(ValueError):
            guess_sequence_type(fasta=self.test_aa_fa)

        with pytest.raises(SystemExit):
            guess_sequence_type(fastx_file="tests/not_testing.fasta")

        # Test with fasta file for amino acid sequences
        seq_type = guess_sequence_type(fastx_file=self.test_aa_fa)
        self.assertEqual("prot", seq_type)

        # Test with fasta file and nucleotide sequences
        seq_type = guess_sequence_type(fastx_file=self.test_nuc_fa)
        self.assertEqual("dna", seq_type)

        # Test with dictionary and a lot of ambiguity characters
        fa_dict = read_fasta_to_dict(self.troublesome_fa)
        seq_type = guess_sequence_type(fasta_dict=fa_dict)
        self.assertEqual('', seq_type)

        # Test with multiple sequence types
        seq_type = guess_sequence_type(fasta_dict={"s1": "AAAACGGGATCGAGCTACG",
                                                   "s2": "GORVVVLMMMNOPLC"})
        self.assertEqual('', seq_type)
        return


if __name__ == '__main__':
    unittest.main()

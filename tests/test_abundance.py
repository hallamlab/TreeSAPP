import unittest


class AbundanceTester(unittest.TestCase):
    def test_strip_file_to_sample_name(self):
        from treesapp.classy import Abundance
        mock_abund = Abundance()
        mock_abund.strip_file_to_sample_name("path/to/test_TarA.1.fastq")
        self.assertEqual("test_TarA", mock_abund.sample_prefix)

        mock_abund.strip_file_to_sample_name("path/to/test_sample1.1.fq")
        self.assertEqual("test_sample1", mock_abund.sample_prefix)

        mock_abund.strip_file_to_sample_name("path/to/test_sample_R1.fq")
        self.assertEqual("test_sample", mock_abund.sample_prefix)
        return


if __name__ == '__main__':
    unittest.main()

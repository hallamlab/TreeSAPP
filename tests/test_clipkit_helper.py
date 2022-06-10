import os
import unittest

from .testing_utils import get_test_data


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.test_fa = get_test_data('PuhA.mfa')
        self.output_fa = 'PuhA.trim.mfa'

    def tearDown(self) -> None:
        if os.path.isfile(self.output_fa):
            os.remove(self.output_fa)

    def test_run(self):
        from treesapp import clipkit_helper
        from clipkit import modes as ck_modes
        ck = clipkit_helper.ClipKitHelper(fasta_in=self.test_fa,
                                          output_dir='./',
                                          mode="smart-gap")
        ck.run()
        ck.compare_original_and_trimmed_multiple_alignments(min_len=200)
        ck.summarise_trimming()
        self.assertTrue(os.path.isfile(self.output_fa))

        ck.mode = ck_modes.TrimmingMode("kpi-smart-gap")
        ck.run()
        ck.compare_original_and_trimmed_multiple_alignments(min_len=200)
        ck.summarise_trimming()
        self.assertTrue(os.path.isfile(self.output_fa))

        ck.mode = ck_modes.TrimmingMode("gappy")
        ck.run()
        ck.compare_original_and_trimmed_multiple_alignments(min_len=200)
        ck.summarise_trimming()
        self.assertTrue(os.path.isfile(self.output_fa))
        return


if __name__ == '__main__':
    unittest.main()

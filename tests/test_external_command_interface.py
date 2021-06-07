import unittest


class ECITester(unittest.TestCase):
    def setUp(self) -> None:
        self.num_threads = 4
        return

    def test_run_apply_async_multiprocessing(self):
        from treesapp import external_command_interface as ts_eci
        from treesapp import utilities

        res = ts_eci.run_apply_async_multiprocessing(func=utilities.median,
                                                     arguments_list=[[[x]] for x in (range(0, 10))],
                                                     num_processes=self.num_threads,
                                                     pbar_desc="test")
        self.assertEqual(10, len(res))
        self.assertEqual(0, min(res))
        self.assertEqual(9, max(res))
        return


if __name__ == '__main__':
    unittest.main()

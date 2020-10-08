import unittest
import pytest


class TreeSAPPTester(unittest.TestCase):
    def test_main(self):
        from treesapp.__main__ import main
        retcode = main(["treesapp", "info"])
        self.assertEqual(0, retcode)
        with pytest.raises(SystemExit):
            main(["treesapp", "package", "-h"])
        return


if __name__ == '__main__':
    unittest.main()

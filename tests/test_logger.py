import unittest
import os
import pytest


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        import logging
        self.logger = logging.getLogger("treesapp")
        self.test_log_file = "./test_log.txt"
        self.logger.handlers.clear()
        return

    def tearDown(self) -> None:
        if os.path.isfile(self.test_log_file):
            os.remove(self.test_log_file)

    def test_logger_name(self):
        from treesapp import logger
        self.assertEqual("treesapp", logger.logger_name())
        return

    def test_prep_logging(self):
        from treesapp import logger

        logger.prep_logging(log_file=self.test_log_file, verbosity=True)
        self.assertTrue(os.path.isfile(self.test_log_file))

        # Reset logging for tests
        self.logger.handlers.clear()
        # Directory doesn't exist, should fail
        with pytest.raises(SystemExit):
            logger.prep_logging(log_file=os.path.join(os.getcwd(), "must", "fail", "test_log.txt"),
                                verbosity=False)
        return


if __name__ == '__main__':
    unittest.main()

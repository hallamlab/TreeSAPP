import logging
from clipkit import clipkit
from clipkit import args_processing

from treesapp import logger


class ClipKitHelper:
    CLIPKIT_MODES = {"smart-gap"}

    def __init__(self, fasta_in: str, mfa_out: str):
        self.logger = logging.getLogger(logger.logger_name())
        self.input = ""
        self.mfa_out = ""

        self.mode = "smart-gap"
        return

    def run(self):
        # clipkit.execute()
        return

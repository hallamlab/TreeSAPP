import logging
import os.path

from clipkit import clipkit as ck
from clipkit import modes as ck_modes

from treesapp import logger


class ClipKitHelper:
    CLIPKIT_MODES = {"smart-gap"}

    def __init__(self, fasta_in: str, mfa_out=None, mode="smart-gap", gap_prop=0.9):
        self.input = fasta_in
        if mfa_out is None:
            prefix, ext = os.path.splitext(fasta_in)
            self.mfa_out = prefix + ".trim" + ext
        else:
            self.mfa_out = mfa_out

        self.logger = logging.getLogger(logger.logger_name())
        self.mode = ck_modes.TrimmingMode(mode)
        self.gap_prop = gap_prop

        self.ff_in = "fasta"
        self.ff_out = "fasta"
        return

    def run(self):

        ck.execute(input_file=self.input,
                   input_file_format=self.ff_in,
                   output_file=self.mfa_out,
                   output_file_format=self.ff_out,
                   gaps=self.gap_prop,
                   complement=False,
                   mode=self.mode,
                   use_log=False)
        return

    def summarise_trimming(self):
        return

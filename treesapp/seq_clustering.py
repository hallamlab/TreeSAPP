import logging
import sys


class Cluster:
    def __init__(self, rep_name):
        self.representative = rep_name
        self.members = list()
        self.lca = ''

    def get_info(self):
        info_string = "Representative: " + str(self.representative) + "\n" + \
                      "LCA: " + self.lca + "\n" + \
                      "Members:\n\t" + "\n\t".join([', '.join(member) for member in self.members]) + "\n"
        return info_string


class BlastAln:
    def __init__(self):
        self.subject = ""
        self.query = ""
        self.start = 0
        self.end = 0
        self.alnlen = 0
        self.gapopen = 0
        self.mismatch = 0
        self.qstart = 0
        self.qend = 0
        self.tstart = 0
        self.tend = 0
        self.evalue = 0.0
        self.pident = 0.0
        self.bits = 0.0

    def load_blast_tab(self, line: str, sep="\t") -> None:
        try:
            fields = line.strip().split(sep)
        except ValueError:
            logging.error("Unable to parse line in from alignment table:\n{}\n".format(line))
            sys.exit(7)

        try:
            self.subject, self.query = fields[0], fields[1]
            self.pident, self.alnlen, self.mismatch, self.gapopen = fields[2:6]
            self.qstart, self.qend, self.tstart, self.tend = fields[6:10]
            self.evalue, self.bits = fields[10], fields[11]
        except ValueError:
            logging.error("Incorrect format for line in alignment file. Twelve were expected, found {}.\n{}\n"
                          "".format(len(fields), line))
            sys.exit(7)

        return
__author__ = 'Connor Morgan-Lang'

import argparse
import logging
import re
import sys
from treesapp.classy import prep_logging


class StockholmSeq:
    def __init__(self, name):
        self.id = name
        self.accession = ""
        self.desc = ""
        self.seq_name = ""

    def tabularize(self):
        fields = [self.id, self.accession, self.seq_name, self.desc]
        return "\t".join(fields)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sto_file", required=True,
                        help="Path to the Stockholm file")
    parser.add_argument("-o", "--tbl_out", required=True,
                        help="Path to the output table")
    args = parser.parse_args()
    return args


def read_stockholm(sto_file: str) -> dict:
    curr_id = ""
    curr_acc = ""
    curr_desc = ""
    sto_dict = dict()
    sto_feature_re = re.compile(r"^#=GF\s([A-Z]{2})\s+(.*)$")
    sto_seq_re = re.compile(r"^#=GS\s(.*)\s+.*$")
    logging.info("Reading stockholm file... ")
    with open(sto_file, 'r', encoding="latin-1") as sto_handler:
        for line in sto_handler:
            if not re.match(r"^#=.*", line):
                curr_id = ""
                curr_acc = ""
                curr_desc = ""
            elif sto_feature_re.match(line):
                key, value = sto_feature_re.match(line).groups()
                if key == "ID":
                    curr_id = value
                elif key == "AC":
                    curr_acc = value
                elif key == "DE":
                    curr_desc = value
                else:
                    pass
            elif sto_seq_re.match(line):
                seq_name = sto_seq_re.match(line).group(1)
                sto_inst = StockholmSeq(curr_id)
                if not curr_id:
                    logging.error("Family ID not set.\n")
                    sys.exit(3)
                sto_inst.desc = curr_desc
                sto_inst.accession = curr_acc
                sto_inst.seq_name = seq_name
                try:
                    sto_dict[curr_id].append(sto_inst)
                except KeyError:
                    sto_dict[curr_id] = [sto_inst]
            else:
                pass
    logging.info("done.\n")
    return sto_dict


def write_sto_table(sto_dict: dict, table_file: str) -> None:
    str_buffer = ""
    str_size = 0
    logging.info("Writing Stockholm table... ")
    with open(table_file, 'w') as tbl_handler:
        for sto_id in sto_dict:
            for sto_inst in sto_dict[sto_id]:  # type: StockholmSeq
                sto_tabs = sto_inst.tabularize()
                str_buffer += sto_tabs + "\n"
                str_size += len(sto_tabs) + 1
                if str_size >= 10E6:
                    tbl_handler.write(str_buffer)
                    str_buffer = ""
                    str_size = 0
        tbl_handler.write(str_buffer)
    logging.info("done.\n")
    return


def main():
    args = get_arguments()
    prep_logging()
    sto_dict = read_stockholm(args.sto_file)
    write_sto_table(sto_dict, args.tbl_out)
    return


main()

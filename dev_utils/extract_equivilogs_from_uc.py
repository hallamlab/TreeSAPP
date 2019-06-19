__author__ = 'Connor Morgan-Lang'


import argparse
import logging
from treesapp.classy import prep_logging
from treesapp.file_parsers import read_uc
from treesapp.classy import Cluster


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--uc", required=True,
                        help="Path to a UCLUST cluster (uc) file")
    parser.add_argument("-s", "--seq_names", required=True,
                        help="Path to a file listing the reference sequence names")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file name, mapping each reference sequence to a list of equivalogs")
    args = parser.parse_args()
    return args


def read_seq_names_list(names_file: str) -> set:
    names_list = set()
    with open(names_file, 'r') as names_handler:
        for line in names_handler:
            if line[0] == '>':
                line = line[1:]
            names_list.add(line.strip())
    return names_list


def get_equivalogs(uc_dict: dict, ref_seq_names: set) -> dict:
    equivalogs = dict()
    for num_id in sorted(uc_dict, key=int):
        cluster = uc_dict[num_id]  # type: Cluster
        member_set = set([n[0] for n in cluster.members])
        member_set.add(cluster.representative)
        shared_names = ref_seq_names.intersection(member_set)
        if len(shared_names) >= 1:
            non_ref_seqs = member_set.difference(ref_seq_names)
            for ref_name in shared_names:
                equivalogs[ref_name] = non_ref_seqs
    return equivalogs


def write_equivalogs(equivalogs: dict, output_file: str) -> None:
    with open(output_file, 'w') as tbl_handler:
        for ref_name in equivalogs:
            tbl_handler.write(ref_name + "\t" + ','.join(equivalogs[ref_name]) + "\n")
    return


def main():
    args = get_arguments()
    prep_logging()
    uc_dict = read_uc(args.uc)
    ref_seq_names = read_seq_names_list(args.seq_names)
    equivalogs = get_equivalogs(uc_dict, ref_seq_names)
    write_equivalogs(equivalogs, args.output)
    return


main()

#!/usr/bin/env python3

import argparse
import sys
import os
import glob

__author__ = 'Connor Morgan-Lang'


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--table_list", required=False, dest="list",
                        help="Text file with list of table to merge, each on a new line.")
    parser.add_argument("-d", "--table_dir", required=False, dest="dir",
                        help="Directory containing only tables to be merged.")
    parser.add_argument("-s", "--separator", default="\t", required=False, dest="sep",
                        help="The field separator used in table to merge. There can only be one! [ DEFAULT = '\\t' ]")
    parser.add_argument("-o", "--output_table", required=True, dest="output",
                        help="Name of the file to write merged tables to.")
    args = parser.parse_args()

    if not (args.list or args.dir):
        parser.error("No tables provided. Please provide either a list (-l) or directory (-d) ot tables.")
        sys.exit(3)
    return args


def read_table_list(list_file):
    tables = []
    with open(list_file) as list_handler:
        for line in list_handler:
            file_path = line.strip()
            if not os.path.isfile(file_path):
                sys.exit("ERROR: " + file_path + " doesn't exist!")
            tables.append(file_path)
    return tables


def fetch_tables_from_dir(table_dir):
    return glob.glob(table_dir + os.sep + "*")


def merge_tables(table_files: list, sep: str, output_table: str):
    master_field_names = []
    final_fields = {}  # Dictionary mapping field names to filed positions for the final output
    data_aggregate = []  # Each list element is a string
    # Identify the superset of fields from all tables
    for tf in table_files:
        tf_handle = open(tf, 'r')
        tf_header = tf_handle.readline().strip().split(sep=sep)
        tf_handle.close()
        for field in tf_header:
            if field not in master_field_names:
                master_field_names.append(field)

    i = 0
    while i < len(master_field_names):
        final_fields[master_field_names[i]] = i
        i += 1
    master_field_names.clear()

    try:
        ot_handler = open(output_table, 'w')
    except IOError:
        sys.exit("Unable to open " + output_table + " for writing.")
    ot_handler.write(sep.join([name for name, pos in sorted(final_fields.items(), key=lambda x: x[1])]) + "\n")

    # Read in data, filling missing columns with 'NA'
    n = 0
    merged_dat = []
    header_map = dict()  # Maps field position to field name, e.g. {0: Sample, 1: Marker: 2: Length}
    for tf in sorted(table_files):
        with open(tf, 'r') as tf_handle:
            tf_header = tf_handle.readline().strip().split(sep=sep)
            while n < len(tf_header):
                header_map[tf_header[n]] = n
                n += 1
            for line in tf_handle:
                data = line.strip().split(sep=sep)
                for name, pos in sorted(final_fields.items(), key=lambda x: x[1]):
                    if name in header_map:
                        merged_dat.append(data[header_map[name]])
                    else:
                        merged_dat.append("NA")
                data_aggregate.append(sep.join(merged_dat))
                merged_dat.clear()
        header_map.clear()
        n = 0
        # Write data to output_table
        ot_handler.write("\n".join(data_aggregate) + "\n")
        data_aggregate.clear()

    ot_handler.close()
    return


def main():
    args = get_options()
    tables = []
    if args.list:
        tables += read_table_list(args.list)
    if args.dir:
        tables += fetch_tables_from_dir(args.dir)
    if not tables:
        sys.exit("ERROR: No tables were found for merging.")
    merge_tables(tables, args.sep, args.output)


main()

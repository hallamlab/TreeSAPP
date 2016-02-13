#!/usr/bin/python

import os
import subprocess
import sys

__author__ = 'Connor Morgan-Lang'


def main():
    # TODO: get this to work
    # Nucleotide test run
    nuc_test_cmd = "./mltreemap.py -i test_data/dna_test.fasta -o nuc_testOut -T 4 --verbose"
    nuc_proc = subprocess.Popen(nuc_test_cmd,
                                shell=True,
                                stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE)
    return_sig_nuc = nuc_proc.returncode
    out, err = nuc_proc.communicate()
    with open("nuc_testOut/test.stdout") as nuc_output:
        nuc_output.write(out)
    with open("nuc_testOut/test.stderr") as nuc_err:
        nuc_err.write(err)
    if return_sig_nuc > 0:
        print "Something went wrong during the DNA input test - look in nuc_testOut/test.stderr for details"
        sys.exit(1)

    # Protein test run
    prot_test_cmd = "./mltreemap.py -i test_data/protein_test.fasta -o prot_testOut -r a -T 4 --verbose"
    prot_proc = subprocess.Popen(prot_test_cmd,
                                 shell=True,
                                 stderr=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
    out, err = prot_proc.communicate()
    return_sig_prot = prot_proc.returncode
    #subprocess.Popen.wait()
    with open("prot_testOut/test.stdout") as prot_output:
        prot_output.write(out)
    with open("prot_testOut/test.stderr") as prot_err:
        prot_err.write(err)
    if return_sig_prot > 0:
        print "Something went wrong during the protein input test - look in prot_testOut/test.stderr for details"
        sys.exit(1)
    os.rmdir("prot_testOut/")

main()

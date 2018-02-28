__author__ = 'Connor Morgan-Lang'

import os
import re
import sys
import Bio
from Bio import Entrez
from urllib import error

from external_command_interface import launch_write_command


def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path_element in os.environ["PATH"].split(os.pathsep):
            path_element = path_element.strip('"')
            exe_file = os.path.join(path_element, program)
            if is_exe(exe_file):
                return exe_file
    return None


def os_type():
    """Return the operating system of the user."""
    x = sys.platform
    if x:

        hits = re.search(r'darwin', x, re.I)
        if hits:
            return 'mac'

        hits = re.search(r'win', x, re.I)
        if hits:
            return 'win'

        hits = re.search(r'linux', x, re.I)
        if hits:
            return 'linux'


def find_executables(args):
    """
    Finds the executables in a user's path to alleviate the requirement of a sub_binaries directory
    :param args: command-line arguments objects
    :return: exec_paths beings the absolute path to each executable
    """
    exec_paths = dict()
    dependencies = ["blastn", "blastx", "blastp", "genewise", "Gblocks", "hmmalign", "raxmlHPC",
                    "trimal", "cmalign", "cmsearch", "makeblastdb", "muscle", "hmmbuild", "cmbuild"]

    # Extra executables necessary for certain modes of TreeSAPP
    try:
        if args.rpkm:
            dependencies += ["bwa", "rpkm"]

        if args.update_tree:
            dependencies.append("usearch")

        if args.consensus:
            dependencies.append("FGS+")
    except AttributeError:
        pass

    if os_type() == "linux":
        args.executables = args.treesapp + "sub_binaries" + os.sep + "ubuntu"
    if os_type() == "mac":
        args.executables = args.treesapp + "sub_binaries" + os.sep + "mac"
    elif os_type() == "win" or os_type() is None:
        sys.exit("ERROR: Unsupported OS")

    for dep in dependencies:
        if is_exe(args.executables + os.sep + dep):
            exec_paths[dep] = str(args.executables + os.sep + dep)
        # For rpkm and potentially other executables that are compiled ad hoc
        elif is_exe(args.treesapp + "sub_binaries" + os.sep + dep):
            exec_paths[dep] = str(args.treesapp + "sub_binaries" + os.sep + dep)
        elif which(dep):
            exec_paths[dep] = which(dep)
        else:
            sys.stderr.write("Could not find a valid executable for " + dep + ". ")
            sys.exit("Bailing out.")

    args.executables = exec_paths
    return args


def reformat_string(string):
    if string and string[0] == '>':
        header = True
    else:
        header = False
    string = re.sub("\[|\]|\(|\)|\/|\\\\|'|<|>", '', string)
    if header:
        string = '>' + string
    string = re.sub("\s|;|,|\|", '_', string)
    if len(string) > 110:
        string = string[0:109]
    while string and string[-1] == '.':
        string = string[:-1]
    return string


class Autovivify(dict):
    """In cases of Autovivify objects, enable the referencing of variables (and sub-variables)
    without explicitly declaring those variables beforehand."""

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def query_entrez_taxonomy(search_term):
    handle = Entrez.esearch(db="Taxonomy",
                            term=search_term,
                            retmode="xml")
    record = Entrez.read(handle)
    try:
        org_id = record["IdList"][0]
        if org_id:
            handle = Entrez.efetch(db="Taxonomy", id=org_id, retmode="xml")
            records = Entrez.read(handle)
            lineage = str(records[0]["Lineage"])
        else:
            return
    except IndexError:
        if 'QueryTranslation' in record.keys():
            # If 'QueryTranslation' is returned, use it for the final Entrez query
            lineage = record['QueryTranslation']
            lineage = re.sub("\[All Names\].*", '', lineage)
            lineage = re.sub('[()]', '', lineage)
            for word in lineage.split(' '):
                handle = Entrez.esearch(db="Taxonomy", term=word, retmode="xml")
                record = Entrez.read(handle)
                try:
                    org_id = record["IdList"][0]
                except IndexError:
                    continue
                handle = Entrez.efetch(db="Taxonomy", id=org_id, retmode="xml")
                records = Entrez.read(handle)
                lineage = str(records[0]["Lineage"])
                if re.search("cellular organisms", lineage):
                    break
        else:
            sys.stderr.write("ERROR: Unable to handle record returned by Entrez.efetch!\n")
            sys.stderr.write("Database = Taxonomy\n")
            sys.stderr.write("term = " + search_term + "\n")
            sys.stderr.write("record = " + str(record) + "\n")
            raise IndexError

    return lineage


def get_lineage(search_term, molecule_type):
    """
    Used to return the NCBI taxonomic lineage of the sequence
    :param: search_term: The NCBI search_term
    :param: molecule_type: "dna", "rrna", "prot", or "tax - parsed from command line arguments
    :return: string representing the taxonomic lineage
    """
    # TODO: fix potential error PermissionError:
    # [Errno 13] Permission denied: '/home/connor/.config/biopython/Bio/Entrez/XSDs'
    # Fixed with `sudo chmod 777 .config/biopython/Bio/Entrez/`
    if not search_term:
        raise AssertionError("ERROR: search_term for Entrez query is empty!\n")
    if float(Bio.__version__) < 1.68:
        # This is required due to a bug in earlier versions returning a URLError
        raise AssertionError("ERROR: version of biopython needs to be >=1.68! " +
                             str(Bio.__version__) + " is currently installed. Exiting now...")
    Entrez.email = "c.morganlang@gmail.com"
    # Test the internet connection:
    try:
        Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
    except error.URLError:
        raise AssertionError("ERROR: Unable to serve Entrez query. Are you connected to the internet?")

    # Determine which database to search using the `molecule_type`
    if molecule_type == "dna" or molecule_type == "rrna" or molecule_type == "ambig":
        database = "nucleotide"
    elif molecule_type == "prot":
        database = "protein"
    elif molecule_type == "tax":
        database = "Taxonomy"
    else:
        sys.stderr.write("Welp. We're not sure how but the molecule type is not recognized!\n")
        sys.stderr.write("Please create an issue on the GitHub page.")
        sys.exit(8)

    # Find the lineage from the search_term ID
    lineage = ""
    ncbi_sequence_databases = ["nucleotide", "protein"]
    handle = None
    if database in ["nucleotide", "protein"]:
        try:
            handle = Entrez.efetch(db=database, id=str(search_term), retmode="xml")
        except error.HTTPError:
            # if molecule_type == "ambig":
                x = 0
                while handle is None and x < len(ncbi_sequence_databases):
                    backup_db = ncbi_sequence_databases[x]
                    if backup_db != database:
                        try:
                            handle = Entrez.efetch(db=backup_db, id=str(search_term), retmode="xml")
                        except error.HTTPError:
                            handle = None
                    x += 1
                if handle is None:
                    sys.stderr.write("\nERROR: Bad Entrez.efetch request and all back-up searches failed for '" +
                                     str(search_term) + "'\n")
                    sys.exit(99)
            # else:
            #     sys.stderr.write("\nWARNING: Bad Entrez.efetch request:\n"
            #                      "id='" + str(search_term) + "' does not exist in database='" + database + "'\n")
        try:
            record = Entrez.read(handle)
        except UnboundLocalError:
            raise UnboundLocalError
        if len(record) >= 1:
            try:
                if "GBSeq_organism" in record[0]:
                    organism = record[0]["GBSeq_organism"]
                    # To prevent Entrez.efectch from getting confused by non-alphanumeric characters:
                    organism = re.sub('[)(\[\]]', '', organism)
                    lineage = query_entrez_taxonomy(organism)
            except IndexError:
                for word in record['QueryTranslation']:
                    lineage = query_entrez_taxonomy(word)
                    print(lineage)
        else:
            # Lineage is already set to "". Just return and move on to the next attempt
            pass
    else:
        try:
            print("Searching taxonomy database for " + search_term)
            lineage = query_entrez_taxonomy(search_term)
        except UnboundLocalError:
            sys.stderr.write("WARNING: Unable to find Entrez taxonomy using organism name:\n\t")
            sys.stderr.write(search_term + "\n")

    return lineage


def remove_dashes_from_msa(fasta_in, fasta_out):
    """
    fasta_out is the new FASTA file written with no dashes (unaligned)
    There are no line breaks in this file, whereas there may have been in fasta_in
    :param fasta_in: Multiply-aligned FASTA file
    :param fasta_out: FASTA file to write
    :return:
    """
    dashed_fasta = open(fasta_in, 'r')
    fasta = open(fasta_out, 'w')
    sequence = ""

    line = dashed_fasta.readline()
    while line:
        if line[0] == '>':
            if sequence:
                fasta.write(sequence + "\n")
                sequence = ""
            fasta.write(line)
        else:
            sequence += re.sub('[-.]', '', line.strip())
        line = dashed_fasta.readline()
    fasta.write(sequence + "\n")
    dashed_fasta.close()
    fasta.close()
    return


def generate_blast_database(args, fasta, molecule, prefix, multiple=True):
    """

    :param args:
    :param fasta: File to make a BLAST database for
    :param molecule: 'prot' or 'nucl' - necessary argument for makeblastdb
    :param prefix: prefix string for the output BLAST database
    :param multiple: Flag indicating the input `fasta` is a MSA. Alignment information is removed prior to makeblastdb
    :return:
    """

    # Remove the multiple alignment information from fasta_replaced_file and write to fasta_mltree
    blastdb_out = prefix + ".fa"
    if multiple:
        if blastdb_out == fasta:
            sys.stderr.write("ERROR: prefix.fa is the same as " + fasta + " and would be overwritten!\n")
            sys.exit(11)
        remove_dashes_from_msa(fasta, blastdb_out)
        blastdb_in = blastdb_out
    else:
        blastdb_in = fasta

    sys.stdout.write("Making the BLAST database for " + blastdb_in + "... ")

    # Format the `makeblastdb` command
    makeblastdb_command = [args.executables["makeblastdb"]]
    makeblastdb_command += ["-in", blastdb_in]
    makeblastdb_command += ["-out", blastdb_out]
    makeblastdb_command += ["-input_type", "fasta"]
    makeblastdb_command += ["-dbtype", molecule]

    # Launch the command
    stdout, makeblastdb_pro_returncode = launch_write_command(makeblastdb_command)

    sys.stdout.write("done\n")
    sys.stdout.flush()

    return stdout, blastdb_out


def clean_lineage_string(lineage):
    bad_strings = ["cellular organisms; ", "delta/epsilon subdivisions; ", "\(miscellaneous\)"]
    for bs in bad_strings:
        lineage = re.sub(bs, '', lineage)
    # filter 'group'
    if re.search('group; ', lineage):
        reconstructed_lineage = ""
        ranks = lineage.split("; ")
        for rank in ranks:
            if not re.search("group$", rank):
                reconstructed_lineage = reconstructed_lineage + str(rank) + '; '
        reconstructed_lineage = re.sub('; $', '', reconstructed_lineage)
        lineage = reconstructed_lineage
    return lineage

__author__ = 'Connor Morgan-Lang'

import os
import re
import sys
import Bio
from Bio import Entrez
from urllib import error
from HMMER_domainTblParser import filter_incomplete_hits, filter_poor_hits, DomainTableParser, format_split_alignments
import time

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
    dependencies = ["prodigal", "hmmbuild", "hmmalign", "hmmsearch", "raxmlHPC", "trimal"]
    # old_dependencies = ["blastn", "blastx", "blastp", "genewise", "Gblocks", "makeblastdb", "muscle"]

    # Extra executables necessary for certain modes of TreeSAPP
    if hasattr(args, "rpkm") and args.rpkm:
        dependencies += ["bwa", "rpkm"]

    if hasattr(args, "update_tree"):
        if args.update_tree:
            dependencies.append("usearch")

    if hasattr(args, "cluster") or hasattr(args, "multiple_alignment") or hasattr(args, "fast"):
        if args.cluster:
            dependencies.append("usearch")
        dependencies.append("mafft")
        if args.fast:
            dependencies.append("FastTree")

    if args.molecule == "rrna":
        dependencies += ["cmalign", "cmsearch", "cmbuild"]

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


def multiple_query_entrez_taxonomy(search_term_set):
    """
    Function for submitting multiple queries using Entrez.efetch to the 'Taxonomy' database.
    :param search_term_set: Inputs are a set of organism names (based off their accession records)
    :return: A dictionary mapping each of the unique organism names in search_term_set to a full taxonomic lineage
    """
    search_term_result_map = dict()
    for search_term in search_term_set:
        search_term_result_map[search_term] = query_entrez_taxonomy(search_term)
    return search_term_result_map


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


def prep_for_entrez_query():
    Entrez.email = "c.morganlang@gmail.com"
    Entrez.tool = "treesapp"
    # Test the internet connection:
    try:
        Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
    except error.URLError:
        raise AssertionError("ERROR: Unable to serve Entrez query. Are you connected to the internet?")
    return


def parse_accessions_from_entrez_xml(record):
    accession = ""
    versioned = ""
    accession_keys = ["GBSeq_locus", "GBSeq_primary-accession"]
    version_keys = ["GBInterval_accession", "GBSeq_accession-version"]
    for accession_key in accession_keys:
        if accession_key in record:
            accession = record[accession_key]
            break
    for version_key in version_keys:
        if version_key in record:
            versioned = record[version_key]
            break
    return accession, versioned


def parse_organism_from_entrez_xml(record):
    organism = ""
    if len(record) >= 1:
        try:
            if "GBSeq_organism" in record:
                organism = record["GBSeq_organism"]
                # To prevent Entrez.efectch from getting confused by non-alphanumeric characters:
                organism = re.sub('[)(\[\]]', '', organism)
        except IndexError:
            sys.stderr.write("WARNING: 'GBSeq_organism' not found in Entrez record.\n")
    else:
        pass
    return organism


def parse_lineage_from_record(record):
    lineage = ""
    if len(record) >= 1:
        try:
            if "GBSeq_organism" in record:
                organism = record["GBSeq_organism"]
                # To prevent Entrez.efectch from getting confused by non-alphanumeric characters:
                organism = re.sub('[)(\[\]]', '', organism)
                lineage = query_entrez_taxonomy(organism)
        except IndexError:
            sys.stderr.write("WARNING: 'GBSeq_organism' not found in Entrez record.\n")
            for word in record['QueryTranslation']:
                lineage = query_entrez_taxonomy(word)
                print(lineage)
    else:
        # Lineage is already set to "". Just return and move on to the next attempt
        pass
    return lineage


def xml_parser(xml_record, term):
    """
    Recursive function for parsing individual xml records
    :param xml_record:
    :param term:
    :return:
    """
    # TODO: Finish this off - would be great for consistently extracting data from xml
    value = None
    if type(xml_record) == str:
        return value
    if term not in xml_record.keys():
        for record in xml_record:
            value = xml_parser(record, term)
            if value:
                return value
            else:
                continue
    else:
        return xml_record[term]
    return value


def get_multiple_lineages(search_term_list, molecule_type, log_file_handler):
    """

    :param search_term_list:
    :param molecule_type: "dna", "rrna", "prot", or "tax - parsed from command line arguments
    :param log_file_handler: A file handler object for the log
    :return: A dictionary mapping accession IDs (keys) to organisms and lineages (values)
    """
    accession_lineage_map = dict()
    all_accessions = set()
    if not search_term_list:
        raise AssertionError("ERROR: search_term for Entrez query is empty!\n")
    if float(Bio.__version__) < 1.68:
        # This is required due to a bug in earlier versions returning a URLError
        raise AssertionError("ERROR: version of biopython needs to be >=1.68! " +
                             str(Bio.__version__) + " is currently installed. Exiting now...")

    # Do some semi-important stuff
    prep_for_entrez_query()

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

    sys.stdout.write("Retrieving Entrez " + database + " records for each reference sequence... ")
    sys.stdout.flush()

    # Must be cautious with this first query since some accessions are not in the Entrez database anymore
    # and return with `urllib.error.HTTPError: HTTP Error 502: Bad Gateway`
    master_records = []
    chunk_size = 60
    log_file_handler.write("\nEntrez.efetch query time for sequences (minutes:seconds):\n")
    for i in range(0, len(search_term_list), chunk_size):
        start_time = time.time()
        chunk = search_term_list[i:i+chunk_size]
        try:
            handle = Entrez.efetch(db=database, id=','.join([str(sid) for sid in chunk]), retmode="xml")
            # for sid in chunk:
            #     handle = Entrez.efetch(db=database, id=sid, retmode="xml")
            master_records += Entrez.read(handle)
        except:
            log_file_handler.write("WARNING: Unable to parse XML data from Entrez.efetch! "
                                   "It is either potentially corrupted or cannot be found in the database.\n")
            log_file_handler.write("Offending accessions from this batch:\n")
            for sid in chunk:
                try:
                    handle = Entrez.efetch(db=database, id=sid, retmode="xml")
                    record = Entrez.read(handle)
                    master_records.append(record[0])
                except:
                    log_file_handler.write("\t" + str(sid) + "\n")
        end_time = time.time()
        hours, remainder = divmod(end_time - start_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        log_file_handler.write("\t" + str(i) + ' - ' + str(i+chunk_size) + "\t" +
                               ':'.join([str(minutes), str(round(seconds, 2))]) + "\n")

    sys.stdout.write("done.\n")
    log_file_handler.write("\n")
    sys.stdout.write("Retrieving lineage information for each sequence from Entrez... ")
    sys.stdout.flush()

    unique_organisms = set()
    # Instantiate the master_records for linking each organism to accessions, and empty fields
    for record in master_records:
        accession, versioned = parse_accessions_from_entrez_xml(record)
        accession_lineage_map[(accession, versioned)] = dict()
        accession_lineage_map[(accession, versioned)]["organism"] = parse_organism_from_entrez_xml(record)
        accession_lineage_map[(accession, versioned)]["lineage"] = ""
        all_accessions.update([accession, versioned])

    for tuple_key in accession_lineage_map.keys():
        unique_organisms.add(accession_lineage_map[tuple_key]["organism"])

    organism_lineage_map = multiple_query_entrez_taxonomy(unique_organisms)
    for tuple_key in accession_lineage_map:
        organism_name = accession_lineage_map[tuple_key]["organism"]
        accession_lineage_map[tuple_key]["lineage"] = organism_lineage_map[organism_name]

    sys.stdout.write("done.\n")

    return accession_lineage_map, all_accessions


def check_lineage(lineage, organism_name):
    """
    Sometimes the NCBI lineage is incomplete.
    Currently, this function uses organism_name to ideally add Species to the lineage
    :param lineage: A semi-colon separated taxonomic lineage
    :param organism_name: Name of the organism. Parsed from the sequence header (usually at the end in square brackets)
    :return: A string with lineage information
    """
    proper_species_re = re.compile("^[A-Z][a-z]+ [a-z]+$")
    if proper_species_re.match(lineage.split("; ")[-1]):
        return lineage
    elif len(lineage.split("; ")) == 7 and proper_species_re.match(organism_name):
        return lineage + "; " + organism_name
    else:
        return lineage


def get_lineage_robust(reference_sequence_list, molecule):
    failed_reference_sequences = list()

    for reference_sequence in reference_sequence_list:
        strikes = 0
        lineage = ""
        while strikes < 3:
            if strikes == 0:
                if reference_sequence.accession:
                    lineage = get_lineage(reference_sequence.accession, molecule)
                else:
                    sys.stderr.write("WARNING: no accession available for Entrez query:\n")
                    reference_sequence.get_info()
                if type(lineage) is str and len(lineage) > 0:
                    # The query was successful
                    strikes = 3
            elif strikes == 1:
                # Unable to determine lineage from the search_term provided,
                # try to parse organism name from description
                if reference_sequence.organism:
                    try:
                        taxon = ' '.join(reference_sequence.organism.split('_')[:2])
                    except IndexError:
                        taxon = reference_sequence.organism
                    lineage = get_lineage(taxon, "tax")
                    if type(lineage) is str and len(lineage) > 0:
                        # The query was successful
                        # try:
                        #     lineage += '; ' + reference_sequence.organism.split('_')[-2]
                        # except IndexError:
                        #     lineage += '; ' + reference_sequence.organism
                        strikes = 3
                else:
                    # Organism information is not available, time to bail
                    strikes += 1
            elif strikes == 2:
                lineage = get_lineage(lineage, "tax")
            strikes += 1
        if not lineage:
            sys.stderr.write("\nWARNING: Unable to find lineage for sequence with following data:\n")
            reference_sequence.get_info()
            lineage = ""
        # TODO: test this
        if reference_sequence.organism:
            lineage = check_lineage(lineage, reference_sequence.organism)
        else:
            reference_sequence.organism = reference_sequence.description
        reference_sequence.lineage = lineage
        failed_reference_sequences.append(reference_sequence)
    return failed_reference_sequences


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
    Entrez.tool = "treesapp"
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
                    # sys.stderr.write("\nWARNING: Bad Entrez.efetch request and all back-up searches failed for '" +
                    #                  str(search_term) + "'\n")
                    return lineage
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
            # sys.stderr.write("WARNING: Searching taxonomy database for '" + search_term + "'\n")
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
    non_standard_names_re = re.compile("group; | cluster; ")
    bad_strings = ["cellular organisms; ", "delta/epsilon subdivisions; ", "\(miscellaneous\)"]
    for bs in bad_strings:
        lineage = re.sub(bs, '', lineage)
    # filter 'group' and 'cluster'
    if non_standard_names_re.search(lineage):
        reconstructed_lineage = ""
        ranks = lineage.split("; ")
        for rank in ranks:
            if not (re.search("group$", rank) or re.search("cluster$", rank)):
                reconstructed_lineage = reconstructed_lineage + str(rank) + '; '
        reconstructed_lineage = re.sub('; $', '', reconstructed_lineage)
        lineage = reconstructed_lineage
    return lineage


def best_match(matches):
    """
    Function for finding the best alignment in a list of HmmMatch() objects
    The best match is based off of the full sequence score
    :param matches: A list of HmmMatch() objects
    :return: The best HmmMatch
    """
    # TODO: Incorporate the alignment intervals to allow for proteins with multiple different functional domains
    # Code currently only permits multi-domains of the same gene
    best_target_hmm = ""
    best_alignment = None
    top_score = 0
    for match in matches:
        # match.print_info()
        if match.full_score > top_score:
            best_alignment = match
            best_target_hmm = match.target_hmm
            top_score = match.full_score
    return best_target_hmm, best_alignment


def parse_domain_tables(args, hmm_domtbl_files, log=None):
    # Check if the HMM filtering thresholds have been set
    if not hasattr(args, "min_e"):
        args.min_e = 0.01
        args.min_acc = 0.6
        args.perc_aligned = 80
    # Print some stuff to inform the user what they're running and what thresholds are being used.
    if args.verbose:
        sys.stdout.write("Filtering HMM alignments using the following thresholds:\n")
        sys.stdout.write("\tMinimum E-value = " + str(args.min_e) + "\n")
        sys.stdout.write("\tMinimum acc = " + str(args.min_acc) + "\n")
        sys.stdout.write("\tMinimum percentage of the HMM covered = " + str(args.perc_aligned) + "%\n")
    sys.stdout.write("Parsing domain tables generated by HMM searches for high-quality matches... ")

    raw_alignments = 0
    seqs_identified = 0
    dropped = 0
    fragmented = 0
    glued = 0
    multi_alignments = 0  # matches of the same query to a different HMM (>1 lines)
    hmm_matches = dict()
    orf_gene_map = dict()

    # TODO: Capture multimatches across multiple domain table files
    for domtbl_file in hmm_domtbl_files:
        rp_marker, reference = re.sub("_domtbl.txt", '', os.path.basename(domtbl_file)).split("_to_")
        domain_table = DomainTableParser(domtbl_file)
        domain_table.read_domtbl_lines()
        distinct_matches, fragmented, glued, multi_alignments, raw_alignments = format_split_alignments(domain_table,
                                                                                                        fragmented,
                                                                                                        glued,
                                                                                                        multi_alignments,
                                                                                                        raw_alignments)
        purified_matches, dropped = filter_poor_hits(args, distinct_matches, dropped)
        complete_gene_hits, dropped = filter_incomplete_hits(args, purified_matches, dropped)
        # for match in complete_gene_hits:
        #     match.genome = reference
        #     if match.target_hmm not in hmm_matches.keys():
        #         hmm_matches[match.target_hmm] = list()
        #     hmm_matches[match.target_hmm].append(match)
        #     seqs_identified += 1

        for match in complete_gene_hits:
            match.genome = reference
            if match.orf not in orf_gene_map:
                orf_gene_map[match.orf] = dict()
            orf_gene_map[match.orf][match.target_hmm] = match
            if match.target_hmm not in hmm_matches.keys():
                hmm_matches[match.target_hmm] = list()

    for orf in orf_gene_map:
        if len(orf_gene_map[orf]) == 1:
            target_hmm = list(orf_gene_map[orf].keys())[0]
            hmm_matches[target_hmm].append(orf_gene_map[orf][target_hmm])
        else:
            optional_matches = [orf_gene_map[orf][target_hmm] for target_hmm in orf_gene_map[orf]]
            target_hmm, match = best_match(optional_matches)
            hmm_matches[target_hmm].append(match)
            multi_alignments += 1
            dropped += (len(optional_matches) - 1)

            if log:
                dropped_annotations = list()
                for optional in optional_matches:
                    if optional.target_hmm != target_hmm:
                        dropped_annotations.append(optional.target_hmm)
                log.write("HMM search annotations for " + orf + ":\n")
                log.write("\tRetained\t" + target_hmm + "\n")
                log.write("\tDropped\t" + ','.join(dropped_annotations) + "\n")

        seqs_identified += 1

    sys.stdout.write("done.\n")

    if seqs_identified == 0 and dropped == 0:
        sys.stderr.write("\tWARNING: No alignments found!\n")
        sys.stderr.write("TreeSAPP is exiting now.\n")
        sys.exit(11)
    if seqs_identified == 0 and dropped > 0:
        sys.stderr.write("\tWARNING: No alignments met the quality cut-offs!\n")
        sys.stderr.write("TreeSAPP is exiting now.\n")
        sys.exit(13)

    sys.stdout.write("\tNumber of markers identified:\n")
    for marker in sorted(hmm_matches):
        sys.stdout.write("\t\t" + marker + "\t" + str(len(hmm_matches[marker])) + "\n")
        # For debugging:
        # for match in hmm_matches[marker]:
        #     match.print_info()
    if args.verbose:
        sys.stdout.write("\tInitial alignments:\t" + str(raw_alignments) + "\n")
        sys.stdout.write("\tAlignments discarded:\t" + str(dropped) + "\n")
        sys.stdout.write("\tFragmented alignments:\t" + str(fragmented) + "\n")
        sys.stdout.write("\tAlignments scaffolded:\t" + str(glued) + "\n")
        sys.stdout.write("\tMulti-alignments:\t" + str(multi_alignments) + "\n")
        sys.stdout.write("\tSequences identified:\t" + str(seqs_identified) + "\n")

    sys.stdout.flush()
    return hmm_matches


def median(lst):
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0


def read_colours_file(args, annotation_file):
    """
    Read annotation data from 'annotation_file' and store it in marker_subgroups under the appropriate
    marker and data_type.
    :param args:
    :param annotation_file:
    :return: A dictionary of lists where each list is populated by tuples with start and end leaves
    """
    try:
        style_handler = open(annotation_file, 'r')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + annotation_file + " for reading!\n")
        sys.exit()

    clusters = dict()
    field_sep = ''
    internal_nodes = True

    line = style_handler.readline()
    # Skip the header
    while line.strip() != "DATA":
        header_fields = line.strip().split(' ')
        if header_fields[0] == "SEPARATOR":
            if header_fields[1] == "SPACE":
                field_sep = ' '
            elif header_fields[1] == "TAB":
                field_sep = '\t'
            else:
                sys.stderr.write("ERROR: Unknown separator used in " + annotation_file + ": " + header_fields[1] + "\n")
                sys.stderr.flush()
                sys.exit()
        line = style_handler.readline()
    # For RGB
    range_line_rgb = re.compile("^(\d+)\|(\d+)" + re.escape(field_sep) +
                                "range" + re.escape(field_sep) +
                                ".*\)" + re.escape(field_sep) +
                                "(.*)$")
    single_node_rgb = re.compile("^(\d+)" + re.escape(field_sep) +
                                 "range" + re.escape(field_sep) +
                                 ".*\)" + re.escape(field_sep) +
                                 "(.*)$")
    lone_node_rgb = re.compile("^(.*)" + re.escape(field_sep) +
                               "range" + re.escape(field_sep) +
                               ".*\)" + re.escape(field_sep) +
                               "(.*)$")

    # For hexadecimal
    range_line = re.compile("^(\d+)\|(\d+)" + re.escape(field_sep) +
                            "range" + re.escape(field_sep) +
                            "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                            "(.*)$")
    single_node = re.compile("^(\d+)" + re.escape(field_sep) +
                             "range" + re.escape(field_sep) +
                             "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                             "(.*)$")
    lone_node = re.compile("^(.*)" + re.escape(field_sep) +
                           "range" + re.escape(field_sep) +
                           "#[0-9A-Za-z]{6}" + re.escape(field_sep) +
                           "(.*)$")

    # Begin parsing the data from 4 columns
    line = style_handler.readline().strip()
    while line:
        if range_line.match(line):
            style_data = range_line.match(line)
            start, end, description = style_data.groups()
            internal_nodes = False
        elif range_line_rgb.match(line):
            style_data = range_line_rgb.match(line)
            start, end, description = style_data.groups()
            internal_nodes = False
        elif single_node.match(line):
            style_data = single_node.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
        elif single_node_rgb.match(line):
            style_data = single_node_rgb.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
        elif lone_node.match(line):
            style_data = lone_node.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
        elif lone_node_rgb.match(line):
            style_data = lone_node_rgb.match(line)
            start, end, description = style_data.group(1), style_data.group(1), style_data.group(2)
        else:
            sys.stderr.write("ERROR: Unrecognized line formatting in " + annotation_file + ":\n")
            sys.stderr.write(line + "\n")
            sys.exit()

        description = style_data.groups()[-1]
        if description not in clusters.keys():
            clusters[description] = list()
        clusters[description].append((start, end))

        line = style_handler.readline().strip()

    style_handler.close()

    if args.verbose:
        sys.stdout.write("\tParsed " + str(len(clusters)) +
                         " clades from " + annotation_file + "\n")

    return clusters, internal_nodes


def convert_outer_to_inner_nodes(clusters, internal_node_map):
    leaf_annotation_map = dict()
    for cluster in clusters.keys():
        if cluster not in leaf_annotation_map:
            leaf_annotation_map[cluster] = list()
        for frond_tips in clusters[cluster]:
            start, end = frond_tips
            # Find the minimum set that includes both start and end
            warm_front = dict()
            # Add all the potential internal nodes
            for inode in internal_node_map:
                clade = internal_node_map[inode]
                if start in clade:
                    warm_front[inode] = clade
            for inode in sorted(warm_front, key=lambda x: len(warm_front[x])):
                if end in warm_front[inode]:
                    leaf_annotation_map[cluster].append(inode)
                    break
    return leaf_annotation_map


def annotate_internal_nodes(args, internal_node_map, clusters):
    """
    A function for mapping the clusters to all internal nodes of the tree.
    It also adds overlapping functional annotations for deep internal nodes and ensures all the leaves are annotated.
    :param args:
    :param internal_node_map: A dictionary mapping the internal nodes (keys) to the leaf nodes (values)
    :param clusters: Dictionary with the cluster names for keys and a tuple containing leaf boundaries as values
    :return: A dictionary of the annotation (AKA group) as keys and internal nodes as values
    """
    annotated_clade_members = dict()
    leaf_group_members = dict()
    leaves_in_clusters = set()

    # Create a dictionary to map the cluster name (e.g. Function, Activity, Class, etc) to all the leaf nodes
    for annotation in clusters.keys():
        if annotation not in annotated_clade_members:
            annotated_clade_members[annotation] = set()
        if annotation not in leaf_group_members:
            leaf_group_members[annotation] = set()
        for i_node in internal_node_map:
            if i_node in clusters[annotation]:
                for leaf in internal_node_map[i_node]:
                    leaf_group_members[annotation].add(leaf)
                    leaves_in_clusters.add(leaf)
        # Find the set of internal nodes that are children of this annotated clade
        for i_node in internal_node_map:
            if leaf_group_members[annotation].issuperset(internal_node_map[i_node]):
                annotated_clade_members[annotation].add(i_node)

    if args.verbose:
        sys.stdout.write("\tCaptured " + str(len(leaves_in_clusters)) + " nodes in clusters.\n")

    return annotated_clade_members, leaves_in_clusters

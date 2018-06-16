__author__ = 'Connor Morgan-Lang'

import sys
import time
import re
import Bio
from Bio import Entrez
from urllib import error


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


def prep_for_entrez_query():
    sys.stdout.write("Preparing Bio.Entrez for NCBI queries... ")
    sys.stdout.flush()
    Entrez.email = "c.morganlang@gmail.com"
    Entrez.tool = "treesapp"
    # Test the internet connection:
    try:
        Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
    except error.URLError:
        raise AssertionError("ERROR: Unable to serve Entrez query. Are you connected to the internet?")
    sys.stdout.write("done.\n")
    return


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
    log_file_handler.write("\nEntrez.efetch query time for accessions (minutes:seconds):\n")
    for i in range(0, len(search_term_list), chunk_size):
        start_time = time.time()
        chunk = search_term_list[i:i+chunk_size]
        try:
            handle = Entrez.efetch(db=database, id=','.join([str(sid) for sid in chunk]), retmode="xml")
            # for sid in chunk:
            #     handle = Entrez.efetch(db=database, id=sid, retmode="xml")
            master_records += Entrez.read(handle)
        # Broad exception clause but THE NUMBER OF POSSIBLE ERRORS IS TOO DAMN HIGH!
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

    start_time = time.time()
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

    end_time = time.time()
    hours, remainder = divmod(end_time - start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    log_file_handler.write("Entrez.efetch query time for lineages (minutes:seconds): ")
    log_file_handler.write(':'.join([str(minutes), str(round(seconds, 2))]) + "\n\n")

    return accession_lineage_map, all_accessions


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


def get_lineage_robust(reference_sequence_list, molecule):
    accession_lineage_map = dict()

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
        accession_lineage_map[reference_sequence.accession] = lineage
    return accession_lineage_map


def verify_lineage_information(accession_lineage_map, all_accessions, fasta_record_objects,
                               taxa_searched, molecule, log_file_handle):
    """
    Function used for parsing records returned by Bio.Entrez.efetch queries and identifying inconsistencies
    between the search terms and the results
    :param accession_lineage_map: A dictionary mapping accession.versionID tuples to taxonomic lineages
    :param taxa_searched: An integer for tracking number of accessions queried (currently number of lineages provided)
    :param molecule: Type of molecule (prot, dna, rrna) used for choosing the Entrez database to query
    :param log_file_handle: A handle for the log file for recording warnings and stats
    :return:
    """
    failed_accession_queries = list()
    if (len(accession_lineage_map.keys()) + taxa_searched) != len(fasta_record_objects):
        # Records were not returned for all sequences. Time to figure out which ones!
        log_file_handle.write("WARNING: Entrez did not return a record for every accession queried.\n")
        log_file_handle.write("Don't worry, though. We'll figure out which ones are missing.\n")
    log_file_handle.write("Entrez.efetch query stats:\n")
    log_file_handle.write("\tDownloaded\t" + str(len(accession_lineage_map.keys())) + "\n")
    log_file_handle.write("\tProvided\t" + str(taxa_searched) + "\n")
    log_file_handle.write("\tTotal\t\t" + str(len(fasta_record_objects)) + "\n\n")

    # Find the lineage searches that failed, add lineages to reference_sequences that were successfully identified
    unambiguous_accession_lineage_map = dict()
    for mltree_id_key in fasta_record_objects.keys():
        reference_sequence = fasta_record_objects[mltree_id_key]
        if not reference_sequence.lineage:
            if reference_sequence.accession in all_accessions:
                taxa_searched += 1
                for tuple_key in accession_lineage_map:
                    accession, versioned = tuple_key
                    if reference_sequence.accession == accession or reference_sequence.accession == versioned:
                        if accession_lineage_map[tuple_key]["lineage"] == "":
                            failed_accession_queries.append(reference_sequence)
                        else:
                            # The query was successful! Add it and increment
                            unambiguous_accession_lineage_map[reference_sequence.accession] = accession_lineage_map[tuple_key]["lineage"]
            else:
                failed_accession_queries.append(reference_sequence)
    # For debugging:
    # print("Currently searched:", taxa_searched)

    # Attempt to find appropriate lineages for the failed accessions (e.g. using organism name as search term)
    # Failing this, lineages will be set to "Unclassified"
    if len(failed_accession_queries) > 0:
        log_file_handle.write("Missed records:\n")
        accession_lineage_map = get_lineage_robust(failed_accession_queries, molecule)
        for reference_sequence in failed_accession_queries:
            lineage = accession_lineage_map[reference_sequence.accession]
            if lineage == "":
                log_file_handle.write("WARNING: Unable to determine the taxonomic lineage for " +
                                      reference_sequence.accession + "\n")
                lineage = "Unclassified"
            taxa_searched += 1
            unambiguous_accession_lineage_map[reference_sequence.accession] = lineage
            log_file_handle.write("Accession=" + reference_sequence.accession + "\t")
            log_file_handle.write("Lineage=" + lineage + "\n")

    if taxa_searched < len(fasta_record_objects.keys()):
        sys.stderr.write("ERROR: Not all sequences (" + str(taxa_searched) + '/'
                         + str(len(fasta_record_objects)) + ") were queried against the NCBI taxonomy database!\n")
        sys.exit(22)

    return fasta_record_objects, unambiguous_accession_lineage_map


def write_accession_lineage_map(mapping_file, accession_lineage_map):
    """
    Function for writing a map of NCBI accession IDs to their respective taxonomic lineages
     using a list of ReferenceSequence objects
    :param mapping_file: Name of a file to write these data
    :param accession_lineage_map: A dictionary mapping accessions to lineages
    :return:
    """
    try:
        map_file_handler = open(mapping_file, 'w')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + mapping_file, " for writing!\n")
        sys.exit()

    for accession in accession_lineage_map:
        map_file_handler.write(accession + "\t" + accession_lineage_map[accession] + "\n")

    map_file_handler.close()
    return


def read_accession_taxa_map(mapping_file):
    """
    A function for reading intermediate files made by write_accession_lineage_map to avoid the time-consuming download

    :param mapping_file:
    :return: accession_lineage_map:  A dictionary mapping accessions to lineages
    """
    try:
        map_file_handler = open(mapping_file, 'r')
    except IOError:
        sys.stderr.write("ERROR: Unable to open " + mapping_file, " for reading!\n")
        sys.exit()

    accession_lineage_map = dict()
    for line in map_file_handler:
        accession, lineage = line.strip().split("\t")
        if accession not in accession_lineage_map:
            accession_lineage_map[accession] = lineage
        else:
            raise AssertionError("ERROR: " + accession + " present in " + mapping_file + " multiple times!")

    map_file_handler.close()
    return accession_lineage_map


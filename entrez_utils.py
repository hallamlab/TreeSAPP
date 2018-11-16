__author__ = 'Connor Morgan-Lang'

import sys
import time
import re
import Bio
import logging
from Bio import Entrez
from urllib import error


def multiple_query_entrez_taxonomy(search_term_set):
    """
    Function for submitting multiple queries using Entrez.efetch to the 'Taxonomy' database.

    :param search_term_set: Inputs are a set of organism names (based off their accession records)
    :return: A dictionary mapping each of the unique organism names in search_term_set to a full taxonomic lineage
    """
    search_term_result_map = dict()
    # TODO: Query Entrez server with multiple tax_ids at once, parse records to map the org_id to tax_id
    for search_term in search_term_set:
        search_term_result_map[search_term] = query_entrez_taxonomy(search_term)
    # TODO: Query Entrez server with multiple org_ids at once, parse records to map the org_id to lineage
    return search_term_result_map


def query_entrez_taxonomy(search_term):
    lineage = ""
    try:
        handle = Entrez.esearch(db="Taxonomy",
                                term=search_term,
                                retmode="xml")
    except error.HTTPError:
        logging.warning("Unable to find the taxonomy for " + search_term + "\n")
        return lineage
    record = Entrez.read(handle)
    try:
        org_id = record["IdList"][0]
        if org_id:
            try:
                handle = Entrez.efetch(db="Taxonomy", id=org_id, retmode="xml")
                records = Entrez.read(handle)
                lineage = str(records[0]["Lineage"])
            except error.HTTPError:
                return lineage
        else:
            return lineage
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
    if not lineage:
        logging.warning("Unable to handle record returned by Entrez.efetch!\n" +
                        "Database = Taxonomy\n" +
                        "term = " + search_term + "\n" +
                        "record = " + str(record) + "\n")
    return lineage


def parse_accessions_from_entrez_xml(record):
    accession = ""
    versioned = ""
    alternatives = list()
    accession_keys = ["GBSeq_locus", "GBSeq_primary-accession"]
    version_keys = ["GBInterval_accession", "GBSeq_accession-version"]
    alternate_keys = ["GBSeq_other-seqids"]
    for accession_key in accession_keys:
        if accession_key in record:
            accession = record[accession_key]
            break
    for version_key in version_keys:
        if version_key in record:
            versioned = record[version_key]
            break
    for alt_key in alternate_keys:
        if alt_key in record:
            for alt in record[alt_key]:
                try:
                    alternatives.append(re.search("\|+(.*)$", alt).group(1))
                except AttributeError:
                    logging.debug("Unable to parse alternative accession from string: '" + str(record[alt_key]) + "'\n")
    return accession, versioned, alternatives


def parse_organism_from_entrez_xml(record):
    organism = ""
    if len(record) >= 1:
        try:
            if "GBSeq_organism" in record:
                organism = record["GBSeq_organism"]
                # To prevent Entrez.efectch from getting confused by non-alphanumeric characters:
                organism = re.sub('[)(\[\]]', '', organism)
        except IndexError:
            logging.warning("'GBSeq_organism' not found in Entrez record.\n")
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
            logging.warning("'GBSeq_organism' not found in Entrez record.\n" +
                            "\n".join([query_entrez_taxonomy(word) for word in record['QueryTranslation']]))
    else:
        # Lineage is already set to "". Just return and move on to the next attempt
        pass
    return lineage


def prep_for_entrez_query():
    logging.info("Preparing Bio.Entrez for NCBI queries... ")
    Entrez.email = "c.morganlang@gmail.com"
    Entrez.tool = "treesapp"
    # Test the internet connection:
    try:
        Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
    except error.URLError:
        logging.error("Unable to serve Entrez query. Are you connected to the internet?")
    logging.info("done.\n")
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


def query_entrez_accessions(search_term_list, molecule_type):
    """

    :param search_term_list:
    :param molecule_type: "dna", "rrna", "prot", or "tax - parsed from command line arguments
    :return: A dictionary mapping accession IDs (keys) to organisms and lineages (values)
    """

    if float(Bio.__version__) < 1.68:
        # This is required due to a bug in earlier versions returning a URLError
        logging.error("Version of biopython needs to be >=1.68! " +
                      str(Bio.__version__) + " is currently installed.\n")
        sys.exit(9)

    # Determine which database to search using the `molecule_type`
    if molecule_type == "dna" or molecule_type == "rrna" or molecule_type == "ambig":
        database = "nucleotide"
    elif molecule_type == "prot":
        database = "protein"
    elif molecule_type == "tax":
        database = "Taxonomy"
    else:
        logging.error("Welp. We're not sure how but the molecule type is not recognized!\n" +
                      "Please create an issue on the GitHub page.")
        sys.exit(9)

    # Must be cautious with this first query since some accessions are not in the Entrez database anymore
    # and return with `urllib.error.HTTPError: HTTP Error 502: Bad Gateway`
    master_records = list()
    durations = list()
    chunk_size = 90

    for i in range(0, len(search_term_list), chunk_size):
        start_time = time.time()
        chunk = search_term_list[i:i+chunk_size]
        try:
            handle = Entrez.efetch(db=database, id=','.join([str(sid) for sid in chunk]), retmode="xml")
            master_records += Entrez.read(handle)
        # Broad exception clause but THE NUMBER OF POSSIBLE ERRORS IS TOO DAMN HIGH!
        except:
            bad_sids = list()
            for sid in chunk:
                try:
                    handle = Entrez.efetch(db=database, id=sid, retmode="xml")
                    record = Entrez.read(handle)
                    master_records.append(record[0])
                except:
                    bad_sids.append("\t" + str(sid))

            logging.warning("Unable to parse XML data from Entrez.efetch! "
                            "Either the XML is corrupted or the query terms cannot be found in the database.\n"
                            "Offending accessions from this batch:\n" + "\n".join(bad_sids) + "\n")
        end_time = time.time()
        hours, remainder = divmod(end_time - start_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        durations.append(str(i) + ' - ' + str(i+chunk_size) + "\t" + ':'.join([str(minutes), str(round(seconds, 2))]))

    logging.debug("Entrez.efetch query time for accessions (minutes:seconds):\n\t" +
                  "\n\t".join(durations) + "\n")
    return master_records


def get_multiple_lineages(search_term_list: list, molecule_type: str):
    if not search_term_list:
        logging.error("Search_term for Entrez query is empty\n")
        sys.exit(9)

    # Do some semi-important stuff
    prep_for_entrez_query()

    logging.info("Retrieving Entrez records for each reference sequence... ")
    accession_lineage_map = dict()
    all_accessions = set()
    unique_organisms = set()
    updated_accessions = dict()
    attempt = 1
    while attempt < 3:
        if len(search_term_list) == 0:
            break

        logging.debug("ATTEMPT " + str(attempt) + ':' +
                      "\n\tNumber of search terms = " + str(len(search_term_list)) + "\n")
        records_batch = query_entrez_accessions(search_term_list, molecule_type)

        # Instantiate the master_records for linking each organism to accessions, and empty fields
        for record in records_batch:
            accession, versioned, alt = parse_accessions_from_entrez_xml(record)
            if not accession and not versioned:
                logging.debug("Neither accession nor accession.version parsed from:\n" + str(record) + "\n")
                continue
            elif accession not in search_term_list and versioned not in search_term_list:
                for alt_key in alt:
                    if alt_key in search_term_list:
                        versioned = alt_key
                    updated_accessions[alt_key] = (accession, versioned)
            accession_lineage_map[(accession, versioned)] = dict()
            accession_lineage_map[(accession, versioned)]["organism"] = parse_organism_from_entrez_xml(record)
            accession_lineage_map[(accession, versioned)]["lineage"] = ""
            all_accessions.update([accession, versioned])

        # Tolerance for failed Entrez accession queries
        rescued = 0
        success = 0
        x = 0
        while x < len(search_term_list):
            if search_term_list[x] in all_accessions:
                search_term_list.pop(x)
                success += 1
            elif search_term_list[x] in updated_accessions:
                search_term_list.pop(x)
                rescued += 1
            else:
                x += 1
        logging.debug("ATTEMPT " + str(attempt) + ":" +
                      "\n\tQueries mapped ideally = " + str(success) +
                      "\n\tQueries mapped with alternative accessions = " + str(rescued) + "\n")
        attempt += 1

    for tuple_key in accession_lineage_map.keys():
        unique_organisms.add(accession_lineage_map[tuple_key]["organism"])

    logging.info("done.\n")

    logging.info("Retrieving lineage information for each sequence from Entrez... ")
    start_time = time.time()
    organism_lineage_map = multiple_query_entrez_taxonomy(unique_organisms)
    for tuple_key in accession_lineage_map:
        accession_lineage_map[tuple_key]["lineage"] = organism_lineage_map[accession_lineage_map[tuple_key]["organism"]]

    logging.info("done.\n")

    end_time = time.time()
    hours, remainder = divmod(end_time - start_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.debug("Entrez.efetch query time for lineages (minutes:seconds):\n" +
                  ':'.join([str(minutes), str(round(seconds, 2))]) + "\n")

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
        logging.error("Welp. We're not sure how but the molecule type is not recognized!\n" +
                      "Please create an issue on the GitHub page.")
        sys.exit(9)

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
            lineage = query_entrez_taxonomy(search_term)
        except UnboundLocalError:
            logging.warning("Unable to find Entrez taxonomy using organism name:\n\t" + search_term + "\n")

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
                    logging.warning("No accession available for Entrez query:\n" +
                                    reference_sequence.get_info())
                if isinstance(lineage, str) and len(lineage) > 0:
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
            logging.warning("Unable to find lineage for sequence with following data:\n" +
                            reference_sequence.get_info())
            lineage = ""
        # TODO: test this
        if reference_sequence.organism:
            lineage = check_lineage(lineage, reference_sequence.organism)
        else:
            reference_sequence.organism = reference_sequence.description
        accession_lineage_map[reference_sequence.accession] = lineage
    return accession_lineage_map


def verify_lineage_information(accession_lineage_map, all_accessions, fasta_record_objects, taxa_searched, molecule):
    """
    Function used for parsing records returned by Bio.Entrez.efetch queries and identifying inconsistencies
    between the search terms and the results

    :param accession_lineage_map: A dictionary mapping accession.versionID tuples to taxonomic lineages
    :param all_accessions:
    :param fasta_record_objects:
    :param taxa_searched: An integer for tracking number of accessions queried (currently number of lineages provided)
    :param molecule: Type of molecule (prot, dna, rrna) used for choosing the Entrez database to query
    :return:
    """
    failed_accession_queries = list()
    if (len(accession_lineage_map.keys()) + taxa_searched) != len(fasta_record_objects):
        # Records were not returned for all sequences. Time to figure out which ones!
        logging.warning("Entrez did not return a record for every accession queried.\n"
                        "Don't worry, though. We'll figure out which ones are missing.\n")
    logging.debug("Entrez.efetch query stats:\n"
                  "\tDownloaded\t" + str(len(accession_lineage_map.keys())) + "\n" +
                  "\tProvided\t" + str(taxa_searched) + "\n" +
                  "\tTotal\t\t" + str(len(fasta_record_objects)) + "\n\n")

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
                        if not reference_sequence.organism and accession_lineage_map[tuple_key]["organism"]:
                            reference_sequence.organism = accession_lineage_map[tuple_key]["organism"]
            else:
                failed_accession_queries.append(reference_sequence)
        else:
            unambiguous_accession_lineage_map[reference_sequence.accession] = reference_sequence.lineage

    # Attempt to find appropriate lineages for the failed accessions (e.g. using organism name as search term)
    # Failing this, lineages will be set to "Unclassified"
    if len(failed_accession_queries) > 0:
        misses_strings = list()
        accession_lineage_map = get_lineage_robust(failed_accession_queries, molecule)
        for reference_sequence in failed_accession_queries:
            lineage = accession_lineage_map[reference_sequence.accession]
            if lineage == "":
                logging.warning("Unable to determine the taxonomic lineage for " +
                                reference_sequence.accession + "\n")
                lineage = "Unclassified"
            taxa_searched += 1
            unambiguous_accession_lineage_map[reference_sequence.accession] = lineage
            misses_strings.append("\tAccession=" + reference_sequence.accession + ", " + "Lineage=" + lineage)
        logging.debug("Recovered records:\n" + "\n".join(misses_strings) + "\n")

    if taxa_searched < len(fasta_record_objects.keys()):
        logging.error("Not all sequences (" + str(taxa_searched) + '/'
                      + str(len(fasta_record_objects)) + ") were queried against the NCBI taxonomy database!\n")
        sys.exit(9)

    # Add the species designation since it is often not included in the sequence record's lineage
    proper_species_re = re.compile("^[A-Z][a-z]+ [a-z]+$")
    for treesapp_id in fasta_record_objects:
        ref_seq = fasta_record_objects[treesapp_id]
        if not ref_seq.lineage:
            try:
                lineage = unambiguous_accession_lineage_map[ref_seq.accession]
            except KeyError:
                logging.error("Lineage information was not retrieved for " + ref_seq.accession + "!\n" +
                              "Please remove the output directory and restart.\n")
                sys.exit(13)
            lr = lineage.split("; ")
            if len(lr) == 7 and proper_species_re.match(ref_seq.organism):
                unambiguous_accession_lineage_map[ref_seq.accession] = lineage + "; " + ref_seq.organism
            elif ref_seq.organism not in lr and len(lr) <= 6 and re.match("^[A-Z][a-z]+$", ref_seq.organism):
                unambiguous_accession_lineage_map[ref_seq.accession] = lineage + "; " + ref_seq.organism
            # print(','.join([lineage, "organism: " + ref_seq.organism,
            #                 "\n", "Final: " + unambiguous_accession_lineage_map[ref_seq.accession]]))

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
        logging.error("Unable to open " + mapping_file, " for writing!\n")
        sys.exit(9)

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
    except OSError:
        logging.error("Unable to open " + mapping_file, " for reading!\n")
        sys.exit(9)

    accession_lineage_map = dict()
    for line in map_file_handler:
        accession, lineage = line.strip().split("\t")
        if accession not in accession_lineage_map:
            accession_lineage_map[accession] = str(lineage)
        else:
            logging.error(accession + " present in " + mapping_file + " multiple times!")
            sys.exit(9)

    map_file_handler.close()
    return accession_lineage_map


def build_entrez_queries(fasta_record_objects: dict):
    """
    Function to create data collections to fulfill entrez query searches

    :param fasta_record_objects: A list of ReferenceSequence objects - lineage information to be filled
    :return: Set containing unique accessions to query Entrez
    """
    num_lineages_provided = 0
    entrez_query_list = set()
    unavailable = list()
    for num_id in fasta_record_objects:
        ref_seq = fasta_record_objects[num_id]
        # Only need to download the lineage information for those sequences that don't have it encoded in their header
        if ref_seq.lineage:
            num_lineages_provided += 1
        else:
            if ref_seq.accession:
                entrez_query_list.add(ref_seq.accession)
            else:
                unavailable.append(ref_seq.description)
    if len(unavailable) > 0:
        logging.warning("Neither accession nor lineage available for:\n\t" +
                        "\n\t".join(unavailable))
    return list(entrez_query_list), num_lineages_provided


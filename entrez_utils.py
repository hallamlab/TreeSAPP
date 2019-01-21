__author__ = 'Connor Morgan-Lang'

import sys
import time
import re
import Bio
import logging
from Bio import Entrez
from urllib import error


def validate_target_db(db_type: str):
    """
    Takes a `db_type` string and matches it with the appropriate database name to be used in an Entrez query
    :param db_type: Molecule or database type
    :return: Proper Entrez database name
    """
    # Determine which database to search using the `db_type`
    if db_type == "dna" or db_type == "rrna" or db_type == "ambig":
        database = "nucleotide"
    elif db_type == "prot":
        database = "protein"
    elif db_type == "tax":
        database = "Taxonomy"
    else:
        logging.error("Welp. We're not sure how but the molecule type is not recognized!\n" +
                      "Please create an issue on the GitHub page.")
        sys.exit(9)

    return database


def tolerant_entrez_query(search_term_list: list, db="Taxonomy", method="fetch", retmode="xml", chunk_size=100):
    """
    Function for performing Entrez-database queries using BioPython's Entrez utilities.
    It is able to break up the complete list of search terms,
     perform the search and if any chunks fail send individual queries for each item in the sub_list.
    :param search_term_list: A list of GenBank accessions, NCBI taxonomy IDs, or organism names
    :param db: Name of the Entrez database to query
    :param method: Either fetch or search corresponding to Entrez.efetch and Entrez.esearch, respectively
    :param retmode: The format of the Entrez records to be returned ('fasta' or 'xml', typically)
    :param chunk_size: Size of the sub_lists for each Entrez query. Fewer than 100 is recommended.
    :return: A list of Entrez records, in the format specific by `retmode`
    """
    read_records = list()
    durations = list()
    failures = list()

    # Check the database name
    if db not in ["nucleotide", "protein", "Taxonomy"]:
        logging.error("Unknown Entrez database '" + db + "'.\n")
        sys.exit(9)

    # Create the sub-lists from `search_term_list` of length `chunk_size`
    for i in range(0, len(search_term_list), chunk_size):
        start_time = time.time()
        sub_list = search_term_list[i:i + chunk_size]
        try:
            if method == "fetch":
                handle = Entrez.efetch(db=db, id=','.join([str(sid) for sid in sub_list]), retmode=retmode)
            else:
                handle = Entrez.esearch(db=db, term=','.join([str(sid) for sid in sub_list]))
            if chunk_size > 1:
                read_records += Entrez.read(handle)
            else:
                read_records.append(Entrez.read(handle))
        # Broad exception clause but THE NUMBER OF POSSIBLE ERRORS IS TOO DAMN HIGH!
        # Try sending individual queries instead to determine problematic query IDs
        except:
            for sid in sub_list:
                try:
                    if method == "fetch":
                        handle = Entrez.efetch(db=db, id=sid, retmode=retmode)
                    else:
                        handle = Entrez.esearch(db=db, term=sid)
                    record = Entrez.read(handle)
                    read_records.append(record[0])
                except:
                    failures.append("\t" + str(sid))

        end_time = time.time()
        duration = end_time - start_time
        hours, remainder = divmod(duration, 3600)
        minutes, seconds = divmod(remainder, 60)
        if duration < 0.8:
            time.sleep(0.5)
        durations.append(str(i) + ' - ' + str(i + chunk_size) + "\t" + ':'.join([str(minutes), str(round(seconds, 2))]))

    logging.debug("Entrez query time for accessions (minutes:seconds):\n\t" +
                  "\n\t".join(durations) + "\n")

    if failures:
        logging.warning("Unable to parse XML data from Entrez! "
                        "Either the XML is corrupted or the query terms cannot be found in the database.\n"
                        "Offending accessions from this batch:\n" + "\n".join(failures) + "\n")
    return read_records, durations, failures


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


def parse_gbseq_info_from_entrez_xml(record, gb_key="GBSeq_organism"):
    """
    Function for pulling out a value for a specific GenBank key from a dictionary
    :param record:
    :param gb_key:
    :return:
    """
    gb_value = ""
    if len(record) >= 1:
        try:
            gb_value = record[gb_key]
            # To prevent Entrez.efectch from getting confused by non-alphanumeric characters:
            gb_value = re.sub('[)(\[\]]', '', gb_value)
        except (IndexError, KeyError):
            logging.warning("'" + gb_key + "' not found in Entrez record.\n")
    return gb_value


def parse_gbseq_info_from_esearch_record(record, gb_key="IdList"):
    gb_value = ""
    if len(record) >= 1:
        try:
            gb_value = record[gb_key][0]
        except (IndexError, KeyError, TypeError):
            return gb_value
    return gb_value


def prep_for_entrez_query():
    """
    Tests checks to ensure the correct version of BioPython is imported,
    sends a test Entrez.efetch query to see if the internet connection is currently stable.
    :return:
    """
    if float(Bio.__version__) < 1.68:
        # This is required due to a bug in earlier versions returning a URLError
        logging.error("Version of biopython needs to be >=1.68! " +
                      str(Bio.__version__) + " is currently installed.\n")
        sys.exit(9)

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


def get_multiple_lineages(search_term_list: list, molecule_type: str):
    """
    Function for retrieving taxonomic lineage information from accession IDs - accomplished in 2 steps:
     1. Query Entrez's Taxonomy database using accession IDs to obtain corresponding organisms
     2. Query Entrez's Taxonomy database using organism names to obtain corresponding TaxIds
     3. Query Entrez's Taxonomy database using TaxIds to obtain corresponding taxonomic lineages
    :param search_term_list:
    :param molecule_type:
    :return:
    """
    if not search_term_list:
        logging.error("Search_term for Entrez query is empty\n")
        sys.exit(9)

    prep_for_entrez_query()
    entrez_db = validate_target_db(molecule_type)

    accession_lineage_map = dict()
    all_accessions = set()
    unique_organisms = set()
    unique_taxa = set()
    updated_accessions = dict()

    ##
    # Step 1: Query Entrez's Taxonomy database using accession IDs to obtain corresponding organisms
    ##
    logging.info("Retrieving taxonomy Entrez records for each accession... ")
    records_batch, durations, org_failures = tolerant_entrez_query(search_term_list, entrez_db)
    logging.info("done.\n")

    # Parse the records returned by tolerant_entrez_query mapping accessions to organism names
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
        tax_organism = parse_gbseq_info_from_entrez_xml(record)
        tax_lineage = check_lineage(parse_gbseq_info_from_entrez_xml(record, "GBSeq_taxonomy"), tax_organism)

        # Add the organism to unique_organisms set for taxonomic lineage querying
        accession_lineage_map[(accession, versioned)]["organism"] = tax_organism
        accession_lineage_map[(accession, versioned)]["tax_id"] = ""
        # If the full taxonomic lineage was not found, then add it to the unique organisms for further querying
        if len(tax_lineage) >= 7:
            accession_lineage_map[(accession, versioned)]["lineage"] = "; ".join(tax_lineage)
        else:
            accession_lineage_map[(accession, versioned)]["lineage"] = ""
            unique_organisms.add(tax_organism)
        all_accessions.update([accession, versioned])

    ##
    # Step 2: Query Entrez's Taxonomy database using organism names to obtain corresponding taxonomic lineages
    ##
    logging.info("Retrieving NCBI taxonomy IDs for each organism... ")
    records_batch, durations, lin_failures = tolerant_entrez_query(list(unique_organisms),
                                                                   "Taxonomy", "search", "xml", 1)
    logging.info("done.\n")
    for record in records_batch:
        try:
            organism = parse_gbseq_info_from_esearch_record(record, 'TranslationStack')['Term']
        except (IndexError, KeyError, TypeError):
            logging.warning("Value for 'TranslationStack' not found in Entrez record:" + str(record) + ".\n" +
                            "Unable to link taxonomy ID to organism.\n")
            continue
        tax_id = parse_gbseq_info_from_esearch_record(record)
        for acc_tuple in accession_lineage_map:
            if re.search(accession_lineage_map[acc_tuple]["organism"], organism):
                accession_lineage_map[acc_tuple]["tax_id"] = tax_id
                unique_taxa.add(tax_id)

    ##
    # Step 3: Fetch the taxonomic lineage for all of the taxonomic IDs
    ##
    logging.info("Retrieving lineage information for each taxonomy ID... ")
    records_batch, durations, lin_failures = tolerant_entrez_query(list(unique_taxa))
    logging.info("done.\n")
    for record in records_batch:
        tax_id = parse_gbseq_info_from_entrez_xml(record, "TaxId")
        for tuple_key in accession_lineage_map:
            if accession_lineage_map[tuple_key]["tax_id"] == tax_id:
                accession_lineage_map[tuple_key]["lineage"] = parse_gbseq_info_from_entrez_xml(record, "Lineage")

    # Report on the tolerance for failed Entrez accession queries
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
    logging.debug("\nQueries mapped ideally = " + str(success) +
                  "\nQueries mapped with alternative accessions = " + str(rescued) + "\n")

    return accession_lineage_map, all_accessions


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
        accession_lineage_map, all_accessions = get_multiple_lineages([ref_seq.accession for ref_seq in failed_accession_queries], molecule)
        for reference_sequence in failed_accession_queries:
            lineage = ""
            for tuple_key in accession_lineage_map.keys():
                if reference_sequence.accession in tuple_key:
                    lineage = accession_lineage_map[tuple_key]["lineage"]
                    break
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


if __name__ == "main":
    tolerant_entrez_query(['12968'])

__author__ = 'Connor Morgan-Lang'

import sys
import time
import re
import Bio
import logging
from Bio import Entrez
from urllib import error
from .utilities import clean_lineage_string


class ReferenceSequence:
    def __init__(self):
        self.accession = ""
        self.description = ""
        self.organism = ""
        self.lineage = ""
        self.short_id = ""
        self.sequence = ""
        self.locus = ""
        self.ncbi_tax = ""
        self.cluster_rep = False
        self.cluster_rep_similarity = 0
        self.cluster_lca = None

    def get_info(self):
        """
        Returns a string with the ReferenceSequence instance's current fields

        :return: str
        """
        info_string = ""
        info_string += "accession = " + self.accession + ", " + "treesapp_id = " + self.short_id + "\n"
        info_string += "organism = " + self.organism + ", " + "NCBI taxid = " + self.ncbi_tax + "\n"
        info_string += "description = " + self.description + ", " + "locus = " + self.locus + "\n"
        info_string += "lineage = " + self.lineage + "\n"
        return info_string


class EntrezRecord(ReferenceSequence):
    def __init__(self, acc, ver):
        super().__init__()
        self.accession = acc
        self.versioned = ver
        self.bitflag = 0  # For monitoring progress during download stage

    def get_info(self):
        info_string = super().get_info()
        info_string += "acc.version = " + self.versioned + ", bitflag = " + str(self.bitflag)
        return info_string

    def tracking_stamp(self):
        self.bitflag = 0
        if self.organism:
            self.bitflag += 1
        if self.ncbi_tax:
            self.bitflag += 2
        if self.lineage:
            self.bitflag += 4
        return


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
                      "Please create an issue on the GitHub page.\n")
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

    if len(search_term_list) == 0:
        return read_records, durations, failures

    # Create the sub-lists from `search_term_list` of length `chunk_size`
    for i in range(0, len(search_term_list), chunk_size):
        start_time = time.time()
        sub_list = search_term_list[i:i + chunk_size]
        try:
            if method == "fetch":
                handle = Entrez.efetch(db=db, id=','.join([str(sid) for sid in sub_list]), retmode=retmode,
                                       api_key="849e32266531ee0cee64c6edbbdcf7b62e09")
            else:
                handle = Entrez.esearch(db=db, term=','.join([str(sid) for sid in sub_list]),
                                        api_key="849e32266531ee0cee64c6edbbdcf7b62e09")
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
                        handle = Entrez.efetch(db=db, id=sid, retmode=retmode,
                                               api_key="849e32266531ee0cee64c6edbbdcf7b62e09")
                    else:
                        handle = Entrez.esearch(db=db, term=sid,
                                                api_key="849e32266531ee0cee64c6edbbdcf7b62e09")
                    record = Entrez.read(handle)
                    read_records.append(record[0])
                except:
                    failures.append("\t" + str(sid))

        end_time = time.time()
        duration = end_time - start_time
        if duration < 0.66:
            time.sleep(0.66 - duration)
            duration = 0.66
        # TODO: allow for the durations to be summed by a number of queries
        hours, remainder = divmod(duration, 3600)
        minutes, seconds = divmod(remainder, 60)
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
    accession_keys = ["GBSeq_primary-accession"]  # NOTE: "GBSeq_locus" found to occassionally map to random accessions
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
                    alternatives.append(re.search(r"\|+(.*)$", alt).group(1))
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


def check_lineage(lineage: str, organism_name: str, verbosity=0):
    """
    Sometimes the NCBI lineage is incomplete.
    Currently, this function uses organism_name to ideally add Species to the lineage
    :param lineage: A semi-colon separated taxonomic lineage
    :param organism_name: Name of the organism. Parsed from the sequence header (usually at the end in square brackets)
    :param verbosity: 1 prints debugging messages
    :return: A list of elements for each taxonomic rank representing the taxonomic lineage
    """
    if verbosity:
        logging.debug("check_lineage():\n\tlineage = '" + lineage + "'\n\torganism = '" + organism_name + "'\n")

    if not lineage:
        return []
    proper_species_re = re.compile("^[A-Z][a-z]+ [a-z]+$")
    lineage_list = clean_lineage_string(lineage).split("; ")
    if proper_species_re.match(lineage_list[-1]):
        if verbosity:
            logging.debug("check_lineage(): Perfect lineage.\n")
    elif len(lineage_list) >= 6 and proper_species_re.match(organism_name):
        if verbosity:
            logging.debug("check_lineage(): Organism name added to complete the lineage.\n")
        lineage_list.append(organism_name)
    elif len(lineage_list) < 6 and organism_name != lineage_list[-1] and re.match("^[A-Z][a-z]+$", organism_name):
        if verbosity:
            logging.debug("check_lineage(): Organism name added to truncated lineage.\n")
        lineage_list.append(organism_name)
    else:
        if verbosity:
            logging.debug("check_lineage(): Bad lineage.\n")
    return lineage_list


def match_file_to_dict(file_handler, key_dict, sep="\t", join_by=0):
    """
    Generator function for mapping a particular field in a file, separated by 'sep', to dictionary keys
    :param file_handler: Opened file object
    :param key_dict: Dictionary to map the selected field to
    :param sep: Field separator
    :param join_by: The field number to search the dictionary for
    :return: Line that matches a key
    """
    for line in file_handler:
        if line.split(sep)[join_by] in key_dict:
            yield line


def map_accession2taxid(query_accession_list, accession2taxid_list):
    er_acc_dict = dict()
    unmapped_queries = list()

    # Create a dictionary for O(1) look-ups, load all the query accessions into unmapped queries
    for e_record in query_accession_list:  # type: EntrezRecord
        try:
            er_acc_dict[e_record.accession].append(e_record)
        except KeyError:
            er_acc_dict[e_record.accession] = [e_record]
        unmapped_queries.append(e_record.accession)

    logging.info("Mapping query accessions to NCBI taxonomy IDs... ")
    for accession2taxid in accession2taxid_list.split(','):
        init_qlen = len(unmapped_queries)
        final_qlen = len(unmapped_queries)
        start = time.time()
        try:
            rosetta_handler = open(accession2taxid, 'r')
        except IOError:
            logging.error("Unable to open '" + accession2taxid + "' for reading.\n")
            sys.exit(13)

        for line_match in match_file_to_dict(rosetta_handler, er_acc_dict):
            try:
                accession, ver, taxid, _ = line_match.strip().split("\t")
            except (ValueError, IndexError):
                logging.warning("Parsing '" + accession2taxid + "' failed.\n")
                break

            try:
                # Update the EntrezRecord elements
                records = er_acc_dict[accession]
                for record in records:
                    if not record.versioned:
                        record.versioned = ver
                    record.ncbi_tax = taxid
                    record.bitflag = 3  # Necessary for downstream filters - indicates taxid has been found
                # Remove accession from unmapped queries
                i = 0
                while i < final_qlen:
                    if unmapped_queries[i] == accession:
                        unmapped_queries.pop(i)
                        final_qlen -= 1
                        break
                    i += 1
                if final_qlen == 0:
                    break
            except KeyError:
                logging.error("Bad key returned by generator.\n")
                sys.exit(13)

        rosetta_handler.close()
        end = time.time()
        logging.debug("Time required to parse '" + accession2taxid + "': " + str(round(end - start, 1)) + "s.\n")
        # Report the number percentage of query accessions mapped
        logging.debug(
            str(round(((init_qlen - final_qlen) * 100 / len(query_accession_list)), 2)) +
            "% of query accessions mapped by " + accession2taxid + ".\n")
    logging.info("done.\n")

    return er_acc_dict


def pull_unmapped_entrez_records(entrez_records: list):
    """
    Prepares a list of accession identifiers for EntrezRecord instances where the bitflag is not equal to 7,
     inferring lineage information was not properly entered.
    :param entrez_records: A list of EntrezRecord instances
    :return: List of EntrezRecords with bitflag != 7
    """
    unmapped_queries = list()
    x = 0
    while x < len(entrez_records):
        e_record = entrez_records[x]  # type: EntrezRecord
        e_record.tracking_stamp()
        if e_record.bitflag < 6:  # Set threshold to 6 as the organism may not have been entered
            unmapped_queries.append(e_record)
            entrez_records.pop(x)
        else:
            x += 1
    return unmapped_queries


def fetch_lineages_from_taxids(entrez_records: list):
    """
    Query Entrez's Taxonomy database for lineages using NCBI taxonomic IDs.
    The TaxId queries are pulled from EntrezRecord instances.
    :param entrez_records: A list of EntrezRecord instances that should have TaxIds in their ncbi_tax element
    :return: entrez_records where successful queries have a populated lineage element
    """
    tax_id_map = dict()

    prep_for_entrez_query()

    # Create a dictionary that will enable rapid look-ups and mapping to EntrezRecord instances
    for e_record in entrez_records:  # type: EntrezRecord
        if e_record.bitflag >= 4:
            # Lineage has been added
            continue
        elif e_record.bitflag < 2:
            # NCBI taxonomy ID has not been added
            logging.debug("Empty NCBI taxonomy ID for incomplete EntrezRecord query:\n" + e_record.get_info() + "\n")
            continue
        taxid = e_record.ncbi_tax
        if taxid and taxid not in tax_id_map:
            tax_id_map[taxid] = []
        tax_id_map[taxid].append(e_record)

    logging.info("Retrieving lineage information for each taxonomy ID... ")
    records_batch, durations, lin_failures = tolerant_entrez_query(list(tax_id_map.keys()))
    logging.info("done.\n")
    for record in records_batch:
        tax_id = parse_gbseq_info_from_entrez_xml(record, "TaxId")
        if len(tax_id) == 0:
            logging.warning("Empty TaxId returned in Entrez XML.\n")
        tax_lineage = parse_gbseq_info_from_entrez_xml(record, "Lineage")
        try:
            for e_record in tax_id_map[tax_id]:
                e_record.lineage = tax_lineage
                e_record.tracking_stamp()
        except KeyError:
            pass
    return entrez_records


def entrez_records_to_accession_set(entrez_records_list: list, bitflag_filter=7):
    query_dict = dict()
    for record in entrez_records_list:  # type: EntrezRecord
        if record.bitflag > bitflag_filter:
            continue
        # Uses the versioned accession ID if available, regular accession otherwise
        if record.versioned:
            query_dict[record.versioned] = record
        elif record.accession:
            query_dict[record.accession] = record
        else:
            continue
    return query_dict


def entrez_records_to_organism_set(entrez_records_list: list, bitflag_filter=7):
    query_dict = dict()
    for record in entrez_records_list:  # type: EntrezRecord
        if record.bitflag > bitflag_filter:

            continue
        if record.organism:
            record.organism = re.sub('[:]', ' ', record.organism)
            record.organism = re.sub(' =.*', '', record.organism)
            record.organism = re.sub(r' \(.*\)', '', record.organism)
            query = record.organism + "[All Names]"
            # Index the accession_lineage_map by organism and map to list of EntrezRecord object
            if query in query_dict:
                query_dict[query].append(record)
            else:
                query_dict[query] = [record]
        else:
            continue
    return query_dict


def entrez_records_to_accession_lineage_map(entrez_records_list):
    # TODO: Remove this necessity. Currently need to reformat and tally accessions like so but its a waste
    # Used for tallying the status of Entrez queries
    success = 0
    rescued = 0
    bad_tax = 0
    bad_org = 0
    failed = 0
    accession_lineage_map = dict()

    for e_record in entrez_records_list:  # type: EntrezRecord
        e_record.tracking_stamp()
        if e_record.bitflag == 0:
            failed += 1
            continue
        # Report on the tolerance for failed Entrez accession queries
        if e_record.bitflag == 7:
            success += 1
        elif e_record.bitflag == 6:
            rescued += 1
        elif 5 >= e_record.bitflag >= 2:
            bad_tax += 1
        elif e_record.bitflag == 1:
            bad_org += 1
        else:
            logging.error("Unexpected bitflag (" + str(e_record.bitflag) + ") encountered for EntrezRecord:\n" +
                          e_record.get_info() + "tax_id = " + str(e_record.ncbi_tax) + "\n")
            sys.exit(19)
        e_record_key = (e_record.accession, e_record.versioned)
        if e_record_key in accession_lineage_map and e_record.lineage != accession_lineage_map[e_record_key]["lineage"]:
            logging.warning(str(e_record_key) + " already present in accession-lineage map with different lineage.\n" +
                            "The most complete lineage will be used.\n")
            if len(e_record.lineage) <= len(accession_lineage_map[e_record_key]["lineage"]):
                continue
        accession_lineage_map[e_record_key] = dict()
        accession_lineage_map[e_record_key]["lineage"] = e_record.lineage
        accession_lineage_map[e_record_key]["organism"] = e_record.organism

    logging.debug("Queries mapped ideally = " + str(success) +
                  "\nQueries with organism unmapped = " + str(bad_org) +
                  "\nQueries with NCBI taxonomy ID unmapped = " + str(bad_tax) +
                  "\nQueries mapped with alternative accessions = " + str(rescued) +
                  "\nQueries that outright failed = " + str(failed) + "\n")

    return accession_lineage_map


def get_multiple_lineages(entrez_query_list: list, molecule_type: str):
    """
    Function for retrieving taxonomic lineage information from accession IDs - accomplished in 2 steps:
     1. Query Entrez's Taxonomy database using accession IDs to obtain corresponding organisms
     2. Query Entrez's Taxonomy database using organism names to obtain corresponding TaxIds
     3. Query Entrez's Taxonomy database using TaxIds to obtain corresponding taxonomic lineages

    :param entrez_query_list: A list of GenBank accession IDs to be mapped to lineages
    :param molecule_type: The type of molecule (e.g. prot, nuc) to be mapped to a proper Entrez database name
    :return: List of EntrezRecord instances
    """
    if not entrez_query_list:
        logging.error("Search_term for Entrez query is empty\n")
        sys.exit(9)

    prep_for_entrez_query()
    entrez_db = validate_target_db(molecule_type)

    ##
    # Step 1: Query Entrez's Taxonomy database using accession IDs to obtain corresponding organisms
    ##
    search_terms = entrez_records_to_accession_set(entrez_query_list, 1)
    logging.info("Retrieving Entrez taxonomy records for each accession... ")
    records_batch, durations, org_failures = tolerant_entrez_query(list(search_terms.keys()), entrez_db)
    logging.info("done.\n")

    # Parse the records returned by tolerant_entrez_query, mapping accessions to organism names
    for record in records_batch:
        e_record = None
        accession, ver, alt = parse_accessions_from_entrez_xml(record)
        if accession in search_terms:
            e_record = search_terms[accession]
        elif ver in search_terms:
            e_record = search_terms[ver]
        else:
            for alt_key in alt:
                if alt_key in search_terms:
                    e_record = search_terms[alt_key]
                    break
        if not e_record:
            logging.warning("Unable to map neither a record' accession or accession.version to an EntrezRecord:\n" +
                            "Accession: '" + str(accession) + "'\n" +
                            "Acc.Version: '" + str(ver) + "'\n" +
                            "Alternatives:" + str(alt) + "\n" +
                            str(record) + "\n")
            continue
        e_record.organism = parse_gbseq_info_from_entrez_xml(record)
        # Entrez replaces special characters with whitespace in organism queries, so doing it here for compatibility
        tax_lineage = check_lineage(parse_gbseq_info_from_entrez_xml(record, "GBSeq_taxonomy"), e_record.organism)

        # If the full taxonomic lineage was not found, then add it to the unique organisms for further querying
        if len(tax_lineage) >= 7 or tax_lineage[-1] == e_record.organism:
            e_record.lineage = "; ".join(tax_lineage)
        e_record.tracking_stamp()

    ##
    # Step 2: Query Entrez's Taxonomy database using organism names to obtain corresponding taxonomic lineages
    ##
    search_terms = entrez_records_to_organism_set(entrez_query_list, 3)
    logging.debug(str(len(search_terms.keys())) + " unique organism queries.\n")
    logging.info("Retrieving NCBI taxonomy IDs for each organism... ")
    records_batch, durations, taxid_failures = tolerant_entrez_query(list(search_terms.keys()),
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
        if not tax_id:
            logging.warning("Entrez returned an empty TaxId for organism '" + organism + "'\n")
        try:
            # This can, and will, lead to multiple accessions being assigned the same tax_id - not a problem, though
            for e_record in search_terms[organism]:
                if e_record.bitflag == 7:
                    continue
                e_record.ncbi_tax = tax_id
                e_record.tracking_stamp()
        except KeyError:
            logging.warning("Unable to map organism '" + organism + "' to an EntrezRecord:\n")
            continue

    entrez_query_list = fetch_lineages_from_taxids(entrez_query_list)

    return entrez_query_list


def verify_lineage_information(accession_lineage_map, fasta_record_objects, taxa_searched):
    """
    Function used for parsing records returned by Bio.Entrez.efetch queries and identifying inconsistencies
    between the search terms and the results

    :param accession_lineage_map: A dictionary mapping accession.versionID tuples to taxonomic lineages
    :param fasta_record_objects:
    :param taxa_searched: An integer for tracking number of accessions queried (currently number of lineages provided)
    :return:
    """
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
    for treesapp_id in sorted(fasta_record_objects.keys()):
        ref_seq = fasta_record_objects[treesapp_id]  # type: EntrezRecord
        ref_seq.tracking_stamp()
        # Could have been set previously, in custom header format for example
        if ref_seq.lineage:
            unambiguous_accession_lineage_map[ref_seq.accession] = clean_lineage_string(ref_seq.lineage)
        else:
            lineage = ""
            for tuple_key in accession_lineage_map:
                accession, versioned = tuple_key
                if ref_seq.accession == accession or ref_seq.accession == versioned:
                    if accession_lineage_map[tuple_key]["lineage"] == "":
                        lineage = "Unclassified"
                    else:
                        # The query was successful! Add it and increment
                        lineage = accession_lineage_map[tuple_key]["lineage"]
                    if not ref_seq.organism and accession_lineage_map[tuple_key]["organism"]:
                        ref_seq.organism = accession_lineage_map[tuple_key]["organism"]

            if not lineage and ref_seq.bitflag == 0:
                logging.error("Lineage information was not retrieved for " + ref_seq.accession + "!\n" +
                              "Please remove the output directory and restart.\n")
                sys.exit(13)
            elif not lineage and ref_seq.bitflag >= 1:
                lineage = "Unclassified"

            ref_seq.lineage = clean_lineage_string("; ".join(check_lineage(lineage, ref_seq.organism)))
            unambiguous_accession_lineage_map[ref_seq.accession] = ref_seq.lineage
        ref_seq.tracking_stamp()
        if ref_seq.bitflag >= 1:
            taxa_searched += 1

    if taxa_searched < len(fasta_record_objects.keys()):
        logging.error("Not all sequences (" + str(taxa_searched) + '/'
                      + str(len(fasta_record_objects)) + ") were queried against the NCBI taxonomy database!\n")
        sys.exit(9)

    return fasta_record_objects, unambiguous_accession_lineage_map


def read_accession_taxa_map(mapping_file):
    """
    A function for reading intermediate files made by write_accession_lineage_map to avoid the time-consuming download

    :param mapping_file:
    :return: accession_lineage_map:  A dictionary mapping accessions to lineages
    """
    try:
        map_file_handler = open(mapping_file, 'r')
    except (IOError, FileNotFoundError):
        logging.error("Unable to open " + mapping_file + " for reading!\n")
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
    Function to create data collections to fulfill entrez query searches.
    The queries are EntrezRecord instances lacking lineage information; they could have either accessions or NCBI taxids

    :param fasta_record_objects: A list of ReferenceSequence objects - lineage information to be filled
    :return: List containing Entrez queries
    """
    num_lineages_provided = 0
    entrez_query_list = list()
    unavailable = list()
    for num_id in fasta_record_objects:
        ref_seq = fasta_record_objects[num_id]  # type: EntrezRecord
        # Only need to download the lineage information for those sequences that don't have it encoded in their header
        if ref_seq.lineage:
            num_lineages_provided += 1
        elif ref_seq.accession or ref_seq.ncbi_tax:
            entrez_query_list.append(ref_seq)
        else:
            unavailable.append(ref_seq.get_info())
    if len(unavailable) > 0:
        logging.warning("Neither accession nor lineage available for:\n\t" +
                        "\n\t".join(unavailable))
    return list(entrez_query_list), num_lineages_provided


def fill_ref_seq_lineages(fasta_record_objects, accession_lineages):
    """
    Adds lineage information from accession_lineages to fasta_record_objects

    :param fasta_record_objects: dict() indexed by TreeSAPP numeric identifiers mapped to ReferenceSequence instances
    :param accession_lineages: a dictionary mapping {accession: lineage}
    :return: None
    """
    for treesapp_id in fasta_record_objects:
        ref_seq = fasta_record_objects[treesapp_id]  # type: EntrezRecord
        if not ref_seq.lineage:
            try:
                lineage = accession_lineages[ref_seq.accession]
            except KeyError:
                logging.error("Lineage information was not retrieved for accession '" + ref_seq.accession + "'.\n" +
                              "Please remove the output directory and restart.\n")
                sys.exit(13)
            # Add the species designation since it is often not included in the sequence record's lineage
            ref_seq.lineage = lineage
        if not ref_seq.organism and ref_seq.lineage:
            ref_seq.organism = ref_seq.lineage.split("; ")[-1]
        else:
            pass
        ref_seq.tracking_stamp()
    return


def load_ref_seqs(fasta_dict, header_registry, ref_seq_dict):
    """
    Function for adding sequences from a fasta-formatted dictionary into dictionary of ReferenceSequence objects

    :param fasta_dict:
    :param header_registry: An optional dictionary of Header objects
    :param ref_seq_dict: A dictionary indexed by arbitrary integers mapping to ReferenceSequence instances
    :return: None
    """
    missing = list()
    if len(header_registry) != len(fasta_dict):
        logging.warning("Number of records in FASTA collection and header list differ.\n" +
                        "Chances are these were short sequences that didn't pass the filter. Carrying on.\n")

    for num_id in sorted(ref_seq_dict.keys(), key=int):
        ref_seq = ref_seq_dict[num_id]
        formatted_header = header_registry[num_id].formatted
        try:
            ref_seq.sequence = fasta_dict[formatted_header]
        except KeyError:
            if len(header_registry) == len(fasta_dict):
                logging.error(formatted_header + " not found in FASTA records due to format incompatibilities.\n")
                sys.exit(21)
            missing.append(str(header_registry[num_id].original))
    if len(missing) > 0:
        logging.debug("The following sequences have been removed from further analyses:\n\t" +
                      "\n\t".join(missing) + "\n")
    return


def map_accessions_to_lineages(query_accession_list: list, molecule: str, accession_to_taxid=None):
    if accession_to_taxid:
        # Determine find the query accessions that are located in the provided accession2taxid file
        entrez_record_dict = map_accession2taxid(query_accession_list, accession_to_taxid)
        entrez_records = []
        for index in entrez_record_dict:
            entrez_records += entrez_record_dict[index]
        # Map lineages to taxids for successfully-mapped query sequences
        fetch_lineages_from_taxids(entrez_records)
        # Use the normal querying functions to obtain lineage information for the unmapped queries
        unmapped_queries = pull_unmapped_entrez_records(entrez_records)
        if len(unmapped_queries) > 0:
            # This tends to be a minority so shouldn't be too taxing
            for e_record in get_multiple_lineages(unmapped_queries, molecule):  # type: EntrezRecord
                try:
                    entrez_record_dict[e_record.accession].append(e_record)
                except KeyError:
                    logging.warning(e_record.accession + " not found in original query list.\n")
                    continue
                entrez_records.append(e_record)
        entrez_record_dict.clear()
        unmapped_queries.clear()
    else:
        entrez_records = get_multiple_lineages(query_accession_list, molecule)
    return entrez_records


if __name__ == "main":
    tolerant_entrez_query(['12968'])

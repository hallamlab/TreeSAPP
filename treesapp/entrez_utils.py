__author__ = 'Connor Morgan-Lang'

import sys
import time
import re
import logging
import csv

from Bio import Entrez
from urllib import error

from treesapp.utilities import get_list_positions, get_field_delimiter
from treesapp.taxonomic_hierarchy import TaxonomicHierarchy, Taxon


class EntrezRecord:
    def __init__(self, acc: str, ver: str):
        """
        Instantiation function for the EntrezRecord class.
        This class is used for querying the Entrez database, typically to retrieve taxonomic lineage information for
        a given sequence deposited in GenBank, or just to get the lineage from a NCBI taxid.

        The remaining class attributes are populated while parsing records returned by various queries.

        :param acc: The accession of the Entrez record
        :param ver: The versioned accession of the Entrez record
        """
        self.accession = acc
        self.versioned = ver
        self.short_id = ""  # A unique alphanumerical TreeSAPP ID
        self.ncbi_tax = ""  # NCBI's taxonomy ID
        self.organism = ""
        self.lineage = ""
        self.taxon_rank = ""  # Taxonomic rank the organism was described to
        self.sequence = ""  # Nucleotide or amino acid sequence
        self.locus = ""
        self.description = ""
        self.cluster_rep = True
        self.cluster_rep_similarity = 0
        self.cluster_lca = None
        self.bitflag = 0  # For monitoring progress during download stage

    def get_info(self) -> str:
        """
        Returns a string with the EntrezRecord instance's current variables

        :return: str
        """
        info_string = "Information for EntrezRecord ID '" + str(self.short_id) + "':\n"
        info_string += "accession = " + self.accession + ", " + "acc.version = " + self.versioned + "\n"
        info_string += "organism = " + str(self.organism) + ", " + "rank resolved = " + self.taxon_rank + "\n"
        info_string += "NCBI taxid = " + str(self.ncbi_tax) + ", " + "bitflag = " + str(self.bitflag) + "\n"
        info_string += "lineage = " + str(self.lineage) + "\n"
        info_string += "description = " + str(self.description) + ", " + "locus = " + str(self.locus) + "\n"
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

    def rebuild_header(self):
        return ' '.join([self.versioned, self.description]).rstrip()


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


def parse_gbseq_info_from_entrez_xml(record: dict, gb_key="GBSeq_organism"):
    """
    Function for pulling out a value for a specific GenBank key from a dictionary

    :param record: A Entrez.Parser.DictionaryElement containing the gb_key
    :param gb_key: A string that refers to a key in a key: value pair in record
    :return: Either a StringElement or ListElement
    """
    gb_value = []
    if len(record) >= 1:
        try:
            gb_value = record[gb_key]
            # # To prevent Entrez.efetch from getting confused by non-alphanumeric characters:
            # try:
            #     gb_value = re.sub(r'[)(\[\]]', '', gb_value)
            # except TypeError:
            #     return gb_value
        except (IndexError, KeyError):
            logging.debug("'" + gb_key + "' not found in Entrez record:\n" + str(record) + "\n")
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

    :return: None
    """

    logging.info("Preparing Bio.Entrez for NCBI queries... ")
    Entrez.email = "c.morganlang@gmail.com"
    Entrez.tool = "treesapp"
    # Test the internet connection:
    try:
        record = Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
    except error.URLError:
        logging.warning("Unable to serve Entrez query. Are you connected to the internet?\n")
        record = None
    logging.info("done.\n")
    return record


def repair_conflict_lineages(t_hierarchy: TaxonomicHierarchy, ref_seq_dict: dict) -> None:
    """
    When some taxon nodes are removed from a TaxonomicHierarchy, the original lineages with those included may persist.

    The purpose of repair_conflict_lineages is to substitute the former nodes causing conflict for the replacement nodes
    in the taxonomic lineages of

    :return: None
    """
    if len(t_hierarchy.conflicts) == 0:
        return

    nodes_replaced_map = t_hierarchy.resolve_conflicts()  # return taxa whose nodes were merged

    for old_taxon, new_taxon in nodes_replaced_map.items():  # type: (Taxon, Taxon)
        taxon_name = old_taxon.prefix_taxon()
        for record in ref_seq_lineage_scanner(ref_seq_dict, taxon_name):  # type: EntrezRecord
            # Find the name of the organism to use for building the new taxonomic lineage
            if t_hierarchy.canonical_prefix.match(record.organism):
                organism_query = record.organism
            elif t_hierarchy.canonical_prefix.match(record.lineage.split(t_hierarchy.lin_sep)[-1]):
                organism_query = record.lineage.split(t_hierarchy.lin_sep)[-1]
            else:
                continue

            # If the current record's most resolved taxon has been substituted
            # swap all taxonomic information to the new representative
            if organism_query == taxon_name:
                record.organism = new_taxon.prefix_taxon()
                organism_query = record.organism

            ref_taxon = t_hierarchy.get_taxon(organism_query)  # type: Taxon
            try:
                record.lineage = t_hierarchy.lin_sep.join([taxon.prefix_taxon() for taxon in ref_taxon.lineage()])
                continue
            except AttributeError:
                logging.warning("Unable to repair the conflicted lineage of record {}, "
                                "'{}'".format(record.accession, record.lineage))
                continue
    return


def ref_seq_lineage_scanner(ref_seq_dict: dict, taxon_name: str) -> list:
    """
    Finds the EntrezRecord instances in ref_seq_dict with lineages that contain a specific taxon (taxon_name).

    :param ref_seq_dict: A dictionary containing unique, numerical TreeSAPP identifiers mapper to EntrezRecord objects
    :param taxon_name: The string to be used for matching EntrezRecord.lineage attributes.
     It should represent a single taxon.
    :return: A list of EntrezRecord instances whose lineage attributes contain taxon_name
    """
    ref_seq_matches = []
    for treesapp_id, ref_seq in ref_seq_dict.items():  # type: (str, EntrezRecord)
        if ref_seq.lineage and taxon_name in ref_seq.lineage:
            ref_seq_matches.append(ref_seq)
    return ref_seq_matches


def sync_record_and_hierarchy_lineages(ref_leaf_nodes: list, records: dict) -> None:
    for leaf_node in ref_leaf_nodes:
        records[leaf_node.number].lineage = leaf_node.lineage  # type: EntrezRecord
    return


def repair_lineages(ref_seq_dict: dict, t_hierarchy: TaxonomicHierarchy) -> None:
    """
    This is used for adding rank prefixes (e.g. d__) to the taxonomic lineages in a dictionary of EntrezRecord instances

    :param ref_seq_dict: A dictionary containing unique, numerical TreeSAPP identifiers mapper to EntrezRecord objects
    :param t_hierarchy: A TaxonomicHierarchy instance
    :return: None
    """
    to_repair = set()  # A set of treesapp_ids that need to be fixed
    unprefixed_lineages = set()
    tmp_lineages = set()

    # Search for any taxon that doesn't have a rank prefix in all the reference sequence lineages
    t_hierarchy.clean_trie = True
    repair_conflict_lineages(t_hierarchy, ref_seq_dict)
    t_hierarchy.build_multifurcating_trie(key_prefix=True)
    # TODO: handle the lineages with rank-prefixes but are absent from the hierarchy
    for treesapp_id in sorted(ref_seq_dict.keys()):  # type: str
        ref_seq = ref_seq_dict[treesapp_id]  # type: EntrezRecord
        if ref_seq.lineage:
            if t_hierarchy.clean_lineage_string(ref_seq.lineage) not in t_hierarchy.trie:
                to_repair.add(treesapp_id)
                unprefixed_lineages.add(ref_seq.lineage)  # It only takes one rank without a prefix to add it
        else:
            ref_seq.lineage = "r__Root"

    prep_for_entrez_query()
    # Build list of entrez queries for EntrezRecords with un-annotated lineages
    entrez_query_list = entrez_records_from_lineages_and_chop(unprefixed_lineages, tmp_lineages,
                                                              t_hierarchy.get_taxon_names())

    while entrez_query_list:
        logging.info("Repairing {0} taxonomic lineages for {1} references.\n".format(len(entrez_query_list),
                                                                                     len(to_repair)))
        # Gather NCBI taxid
        o_search_terms = entrez_records_to_organism_set(entrez_query_list, 3)
        fetch_taxids_from_organisms(o_search_terms)
        # Fetch the lineages for each NCBI taxid
        fetch_lineages_from_taxids(entrez_query_list, t_hierarchy)

        # Remove lineages from unprefixed_lineages if all ranks are repaired
        while tmp_lineages:
            lineage = tmp_lineages.pop()
            if not t_hierarchy.project_lineage(lineage_str=lineage):
                unprefixed_lineages.add(lineage)
        # Update entrez_query_list with remaining lineages in unprefixed_lineages
        entrez_query_list = entrez_records_from_lineages_and_chop(unprefixed_lineages, tmp_lineages,
                                                                  t_hierarchy.get_taxon_names())

    # Add rank prefixes to the broken lineages
    while to_repair:
        ref_seq = ref_seq_dict[to_repair.pop()]  # type: EntrezRecord
        ref_lineage = t_hierarchy.get_prefixed_lineage_from_bare(ref_seq.lineage)
        if not ref_lineage:
            t_hierarchy.so_long_and_thanks_for_all_the_fish("Unable to find any ranks from lineage '{}'"
                                                            " in taxonomic hierarchy.\n".format(ref_seq.lineage) +
                                                            "Consider providing these data in a table or"
                                                            " removing this sequence from your analysis.\n")
        ref_seq.lineage = ref_lineage

        if len(to_repair) == 0:
            logging.info("done.\n")

    t_hierarchy.root_domains(root=t_hierarchy.find_root_taxon())
    for treesapp_id in sorted(ref_seq_dict.keys()):  # type: str
        e_record = ref_seq_dict[treesapp_id]  # type: EntrezRecord
        e_record.lineage = t_hierarchy.check_lineage(e_record.lineage, e_record.organism)

    return


def fill_entrez_record_taxon_rank(entrez_record_map: dict, t_hierarchy: TaxonomicHierarchy) -> None:
    for treesapp_id, ref_seq in entrez_record_map.items():  # type: (str, EntrezRecord)
        if not ref_seq.taxon_rank:
            ref_seq.taxon_rank = t_hierarchy.resolved_to(ref_seq.lineage)
    return


def fill_ref_seq_lineages(entrez_record_map: dict, accession_lineages: dict, complete=True) -> None:
    """
    Adds lineage information from a dictionary (accession_lineages) to all of the EntrezRecord values in a dictionary
    indexed by their respective numerical TreeSAPP identifiers (entrez_record_map).

    :param entrez_record_map: A dictionary indexed by TreeSAPP numeric identifiers mapped to EntrezRecord instances
    :param accession_lineages: A dictionary mapping {accession: lineage}. acc.version format also accepted
    :param complete: Boolean indicating whether the accession_lineage contains lineages for all records or not.
    If True, iteration will continue when accession_lineages is missing an accession.
    :return: None
    """
    lineage_added = 0
    for treesapp_id in entrez_record_map:
        ref_seq = entrez_record_map[treesapp_id]  # type: EntrezRecord
        if not ref_seq.lineage:
            try:
                lineage = accession_lineages[ref_seq.accession]
            except KeyError:
                if ref_seq.versioned in accession_lineages:
                    lineage = accession_lineages[ref_seq.versioned]
                elif not complete:
                    continue
                else:
                    logging.error("Lineage information not retrieved for, or could not be mapped to, accession '{}'.\n"
                                  "Please remove the output directory and restart.\n".format(ref_seq.accession))
                    sys.exit(13)
            # Add the species designation since it is often not included in the sequence record's lineage
            ref_seq.lineage = lineage
            lineage_added += 1
        if not ref_seq.organism and ref_seq.lineage:
            ref_seq.organism = ref_seq.lineage.split("; ")[-1]
        else:
            pass
        # TODO: Come up with a better way of removing the organism name after filling the organism attribute
        if len(ref_seq.lineage.split("; ")) > 8 and not re.match(r"[a-z]__.*", ref_seq.organism):
            ref_seq.lineage = "; ".join(ref_seq.lineage.split("; ")[:-1])
        ref_seq.tracking_stamp()

    if lineage_added == 0:
        logging.debug("No lineages from the accession map were added to the EntrezRecord attributes.\n")

    return


def entrez_record_snapshot(entrez_records: dict) -> dict:
    er_snaps = dict()
    for index in entrez_records:
        ref_seq = entrez_records[index]  # type: EntrezRecord
        if ref_seq.cluster_rep:
            er_snaps[id(ref_seq)] = ref_seq
    return er_snaps


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


def map_accession2taxid(query_accessions: list, accession2taxid_list: str) -> dict:
    """
    Maps NCBI accessions to taxonomy IDs via NCBI .accession2taxid files

    :param query_accessions: A list of EntrezRecord instances with accessions that need to be mapped to NCBI taxids
    :param accession2taxid_list: A comma-separated list of files
    :return: A dictionary mapping sequence accessions to EntrezRecord instances
    """
    er_acc_dict = dict()
    unmapped_queries = list()

    # Create a dictionary for O(1) look-ups, load all the query accessions into unmapped queries
    for e_record in query_accessions:  # type: EntrezRecord
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
            str(round(((init_qlen - final_qlen) * 100 / len(query_accessions)), 2)) +
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


def fetch_lineages_from_taxids(entrez_records: list, t_hierarchy=None) -> None:
    """
    Query Entrez's Taxonomy database for lineages using NCBI taxonomic IDs.
    The TaxId queries are pulled from EntrezRecord instances.
    The lineage from each successful Entrez query is fed into the TaxonomicHierarchy and can be subsequently queried

    :param entrez_records: A list of EntrezRecord instances that should have TaxIds in their ncbi_tax element
    :param t_hierarchy: A TaxonomicHierarchy instance
    :return: None
    """
    tax_id_map = dict()
    pulled_tax_ids = set()
    if not t_hierarchy:
        t_hierarchy = TaxonomicHierarchy()
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
            continue
        pulled_tax_ids.add(tax_id)
        tax_lineage = parse_gbseq_info_from_entrez_xml(record, "Lineage")
        tax_organism = parse_gbseq_info_from_entrez_xml(record, "ScientificName")
        tax_rank = parse_gbseq_info_from_entrez_xml(record, "Rank")
        lineage_ex = parse_gbseq_info_from_entrez_xml(record, "LineageEx")
        if not lineage_ex:
            logging.debug("Unable to find taxonomic ranks for organism '{0}' in record:\n"
                          "{1}\n.".format(tax_organism, record))
            continue

        # If the organism name isn't the last element of the lineage, add it as well as its rank
        if tax_lineage.split("; ")[-1] != tax_organism:
            tax_lineage += "; " + tax_organism
            lineage_ex += [{"ScientificName": tax_organism, "Rank": tax_rank}]

        taxon = t_hierarchy.feed(tax_lineage, lineage_ex)  # type: Taxon
        # We don't want to begin accumulating coverage at this stage
        for t in taxon.lineage():
            t.coverage = 0
        lineage_anno = t_hierarchy.emit(taxon.prefix_taxon(), True)

        try:
            for e_record in tax_id_map[tax_id]:  # type: EntrezRecord
                e_record.lineage = lineage_anno
                e_record.organism = tax_organism
                e_record.taxon_rank = t_hierarchy.resolved_to(lineage_anno)
                e_record.tracking_stamp()
        except KeyError:
            continue

    if len(pulled_tax_ids.symmetric_difference(set(tax_id_map.keys()))) > 0:
        dl_taxids = set(tax_id_map.keys())
        logging.debug("The following NCBI taxids are unique to the queries:\n{}\n"
                      "The following NCBI taxids are unique to the downloads:\n{}\n"
                      "".format(", ".join(dl_taxids.difference(pulled_tax_ids)),
                                ", ".join(pulled_tax_ids.difference(dl_taxids))))
    return


def entrez_records_to_accession_set(entrez_records_list: list, bitflag_filter=7):
    query_dict = dict()
    for record in entrez_records_list:  # type: EntrezRecord
        if record.bitflag > bitflag_filter:
            continue
        # Uses the accession ID if available, versioned accession otherwise
        if record.accession:
            query_dict[record.accession] = record
        elif record.versioned:
            query_dict[record.versioned] = record
        else:
            continue
    return query_dict


def entrez_records_from_lineages_and_chop(lineages: set, chopped_lineages: set, skip_dict: set):
    """
    With the taxonomic lineages in lineages, a set of unique lineages is created and added to the set chopped_lineages
    while the lineages is emptied (via popping the elements).

    :param lineages: A set of taxonomic lineages that need to be searched for in Entrez
    :param chopped_lineages: A set to add the lineages with their most resolved taxon removed
    :param skip_dict: A set containing organism names that don't need their lineages downloaded
    :return: A list of EntrezRecord instances
    """
    entrez_query_list = []
    unique_queries = set()
    while lineages:
        taxa = lineages.pop().split("; ")
        while len(taxa) > 1 and taxa[-1] in skip_dict:
            taxa = taxa[:-1]
        if taxa[-1] and taxa[-1] not in unique_queries:
            er = EntrezRecord("NA", "NA")  # Mock EntrezRecord, with no real accession or acc.version
            er.organism = taxa[-1]  # Set the organism to the most resolved taxon
            er.tracking_stamp()
            entrez_query_list.append(er)
            unique_queries.add(taxa[-1])
        if len(taxa) > 0:  # Don't add the chopped lineage if it was already the deepest element
            chopped_lineages.add("; ".join(taxa[:-1]))
    return entrez_query_list


def entrez_records_to_organism_set(entrez_records_list: list, bitflag_filter=7) -> dict:
    """
    Taking a list of EntrezRecord instances, the unique set of organism names is compiled and formatted for
    querying the Entrez database (requires '[All Names]' appended to it).
    Additionally, a filter (bitflag_filter) can be applied that will skip adding queries from EntrezRecords with a
    bitflag property value greater than the bitflag_filter.
    This is to skip over records that already have the relevant information.

    :param entrez_records_list: A list of EntrezRecord instances
    :param bitflag_filter: An integer used for skipping complete EntrezRecords
    :return: Dictionary with organism names mapped to EntrezRecords
    """
    query_dict = dict()
    for record in entrez_records_list:  # type: EntrezRecord
        if record.bitflag > bitflag_filter:
            continue

        if record.organism:
            record.organism = re.sub(r"[a-z]__", '', record.organism)
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


def fetch_taxids_from_organisms(search_terms: dict) -> None:
    """
    This function uses the keys of the dictionary to submit queries to the Entrez database, in hopes of retrieving
    NCBI taxonomy IDs (those unique, arbitrary numbers).
    For each record that is returned, all EntrezRecords mapped to record's organism name in search_terms are updated
    with the NCBI taxid (to their ncbi_tax variable) and their bitflag is recalculated.

    :param search_terms: A dictionary containing organism names as keys and lists of EntrezRecords as values
    :return: None
    """
    logging.debug(str(len(search_terms.keys())) + " unique organism queries.\n")
    logging.info("Retrieving NCBI taxonomy IDs for each organism... ")
    records_batch, durations, taxid_failures = tolerant_entrez_query(list(search_terms.keys()),
                                                                     "Taxonomy", "search", "xml", 1)
    logging.info("done.\n")

    for record in records_batch:
        try:
            organism = parse_gbseq_info_from_esearch_record(record, 'TranslationStack')['Term']
        except (IndexError, KeyError, TypeError):
            logging.debug("Value for 'TranslationStack' not found in Entrez record."
                          " It is likely this organism name doesn't exist in Entrez's taxonomy database.\n" +
                          "Unable to link taxonomy ID to organism.\nRecord:\n{}\n".format(record))
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
    return


def entrez_records_to_accession_lineage_map(entrez_records: list):
    # TODO: Remove this necessity. Currently need to reformat and tally accessions like so but its a waste
    # Used for tallying the status of Entrez queries
    success = 0
    rescued = 0
    bad_tax = 0
    bad_org = 0
    failed = 0
    accession_lineage_map = dict()

    for e_record in entrez_records:  # type: EntrezRecord
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


def get_multiple_lineages(entrez_query_list: list, t_hierarchy: TaxonomicHierarchy, molecule_type: str) -> None:
    """
    Function for retrieving taxonomic lineage information from accession IDs - accomplished in 3 steps:
     1. Query Entrez's Taxonomy database using accession IDs to obtain corresponding organisms
     2. Query Entrez's Taxonomy database using organism names to obtain corresponding TaxIds
     3. Query Entrez's Taxonomy database using TaxIds to obtain corresponding taxonomic lineages

    :param entrez_query_list: A list of EntrezRecord instances with accession IDs to be mapped to lineages
    :param t_hierarchy: A TaxonomicHierarchy instance
    :param molecule_type: The type of molecule (e.g. prot, nuc) to be mapped to a proper Entrez database name
    :return: None
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
        # Try to find the EntrezRecord matching the accession or accession.version
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
            logging.warning("Unable to map neither a record's accession nor accession.version to an EntrezRecord:\n" +
                            "Accession: '" + str(accession) + "'\n" +
                            "Acc.Version: '" + str(ver) + "'\n" +
                            "Alternatives:" + str(alt) + "\n" +
                            str(record) + "\n")
            continue
        e_record.organism = parse_gbseq_info_from_entrez_xml(record)
        e_record.tracking_stamp()

    ##
    # Step 2: Query Entrez's Taxonomy database using organism names to obtain corresponding taxonomic lineages
    ##
    o_search_term_map = entrez_records_to_organism_set(entrez_query_list, 3)
    fetch_taxids_from_organisms(o_search_term_map)

    fetch_lineages_from_taxids(entrez_query_list, t_hierarchy)

    return


def verify_lineage_information(accession_lineage_map: dict, entrez_record_map: dict,
                               t_hierarchy: TaxonomicHierarchy, taxa_searched: int) -> None:
    """
    Function used for parsing records returned by Bio.Entrez.efetch queries and identifying inconsistencies
    between the search terms and the results

    :param accession_lineage_map: A dictionary mapping accession.versionID tuples to taxonomic lineages
    :param entrez_record_map: A dictionary of EntrezRecord instances indexed by their unique TreeSAPP numerical IDs
    :param t_hierarchy: A TaxonomicHierarchy instance, that by this point should be fully populated
    :param taxa_searched: An integer for tracking number of accessions queried (currently number of lineages provided)
    :return: None
    """
    if (len(accession_lineage_map.keys()) + taxa_searched) != len(entrez_record_map):
        # Records were not returned for all sequences. Time to figure out which ones!
        logging.warning("Entrez did not return a record for every accession queried.\n"
                        "Don't worry, though. We'll figure out which ones are missing.\n")
    logging.debug("Entrez.efetch query stats:\n"
                  "\tDownloaded\t" + str(len(accession_lineage_map.keys())) + "\n" +
                  "\tProvided\t" + str(taxa_searched) + "\n" +
                  "\tTotal\t\t" + str(len(entrez_record_map)) + "\n\n")

    # Find the lineage searches that failed, add lineages to reference_sequences that were successfully identified
    for treesapp_id in sorted(entrez_record_map.keys()):
        ref_seq = entrez_record_map[treesapp_id]  # type: EntrezRecord
        ref_seq.tracking_stamp()
        # Could have been set previously, in custom header format for example
        if not ref_seq.lineage:
            lineage = ""
            for tuple_key in accession_lineage_map:
                accession, versioned = tuple_key
                if ref_seq.accession == accession or ref_seq.accession == versioned:
                    # if accession_lineage_map[tuple_key]["lineage"] == "":
                    #     lineage = "r__Unclassified"
                    # else:
                    #     # The query was successful! Add it and increment
                    lineage = accession_lineage_map[tuple_key]["lineage"]
                    if not ref_seq.organism and accession_lineage_map[tuple_key]["organism"]:
                        ref_seq.organism = accession_lineage_map[tuple_key]["organism"]

            if not lineage and ref_seq.bitflag == 0:
                logging.error("Lineage information was not retrieved for " + ref_seq.accession + "!\n" +
                              "Please remove the output directory and restart.\n")
                sys.exit(13)
            # elif not lineage and ref_seq.bitflag >= 1:
            #     lineage = "r__Unclassified"
        else:
            lineage = ref_seq.lineage

        ref_seq.lineage = t_hierarchy.check_lineage(lineage, ref_seq.organism)
        ref_seq.tracking_stamp()
        ref_seq.taxon_rank = t_hierarchy.resolved_to(ref_seq.lineage)
        if ref_seq.bitflag >= 1:
            taxa_searched += 1

    if taxa_searched < len(entrez_record_map.keys()):
        logging.error("Some sequences ({}/{}) were not used to query Entrez's taxonomy database!\n"
                      "".format(len(entrez_record_map)-taxa_searched, len(entrez_record_map)))
        sys.exit(9)

    return


def accession_lineage_map_from_entrez_records(ref_seq_map: dict) -> dict:
    """

    :param ref_seq_map: A dictionary of EntrezRecord instances indexed by their unique TreeSAPP numerical IDs
    :return: A dictionary mapping unique accessions to their taxonomic lineage. Will be written to accession_lineage_map
    """
    unambiguous_accession_lineage_map = dict()
    for treesapp_id in sorted(ref_seq_map.keys()):
        ref_seq = ref_seq_map[treesapp_id]  # type: EntrezRecord
        unambiguous_accession_lineage_map[ref_seq.accession] = ref_seq.lineage

    return unambiguous_accession_lineage_map


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
            logging.error("Accession '{}' present in {} multiple times!\n".format(accession, mapping_file))
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


def load_ref_seqs(fasta_dict: dict, header_registry: dict, ref_seq_dict: dict):
    """
    Function for adding sequences from a fasta-formatted dictionary into dictionary of ReferenceSequence objects

    :param fasta_dict: A fasta-formatted dictionary from a FASTA instance
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


def map_accessions_to_lineages(query_accession_list: list, t_hierarchy: TaxonomicHierarchy,
                               molecule: str, accession_to_taxid=None) -> None:
    if accession_to_taxid:
        # Determine find the query accessions that are located in the provided accession2taxid file
        entrez_record_dict = map_accession2taxid(query_accession_list, accession_to_taxid)
        entrez_records = []
        for index in entrez_record_dict:
            entrez_records += entrez_record_dict[index]
        # Map lineages to taxids for successfully-mapped query sequences
        fetch_lineages_from_taxids(entrez_records=entrez_records, t_hierarchy=t_hierarchy)
        # Use the normal querying functions to obtain lineage information for the unmapped queries
        unmapped_queries = pull_unmapped_entrez_records(entrez_records)
        if len(unmapped_queries) > 0:
            # This tends to be a minority so shouldn't be too taxing
            get_multiple_lineages(unmapped_queries, t_hierarchy, molecule)
            for e_record in unmapped_queries:  # type: EntrezRecord
                try:
                    entrez_record_dict[e_record.accession].append(e_record)
                except KeyError:
                    logging.warning(e_record.accession + " not found in original query list.\n")
                    continue
                entrez_records.append(e_record)
        entrez_record_dict.clear()
        unmapped_queries.clear()
    else:
        get_multiple_lineages(query_accession_list, t_hierarchy, molecule)
    return


class Lineage:
    def __init__(self, lin_sep="; "):
        self.Organism = None
        self.Lineage = None
        self.Domain = None
        self.Phylum = None
        self.Class = None
        self.Order = None
        self.Family = None
        self.Genus = None
        self.Species = None
        self.rank_attributes = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        self.lin_sep = lin_sep

    def verify_rank_occupancy(self) -> None:
        """
        Checks each rank in the lineage and ensures that a string representing a taxon exists.
        If a rank is empty, it truncates the lineage at the preceding rank (which is occupied).

        :return: None
        """
        if not self.Lineage:
            return

        i = 0
        ranks = self.Lineage.split(self.lin_sep)  # Create a list from which the lineage can be rebuilt
        while i < len(ranks):  # type: str
            if not ranks[i]:
                break
            pn = ranks[i].split('__')
            if len(pn) == 2:
                if len(pn[1]) == 0:
                    break
            i += 1

        while i < len(ranks):
            ranks.pop(i)
        self.Lineage = self.lin_sep.join(ranks)

        return

    def build_lineage(self, add_organism=False) -> str:
        self.domain_check()
        self.ensure_prefix()
        self.lineage_format_check()
        if not self.Lineage:
            taxa = []
            # Cut the lineage at the first empty rank
            for rank in self.rank_attributes:
                taxon = self.__dict__[rank]
                taxa.append(taxon)
            self.Lineage = self.lin_sep.join(taxa)

        self.verify_rank_occupancy()

        if not self.Lineage:
            logging.warning("Taxonomic lineage information was found in neither lineage nor taxonomic rank fields.\n")
            return ""

        if add_organism and self.Organism:
            if self.Lineage.split(self.lin_sep)[-1] != self.Organism:
                self.Lineage += self.lin_sep + self.Organism

        return self.Lineage

    def ensure_prefix(self) -> None:
        """
        Ensures the value for each rank attribute (taxon) has the appropriate prefix and prepends it missing
        For example, the prefix for self.Domain is 'd__' so a taxon might be 'd__Bacteria'

        :return: None
        """
        for rank in self.rank_attributes:
            prefix = rank.lower()[0] + "__"
            taxon = self.__dict__[rank]
            # Check the prefix if there is value for the attribute
            if taxon and not re.search(r"^" + re.escape(prefix), taxon):
                self.__dict__[rank] = prefix + taxon

        # Add a no_rank prefix to the Organism name
        if self.Organism:
            if not re.match(r"^[a-z]__.*", self.Organism):
                self.Organism = "n__" + self.Organism
        return

    def lineage_format_check(self) -> bool:
        """
        The Lineage.Lineage attribute's format is checked to determine whether:

1. it is populated
2. the separator matches the Lineage object's lineage separator

        :return: A boolean representing whether the lineage's format was modified
        """
        if not self.Lineage:
            return False

        try:
            # Test whether the lineage separator (self.lin_sep) exists in self.Lineage
            if self.Lineage.find(self.lin_sep) >= 0:
                return False
            else:
                # Does a separator exist?
                split_lin = re.split(r'[,|;"]+', self.Lineage)
                # TODO: Find the correct separator being used
                if len(split_lin) > 1:
                    # Replace the current separator with self.lin_sep
                    self.Lineage = self.lin_sep.join(split_lin)
                    return True
                else:
                    return False

        except ValueError:
            logging.error("Unable to split the lineage '{}' by the separator {}.\n".format(self.Lineage, self.lin_sep))
            raise ValueError()

    def domain_check(self) -> None:
        """
        Checks for whether or not the Domain is valid.
        Options are Bacteria, Archaea, Eukaryota and Viruses

        :return: None
        """
        if not self.Domain:
            return
        else:
            if re.sub(r"[a-z]__", '', self.Domain) not in ["Bacteria", "Archaea", "Eukaryota", "Viruses"]:
                logging.error("I've seen some taxonomic domains in my time and, friend, '{}' isn't one of them.\n"
                              "This may indicate a problem parsing your seqs2lineage file.\n".format(self.Domain))
                sys.exit(19)
            return


def read_seq_taxa_table(seq_names_to_taxa: str) -> dict:
    """
    Reads a table containing sequence names (of ORFs, contigs, etc.) in the first column and lineage information in the
    subsequent columns. Different types of information are indicated by their column name
    (you can name the first column what you want :)).

    Accepted columns names (case-insensitive):
    'Organism', 'Lineage', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'

    :param seq_names_to_taxa: Path to file containing the sequence name table
    :return: Dictionary with sequence names as keys and Lineage instances as values
    """
    seq_lineage_map = dict()
    sep = get_field_delimiter(seq_names_to_taxa)

    try:
        handler = open(seq_names_to_taxa, 'r', newline='')
    except IOError:
        logging.error("Unable to open '" + seq_names_to_taxa + "' for reading!\n")
        sys.exit(3)
    tbl_reader = csv.reader(handler, delimiter=sep)
    # Set up the named tuple that will be used for storing lineage information
    header_names = ["Organism", "Lineage", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    fields = next(tbl_reader)
    field_positions = get_list_positions(fields, header_names)
    if not field_positions:
        logging.error("Unable to read headers from sequence-taxa table file {}."
                      " The table must have some combination of the following column names:\n{}\n"
                      "".format(seq_names_to_taxa, ','.join(header_names)))
        sys.exit(3)

    try:
        for row in tbl_reader:
            seq_name = row[0]
            if seq_name[0] == '>':
                seq_name = seq_name[1:]
            seq_lin_info = Lineage()
            for field_name in field_positions:
                seq_lin_info.__dict__[field_name] = row[field_positions[field_name]].strip()
            seq_lineage_map[seq_name] = seq_lin_info
    except (csv.Error, IndexError) as e:
        logging.error("Reading file '{}', line {}:\n{}\n".format(seq_names_to_taxa, tbl_reader.line_num, e))
        sys.exit(3)

    handler.close()
    return seq_lineage_map


def map_orf_lineages(seq_lineage_tbl: str, header_registry: dict, refpkg_name=None) -> (dict, list):
    """
    The classified sequences have a signature at the end (|RefPkg_name|start_stop) that needs to be removed
    Iterates over the dictionary of sequence names and attempts to match those with headers in the registry.
    If a match is found the header is assigned the corresponding lineage in seq_lineage_map.

    :param seq_lineage_tbl: Path to file containing the sequence name table
    :param header_registry: A dictionary mapping numerical TreeSAPP identifiers to Header instances
    :param refpkg_name: The reference package's name
    :return: A dictionary mapping each classified sequence to a lineage and list of TreeSAPP IDs that were mapped
    """
    logging.info("Mapping assigned sequences to provided taxonomic lineages... ")
    seq_lineage_map = read_seq_taxa_table(seq_lineage_tbl)
    classified_seq_lineage_map = dict()
    treesapp_nums = list(header_registry.keys())
    mapped_treesapp_nums = []
    for seq_name, lineage in seq_lineage_map.items():  # type: (str, Lineage)
        # Its slow to perform so many re.search's but without having a guaranteed ORF pattern
        # we can't use hash-based data structures to bring it to O(N)
        parent_re = re.compile(seq_name)
        x = 0
        while x < len(treesapp_nums):
            header = header_registry[treesapp_nums[x]]
            original = header.original
            assigned_seq_name = re.sub(r"\|{0}\|\d+_\d+.*".format(refpkg_name), '', original)
            if parent_re.search(assigned_seq_name):
                classified_seq_lineage_map[header.first_split] = lineage.build_lineage(add_organism=True)
                mapped_treesapp_nums.append(treesapp_nums.pop(x))
            else:
                x += 1
        if len(treesapp_nums) == 0:
            logging.info("done.\n")
            return classified_seq_lineage_map, mapped_treesapp_nums
    logging.debug("Unable to find parent for " + str(len(treesapp_nums)) + " ORFs in sequence-lineage map:\n" +
                  "\n".join([header_registry[n].original for n in treesapp_nums]) + "\n")

    logging.info("done.\n")

    if len(mapped_treesapp_nums) == 0:
        logging.error("Unable to match any sequence names in {}.\n".format(seq_lineage_tbl))
        sys.exit(13)

    return classified_seq_lineage_map, mapped_treesapp_nums


def main():
    th = TaxonomicHierarchy()
    prep_for_entrez_query()
    tolerant_entrez_query(['12968'])
    er_vparadoxus = EntrezRecord(acc="WP_042579442", ver="WP_042579442.1")
    er_pmarinus = EntrezRecord(acc="WP_075487081", ver="WP_075487081.1")
    er_dict = {"1": er_vparadoxus, "2": er_pmarinus}
    get_multiple_lineages(list(er_dict.values()), th, "prot")
    alm = entrez_records_to_accession_lineage_map(list(er_dict.values()))
    repair_lineages(er_dict, th)
    verify_lineage_information(accession_lineage_map=alm, entrez_record_map=er_dict, t_hierarchy=th, taxa_searched=2)
    print(er_vparadoxus.get_info(),
          er_pmarinus.get_info())
    return


if __name__ == "__main__":
    main()

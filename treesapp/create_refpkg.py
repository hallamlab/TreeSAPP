#!/usr/bin/env python3

__author__ = "Connor Morgan-Lang"
__maintainer__ = "Connor Morgan-Lang"

try:
    import argparse
    import logging
    import sys
    import os
    import shutil
    import re
    import traceback

    from time import gmtime, strftime, sleep

    from . import utilities
    from . import fasta
    from . import classy
    from . import file_parsers
    from . import entrez_utils
    from .wrapper import run_odseq, run_mafft
    from .external_command_interface import launch_write_command
    from .entish import annotate_partition_tree
    from .lca_calculations import megan_lca, clean_lineage_list

except ImportError:
    sys.stderr.write("Could not load some user defined module functions:\n")
    sys.stderr.write(str(traceback.print_exc(10)))
    sys.exit(13)


def generate_cm_data(args, unaligned_fasta):
    """
    Using the input unaligned FASTA file:
     1. align the sequences using cmalign against a reference Rfam covariance model to generate a Stockholm file
     2. use the Stockholm file (with secondary structure annotated) to build a covariance model
     3. align the sequences using cmalign against a reference Rfam covariance model to generate an aligned fasta (AFA)
    :param args: command-line arguments objects
    :param unaligned_fasta:
    :return:
    """
    logging.info("Running cmalign to build Stockholm file with secondary structure annotations... ")

    cmalign_base = [args.executables["cmalign"],
                    "--mxsize", str(3084),
                    "--informat", "FASTA",
                    "--cpu", str(args.num_threads)]
    # First, generate the stockholm file
    cmalign_sto = cmalign_base + ["-o", args.code_name + ".sto"]
    cmalign_sto += [args.rfam_cm, unaligned_fasta]

    stdout, cmalign_pro_returncode = launch_write_command(cmalign_sto)

    if cmalign_pro_returncode != 0:
        logging.error("cmalign did not complete successfully for:\n" + ' '.join(cmalign_sto) + "\n")
        sys.exit(13)

    logging.info("done.\n")
    logging.info("Running cmbuild... ")

    # Build the CM
    cmbuild_command = [args.executables["cmbuild"]]
    cmbuild_command += ["-n", args.code_name]
    cmbuild_command += [args.code_name + ".cm", args.code_name + ".sto"]

    stdout, cmbuild_pro_returncode = launch_write_command(cmbuild_command)

    if cmbuild_pro_returncode != 0:
        logging.error("cmbuild did not complete successfully for:\n" +
                      ' '.join(cmbuild_command) + "\n")
        sys.exit(13)
    os.rename(args.code_name + ".cm", args.final_output_dir + os.sep + args.code_name + ".cm")
    if os.path.isfile(args.final_output_dir + os.sep + args.code_name + ".sto"):
        logging.warning("Overwriting " + args.final_output_dir + os.sep + args.code_name + ".sto\n")
        os.remove(args.final_output_dir + os.sep + args.code_name + ".sto")
    shutil.move(args.code_name + ".sto", args.final_output_dir)

    logging.info("done.\n")
    logging.info("Running cmalign to build MSA... ")

    # Generate the aligned FASTA file which will be used to build the BLAST database and tree with RAxML
    aligned_fasta = args.code_name + ".fa"
    cmalign_afa = cmalign_base + ["--outformat", "Phylip"]
    cmalign_afa += ["-o", args.code_name + ".phy"]
    cmalign_afa += [args.rfam_cm, unaligned_fasta]

    stdout, cmalign_pro_returncode = launch_write_command(cmalign_afa)

    if cmalign_pro_returncode != 0:
        logging.error("cmalign did not complete successfully for:\n" + ' '.join(cmalign_afa) + "\n")
        sys.exit(13)

    # Convert the Phylip file to an aligned FASTA file for downstream use
    seq_dict = file_parsers.read_phylip_to_dict(args.code_name + ".phy")
    fasta.write_new_fasta(seq_dict, aligned_fasta)

    logging.info("done.\n")

    return aligned_fasta


def create_new_ref_fasta(out_fasta, ref_seq_dict, dashes=False):
    """
    Writes a new FASTA file using a dictionary of ReferenceSequence class objects

    :param out_fasta: Name of the FASTA file to write to
    :param ref_seq_dict: Dictionary containing ReferenceSequence objects, numbers are keys
    :param dashes: Flag indicating whether hyphens should be retained in sequences
    :return:
    """
    out_fasta_handle = open(out_fasta, "w")
    num_seqs_written = 0

    for treesapp_id in sorted(ref_seq_dict, key=int):
        ref_seq = ref_seq_dict[treesapp_id]
        if dashes is False:
            sequence = re.sub('[-.]', '', ref_seq.sequence)
        else:
            # sequence = re.sub('\.', '', ref_seq.sequence)
            sequence = ref_seq.sequence
        out_fasta_handle.write(">" + ref_seq.short_id + "\n" + sequence + "\n")
        num_seqs_written += 1

    out_fasta_handle.close()

    if num_seqs_written == 0:
        logging.error("No sequences written to " + out_fasta + ".\n" +
                      "The headers in your input file are probably not accommodated in the regex patterns used. " +
                      "Function responsible: get_header_format. Please make an issue on the GitHub page.\n")
        sys.exit(5)

    return


def regenerate_cluster_rep_swaps(args, cluster_dict, fasta_replace_dict):
    """
    Function to regenerate the swappers dictionary with the original headers as keys and
    the new header (swapped in the previous attempt based on USEARCH's uc file) as a value
    :param args: command-line arguments objects
    :param cluster_dict: Dictionary where keys are centroid headers and values are headers of identical sequences
    :param fasta_replace_dict: Immature (lacking sequences) dictionary with header information parsed from tax_ids file
    :return:
    """
    swappers = dict()
    if args.verbose:
        sys.stderr.write("Centroids with identical sequences in the unclustered input file:\n")
    for rep in sorted(cluster_dict):
        matched = False
        subs = cluster_dict[rep]
        # If its entry in cluster_dict == 0 then there were no identical
        # sequences and the header could not have been swapped
        if len(subs) >= 1:
            # If there is the possibility the header could have been swapped,
            # check if the header is in fasta_replace_dict
            for treesapp_id in fasta_replace_dict:
                if matched:
                    break
                ref_seq = fasta_replace_dict[treesapp_id]
                # If the accession from the tax_ids file is the same as the representative
                # this one has not been swapped for an identical sequence's header since it is in use
                if re.search(ref_seq.accession, rep):
                    if args.verbose:
                        sys.stderr.write("\tUnchanged: " + rep + "\n")
                        matched = True
                    break
                # The original representative is no longer in the reference sequences
                # so it was replaced, with this sequence...
                for candidate in subs:
                    if rep in swappers or matched:
                        break

                    # parse the accession from the header
                    header_format_re, header_db, header_molecule = fasta.get_header_format(candidate, args.code_name)
                    sequence_info = header_format_re.match(candidate)
                    if sequence_info:
                        candidate_acc = sequence_info.group(1)
                    else:
                        logging.error("Unable to handle header: " + candidate + "\n")
                        sys.exit(13)

                    # Now compare...
                    if candidate_acc == ref_seq.accession:
                        if args.verbose:
                            sys.stderr.write("\tChanged: " + candidate + "\n")
                        swappers[rep] = candidate
                        matched = True
                        break
            sys.stderr.flush()
    return swappers


def finalize_cluster_reps(cluster_dict: dict, refseq_objects, header_registry):
    """
        Transfer information from the cluster data (representative sequence, identity and cluster taxonomic LCA) to the
    dictionary of ReferenceSequence objects. The sequences not representing a cluster will have their `cluster_rep`
    flags remain *False* so as to not be analyzed further.

    :param cluster_dict:
    :param refseq_objects:
    :param header_registry: A list of Header() objects, each used to map various header formats to each other
    :return: Dictionary of ReferenceSequence objects with complete clustering information
    """
    logging.debug("Finalizing representative sequence clusters... ")
    # Create a temporary dictionary mapping formatted headers to TreeSAPP numeric IDs
    tmp_dict = dict()
    for treesapp_id in header_registry:
        tmp_dict[header_registry[treesapp_id].original] = treesapp_id

    for cluster_id in sorted(cluster_dict, key=int):
        cluster_info = cluster_dict[cluster_id]
        treesapp_id = tmp_dict[cluster_info.representative]
        refseq_objects[treesapp_id].cluster_rep_similarity = '*'
        refseq_objects[treesapp_id].cluster_rep = True
        refseq_objects[treesapp_id].cluster_lca = cluster_info.lca

    logging.debug("done.\n")
    return refseq_objects


def present_cluster_rep_options(cluster_dict, refseq_objects, header_registry, important_seqs=None, each_lineage=False):
    """
    Present the headers of identical sequences to user for them to decide on representative header

    :param cluster_dict: dictionary from read_uc(uc_file)
    :param refseq_objects:
    :param header_registry: A list of Header() objects, each used to map various header formats to each other
    :param important_seqs: If --guarantee is provided, a dictionary mapping headers to seqs from format_read_fasta()
    :param each_lineage: If set to True, each candidate's lineage is shown as well as the cluster's LCA
    :return:
    """
    if not important_seqs:
        important_seqs = set()
    candidates = dict()
    for cluster_id in sorted(cluster_dict, key=int):
        cluster_info = cluster_dict[cluster_id]
        acc = 1
        candidates.clear()
        for num_id in sorted(refseq_objects, key=int):
            if header_registry[num_id].original == cluster_info.representative:
                refseq_objects[num_id].cluster_rep_similarity = '*'
                candidates[str(acc)] = refseq_objects[num_id]
                acc += 1
                break
        if acc != 2:
            raise AssertionError("Unable to find " + cluster_info.representative + " in ReferenceSequence objects!")

        if len(cluster_info.members) >= 1 and cluster_info.representative not in important_seqs:
            for cluster_member_info in cluster_info.members:
                for treesapp_id in sorted(refseq_objects, key=int):
                    formatted_header = header_registry[treesapp_id].original
                    if formatted_header == cluster_member_info[0]:
                        refseq_objects[treesapp_id].cluster_rep_similarity = cluster_member_info[1]
                        candidates[str(acc)] = refseq_objects[treesapp_id]
                        acc += 1
                        break
            sys.stderr.write("Sequences in '" + cluster_info.lca + "' cluster:\n")
            for num in sorted(candidates.keys(), key=int):
                sys.stderr.write("\t" + num + ". ")
                sys.stderr.write('\t'.join([candidates[num].organism + " | " + candidates[num].accession + "\t",
                                            str(len(candidates[num].sequence)) + "bp or aa",
                                            str(candidates[num].cluster_rep_similarity)]))
                if each_lineage:
                    sys.stderr.write("\t" + "(lineage = " + candidates[num].lineage + ")")
                sys.stderr.write("\n")
            sys.stderr.flush()

            best = input("Number of the best representative? ")
            # Useful for testing - no need to pick which sequence name is best!
            # best = str(1)
            while best not in candidates.keys():
                best = input("Invalid number. Number of the best representative? ")
            candidates[best].cluster_rep = True
            candidates[best].cluster_lca = cluster_info.lca
        else:
            refseq_objects[num_id].cluster_rep = True
            refseq_objects[num_id].cluster_lca = cluster_info.lca

    return refseq_objects


def reformat_headers(header_dict):
    """
    Imitate format_read_fasta header name reformatting
    :param header_dict: Dictionary of old header : new header key : value pairs
    :return:
    """
    swappers = dict()

    for old, new in header_dict.items():
        swappers[utilities.reformat_string(old)] = utilities.reformat_string(new)
    return swappers


def screen_filter_taxa(fasta_records: dict, screen_strs="", filter_strs="", guarantees=None) -> dict:
    """
    Searches the ReferenceSequence lineages in fasta_records for taxa that are either
    not meant to be retained (controlled by screen_strs) or supposed to be filtered out (controlled by filter_strs).

    :param fasta_records: A dictionary mapping `treesapp_id`s (integers) to ReferenceSequence objects
    :param screen_strs: The comma-separated string of taxa that should be RETAINED
    :param filter_strs: The comma-separated string of taxa that should be REMOVED
    :param guarantees: Optional set of numerical treesapp IDs that should not be removed from fasta_records
    :return: A new fasta_records dictionary with only the target taxa
    """
    fasta_replace_dict = dict()
    if not screen_strs and not filter_strs:
        return fasta_records

    if screen_strs:
        screen_terms = screen_strs.split(',')
    else:
        screen_terms = []
    if filter_strs:
        filter_terms = filter_strs.split(',')
    else:
        filter_terms = []

    num_filtered = 0
    num_screened = 0
    saved = set()
    for treesapp_id in fasta_records:
        screen_pass = False
        filter_pass = True
        ref_seq = fasta_records[treesapp_id]  # type: entrez_utils.ReferenceSequence
        # Screen
        if len(screen_terms) > 0:
            for term in screen_terms:
                # If any term is found in the lineage, it will pass... unless it fails the filter
                if re.search(term, ref_seq.lineage):
                    screen_pass = True
                    break
        else:
            screen_pass = True
        # Filter
        if len(filter_terms) > 0:
            for term in filter_terms:
                if re.search(term, ref_seq.lineage):
                    filter_pass = False

        if filter_pass and screen_pass:
            fasta_replace_dict[treesapp_id] = ref_seq
        elif guarantees and treesapp_id in guarantees:
            saved.add(treesapp_id)
            fasta_replace_dict[treesapp_id] = ref_seq
        else:
            if screen_pass is False:
                num_screened += 1
            if filter_pass is False:
                num_filtered += 1

    logging.debug('\t' + str(num_screened) + " sequences removed after failing screen.\n" +
                  '\t' + str(num_filtered) + " sequences removed after failing filter.\n" +
                  '\t' + str(len(fasta_replace_dict)) + " sequences retained.\n")
    if saved:
        logging.debug('\t' + str(len(saved)) + " guaranteed sequences saved from taxonomic filtering.\n")

    return fasta_replace_dict


def remove_by_truncated_lineages(fasta_records, min_taxonomic_rank, guarantees=None):
    rank_depth_map = {'k': 1, 'p': 2, 'c': 3, 'o': 4, 'f': 5, 'g': 6, 's': 7}
    min_depth = rank_depth_map[min_taxonomic_rank]
    if min_taxonomic_rank == 'k':
        return fasta_records

    num_removed = 0
    fasta_replace_dict = dict()

    for treesapp_id in fasta_records:
        ref_seq = fasta_records[treesapp_id]
        if guarantees and treesapp_id in guarantees:
            fasta_replace_dict[treesapp_id] = ref_seq
            continue
        if len(ref_seq.lineage.split("; ")) < min_depth:
            num_removed += 1
        elif re.search("^unclassified", ref_seq.lineage.split("; ")[min_depth-1], re.IGNORECASE):
            num_removed += 1
        else:
            fasta_replace_dict[treesapp_id] = ref_seq

    logging.debug('\t' + str(num_removed) + " sequences removed with truncated taxonomic lineages.\n" +
                  '\t' + str(len(fasta_replace_dict) - num_removed) + " sequences retained for building tree.\n")

    return fasta_replace_dict


def order_dict_by_lineage(fasta_object_dict):
    """
    Re-order the fasta_record_objects by their lineages (not phylogenetic, just alphabetical sort)
    Remove the cluster members since they will no longer be used

    :param fasta_object_dict: A dictionary mapping `treesapp_id`s (integers) to ReferenceSequence objects
    :return: An ordered, filtered version of the input dictionary
    """
    # Create a new dictionary with lineages as keys
    logging.debug("Re-enumerating the reference sequences in taxonomic order... ")
    lineage_dict = dict()
    sorted_lineage_dict = dict()
    for treesapp_id in fasta_object_dict:
        ref_seq = fasta_object_dict[treesapp_id]
        # Skip the redundant sequences that are not cluster representatives
        if not ref_seq.cluster_rep:
            continue
        try:
            lineage_dict[utilities.clean_lineage_string(ref_seq.lineage)].append(ref_seq)
        except KeyError:
            lineage_dict[utilities.clean_lineage_string(ref_seq.lineage)] = [ref_seq]

    # Now re-write the fasta_object_dict, but the numeric keys are now sorted by lineage
    #  AND it doesn't contain redundant fasta objects
    num_key = 1
    for lineage in sorted(lineage_dict.keys(), key=str):
        for ref_seq in lineage_dict[lineage]:
            if ref_seq.cluster_rep:
                # Replace the treesapp_id object
                code = '_'.join(ref_seq.short_id.split('_')[1:])
                ref_seq.short_id = str(num_key) + '_' + code
                sorted_lineage_dict[str(num_key)] = ref_seq
                num_key += 1

    logging.debug("done.\n")
    return sorted_lineage_dict


def threshold(lst, confidence="low"):
    """

    :param lst:
    :param confidence:
    :return:
    """
    if confidence == "low":
        # Majority calculation
        index = round(len(lst)*0.51)-1
    elif confidence == "medium":
        # >=75% of the list is reported
        index = round(len(lst)*0.75)-1
    else:
        # confidence is "high" and >=90% of the list is reported
        index = round(len(lst)*0.9)-1
    return sorted(lst, reverse=True)[index]


def estimate_taxonomic_redundancy(lineages: list):
    """
    Determines the lowest/shallowest taxonomic rank with reasonable coverage that should allow for quality assignments
    :param lineages: A list of taxonomic lineages
    :return: str
    """
    # TODO: Factor proximity of leaves in the tree into this measure
    # For instance, if the two or so species of the same genus are in the tree,
    # are they also beside each other in the same clade or are they located in different clusters?
    lowest_reliable_rank = "Species"
    rank_depth_map = {1: "Kingdoms", 2: "Phyla", 3: "Classes", 4: "Orders", 5: "Families", 6: "Genera", 7: "Species"}
    taxa_counts = dict()
    for depth in rank_depth_map:
        name = rank_depth_map[depth]
        taxa_counts[name] = dict()
    for lineage in lineages:
        taxa = lineage.split('; ')
        position = 1
        while position < len(taxa) and position < 8:
            if taxa[position] not in taxa_counts[rank_depth_map[position]]:
                taxa_counts[rank_depth_map[position]][taxa[position]] = 0
            taxa_counts[rank_depth_map[position]][taxa[position]] += 1
            position += 1
    for depth in rank_depth_map:
        rank = rank_depth_map[depth]
        redundancy = list()
        for taxon in taxa_counts[rank]:
            redundancy.append(taxa_counts[rank][taxon])
        if threshold(redundancy, "medium") == 1:
            lowest_reliable_rank = rank
            break

    logging.debug("Lowest reliable rank for taxonomic classification is: " + lowest_reliable_rank + "\n")

    return lowest_reliable_rank


def summarize_reference_taxa(reference_dict: dict, cluster_lca=False):
    """
    Function for enumerating the representation of each taxonomic rank within the finalized reference sequences
    :param reference_dict: A dictionary holding ReferenceSequence objects indexed by their unique numerical identifier
    :param cluster_lca: Boolean specifying whether a cluster's LCA should be used for calculation or not
    :return: A formatted, human-readable string stating the number of unique taxa at each rank
    """
    taxonomic_summary_string = ""
    # Not really interested in Cellular Organisms or Strains.
    rank_depth_map = {0: "Kingdoms", 1: "Phyla", 2: "Classes", 3: "Orders", 4: "Families", 5: "Genera", 6: "Species"}
    taxa_counts = dict()
    unclassifieds = 0

    for depth in rank_depth_map:
        name = rank_depth_map[depth]
        taxa_counts[name] = set()
    for num_id in sorted(reference_dict.keys(), key=int):
        if cluster_lca and reference_dict[num_id].cluster_lca:
            lineage = reference_dict[num_id].cluster_lca
        else:
            lineage = reference_dict[num_id].lineage

        if re.search("unclassified", lineage, re.IGNORECASE):
            unclassifieds += 1

        position = 0
        # Ensure the root/ cellular organisms designations are stripped
        taxa = utilities.clean_lineage_string(lineage).split('; ')
        while position < len(taxa) and position < 7:
            taxa_counts[rank_depth_map[position]].add(taxa[position])
            position += 1

    taxonomic_summary_string += "Number of unique lineages:\n"
    for depth in rank_depth_map:
        rank = rank_depth_map[depth]
        buffer = " "
        while len(rank) + len(str(len(taxa_counts[rank]))) + len(buffer) < 12:
            buffer += ' '
        taxonomic_summary_string += "\t" + rank + buffer + str(len(taxa_counts[rank])) + "\n"
    # Report number of "Unclassified" lineages
    taxonomic_summary_string += "Unclassified lineages account for " +\
                                str(unclassifieds) + '/' + str(len(reference_dict.keys())) + ' (' +\
                                str(round(float(unclassifieds*100)/len(reference_dict.keys()), 1)) + "%) references.\n"

    return taxonomic_summary_string


def write_tax_ids(fasta_replace_dict, tax_ids_file, taxa_lca=False):
    """
    Write the number, organism and accession ID, if possible
    :param fasta_replace_dict: Dictionary mapping numbers (internal treesapp identifiers) to ReferenceSequence objects
    :param tax_ids_file: The name of the output file
    :param taxa_lca: Flag indicating whether a cluster's lineage is just the representatives or the LCA of all members
    :return: Nothing
    """

    tree_taxa_string = ""
    warning_string = ""
    no_lineage = list()

    for treesapp_id in sorted(fasta_replace_dict.keys(), key=int):
        # Definitely will not uphold phylogenetic relationships but at least sequences
        # will be in the right neighbourhood rather than ordered by their position in the FASTA file
        reference_sequence = fasta_replace_dict[treesapp_id]
        if taxa_lca:
            lineage = reference_sequence.cluster_lca
        else:
            lineage = reference_sequence.lineage
        if not lineage:
            no_lineage.append(reference_sequence.accession)
            lineage = ''

        tree_taxa_string += "\t".join([str(treesapp_id),
                                      reference_sequence.organism + " | " + reference_sequence.accession,
                                       lineage]) + "\n"

    # Write the tree_taxa_string to the tax_ids file
    tree_tax_list_handle = open(tax_ids_file, "w")
    tree_tax_list_handle.write(tree_taxa_string)
    tree_tax_list_handle.close()

    if len(no_lineage) > 0:
        warning_string += str(len(no_lineage)) + " reference sequences did not have an associated lineage!\n\t"
        warning_string += "\n\t".join(no_lineage)

    return warning_string


def read_tax_ids(tree_taxa_list):
    """
    Reads the taxonomy and accession ID affiliated with each sequence number.
    This information is used to avoid horrible manual work if the pipeline is ran multiple times
    :param tree_taxa_list: The name of the tax_ids file to read
    :return:
    """
    try:
        tree_tax_list_handle = open(tree_taxa_list, 'r')
    except IOError:
        logging.error("Unable to open taxa list file '" + tree_taxa_list + "' for reading!\n")
        sys.exit(13)
    fasta_replace_dict = dict()
    line = tree_tax_list_handle.readline()
    while line:
        fields = line.strip().split("\t")
        if len(fields) == 3:
            treesapp_id, seq_info, lineage = fields
        else:
            treesapp_id, seq_info = fields
            lineage = ""
        ref_seq = entrez_utils.ReferenceSequence()
        try:
            ref_seq.organism = seq_info.split(" | ")[0]
            ref_seq.accession = seq_info.split(" | ")[1]
            ref_seq.lineage = lineage
        except IndexError:
            ref_seq.organism = seq_info
        fasta_replace_dict[treesapp_id] = ref_seq
        line = tree_tax_list_handle.readline()
    tree_tax_list_handle.close()

    return fasta_replace_dict


def update_build_parameters(param_file, marker_package: classy.MarkerBuild):
    """
    Function to update the data/tree_data/ref_build_parameters.tsv file with information on this new reference sequence
    Format of file is:
     "\t".join(["name","code","molecule","sub_model","cluster_identity","ref_sequences","tree-tool","poly-params",
     "lowest_reliable_rank","last_updated","description"])

    :param param_file: Path to the ref_build_parameters.tsv file used by TreeSAPP for storing refpkg metadata
    :param marker_package: A MarkerBuild instance
    :return: None
    """
    with open(param_file) as param_handler:
        param_lines = param_handler.readlines()
    # Read the the original build parameters file, extracting useful information from the old version if present
    updated_lines = []
    for line in param_lines:
        fields = line.strip("\n").split("\t")
        if fields[0] != marker_package.cog:
            updated_lines.append("\t".join(fields))
        else:
            if not marker_package.denominator or marker_package.denominator == "Z1111":
                marker_package.denominator = fields[1]
            if not marker_package.description:
                marker_package.description = fields[-1]

    # Add the new reference package's build parameters to the list of others
    marker_package.update = strftime("%d_%b_%Y", gmtime())
    build_list = [marker_package.cog, marker_package.denominator, marker_package.molecule, marker_package.model,
                  marker_package.kind, str(marker_package.pid), str(marker_package.num_reps), marker_package.tree_tool,
                  ','.join([str(param) for param in marker_package.pfit]),
                  marker_package.lowest_confident_rank, marker_package.update, marker_package.description]
    updated_lines.append("\t".join(build_list))

    try:
        params = open(param_file, 'w')
    except IOError:
        logging.error("Unable to open " + param_file + "for appending.\n")
        sys.exit(13)
    params.write("\n".join(updated_lines) + "\n")
    params.close()

    return


def parse_model_parameters(placement_trainer_file):
    """
    Returns the model parameters on the line formatted like
     'Regression parameters = (m,b)'
    in the file placement_trainer_results.txt
    :return: tuple
    """
    trainer_result_re = re.compile(r"^Regression parameters = \(([0-9,.-]+)\)$")
    try:
        trainer_handler = open(placement_trainer_file, 'r')
    except IOError:
        logging.error("Unable to open '" + placement_trainer_file + "' for reading!\n")
        sys.exit(3)
    params = None
    for line in trainer_handler:
        match = trainer_result_re.match(line)
        if match:
            params = match.group(1).split(',')
    trainer_handler.close()
    return params


def update_tax_ids_with_lineage(args, tree_taxa_list):
    tax_ids_file = args.treesapp + os.sep + "data" + os.sep + "tree_data" + os.sep + "tax_ids_%s.txt" % args.code_name
    if not os.path.exists(tax_ids_file):
        logging.error("Unable to find " + tax_ids_file + "!\n")
        raise FileNotFoundError
    else:
        fasta_replace_dict = read_tax_ids(tax_ids_file)
        # Determine how many sequences already have lineage information:
        lineage_info_complete = 0
        for treesapp_id in fasta_replace_dict:
            ref_seq = fasta_replace_dict[treesapp_id]
            if ref_seq.lineage:
                lineage_info_complete += 1
        # There are some that are already complete. Should they be over-written?
        if lineage_info_complete >= 1:
            overwrite_lineages = input(tree_taxa_list + " contains some sequences with complete lineages. "
                                                        "Should they be over-written? [y|n] ")
            while overwrite_lineages != "y" and overwrite_lineages != "n":
                overwrite_lineages = input("Incorrect response. Please input either 'y' or 'n'. ")
            if overwrite_lineages == 'y':
                ref_seq_dict = dict()
                for treesapp_id in fasta_replace_dict:
                    ref_seq = fasta_replace_dict[treesapp_id]
                    if ref_seq.lineage:
                        ref_seq.lineage = ""
                    ref_seq_dict[treesapp_id] = ref_seq
        # write_tax_ids(args, fasta_replace_dict, tax_ids_file, args.molecule)
    return


def remove_outlier_sequences(fasta_record_objects, od_seq_exe, mafft_exe, output_dir="./outliers", num_threads=2):
    od_input = output_dir + "od_input.fasta"
    od_output = output_dir + "outliers.fasta"
    outlier_names = list()
    tmp_dict = dict()

    outlier_test_fasta_dict = order_dict_by_lineage(fasta_record_objects)

    logging.info("Detecting outlier reference sequences... ")
    create_new_ref_fasta(od_input, outlier_test_fasta_dict)
    od_input_m = '.'.join(od_input.split('.')[:-1]) + ".mfa"
    # Perform MSA with MAFFT
    run_mafft(mafft_exe, od_input, od_input_m, num_threads)
    # Run OD-seq on MSA to identify outliers
    run_odseq(od_seq_exe, od_input_m, od_output, num_threads)
    # Remove outliers from fasta_record_objects collection
    outlier_seqs = fasta.read_fasta_to_dict(od_output)
    for seq_num_id in fasta_record_objects:
        ref_seq = fasta_record_objects[seq_num_id]
        tmp_dict[ref_seq.short_id] = ref_seq

    for seq_name in outlier_seqs:
        ref_seq = tmp_dict[seq_name]
        ref_seq.cluster_rep = False
        outlier_names.append(ref_seq.accession)

    logging.info("done.\n")
    logging.debug(str(len(outlier_seqs)) + " outlier sequences detected and discarded.\n\t" +
                  "\n\t".join([outseq for outseq in outlier_names]) + "\n")

    return fasta_record_objects


def guarantee_ref_seqs(cluster_dict: dict, important_seqs: set) -> dict:
    """
    Ensures all "guaranteed sequences" are representative sequences, swapping non-guaranteed sequences for the
    Cluster.representative where necessary.
    Also makes sure all guaranteed sequences are accounted for.
    :param cluster_dict:
    :param important_seqs:
    :return:
    """
    num_swaps = 0
    important_finds = set()
    nonredundant_guarantee_cluster_dict = dict()  # Will be used to replace cluster_dict
    expanded_cluster_id = 0
    for cluster_id in sorted(cluster_dict, key=int):
        cluster_inst = cluster_dict[cluster_id]  # type: classy.Cluster
        representative = cluster_inst.representative
        if len(cluster_inst.members) == 0:
            if representative in important_seqs:
                important_finds.add(representative)
            nonredundant_guarantee_cluster_dict[expanded_cluster_id] = cluster_inst
        else:
            contains_important_seq = False
            # The case where a member of a cluster is a guaranteed sequence, but not the representative
            x = 0
            while x < len(cluster_inst.members):
                member = cluster_inst.members[x]
                if member[0] in important_seqs:
                    nonredundant_guarantee_cluster_dict[expanded_cluster_id] = classy.Cluster(member[0])
                    nonredundant_guarantee_cluster_dict[expanded_cluster_id].members = []
                    nonredundant_guarantee_cluster_dict[expanded_cluster_id].lca = cluster_inst.lca
                    expanded_cluster_id += 1
                    contains_important_seq = True
                    important_finds.add(member[0])
                    cluster_inst.members.pop(x)
                else:
                    x += 1
            if representative in important_seqs:
                # So there is no opportunity for the important representative sequence to be swapped, clear members
                nonredundant_guarantee_cluster_dict[expanded_cluster_id] = cluster_inst
                important_finds.add(representative)
            elif contains_important_seq and representative not in important_seqs:
                num_swaps += 1
            else:
                nonredundant_guarantee_cluster_dict[expanded_cluster_id] = cluster_inst
        expanded_cluster_id += 1

    # Some final accounting - in case the header formats are altered!
    if important_seqs.difference(important_finds):
        logging.error(str(len(important_finds)) + '/' + str(len(important_seqs)) +
                      " guaranteed sequences found in cluster output file. The following are missing:\n" +
                      ", ".join([vis for vis in list(important_seqs.difference(important_finds))]) + "\n")
        sys.exit(7)
    logging.debug(str(num_swaps) + " former representative sequences were succeeded by 'guaranteed-sequences'.\n")

    return nonredundant_guarantee_cluster_dict


def cluster_lca(cluster_dict: dict, fasta_record_objects, header_registry: dict):
    # Create a temporary dictionary for faster mapping
    formatted_to_num_map = dict()
    for num_id in fasta_record_objects:
        formatted_to_num_map[header_registry[num_id].original] = num_id

    lineages = list()
    for cluster_id in sorted(cluster_dict, key=int):
        cluster_inst = cluster_dict[cluster_id]  # type: classy.Cluster
        members = [cluster_inst.representative]
        # format of member list is: [header, identity, member_seq_length/representative_seq_length]
        members += [member[0] for member in cluster_inst.members]
        # Create a lineage list for all sequences in the cluster
        for member in members:
            try:
                num_id = formatted_to_num_map[member]
                lineages.append(fasta_record_objects[num_id].lineage)
            except KeyError:
                logging.warning("Unable to map '" + str(member) + "' to a TreeSAPP numeric ID.\n")
        cleaned_lineages = clean_lineage_list(lineages)
        cluster_inst.lca = megan_lca(cleaned_lineages)
        # For debugging
        # if len(lineages) != len(cleaned_lineages) and len(lineages) > 1:
        #     print("Before:")
        #     for l in lineages:
        #         print(l)
        #     print("After:")
        #     for l in cleaned_lineages:
        #         print(l)
        #     print("LCA:", cluster_inst.lca)
        lineages.clear()
    formatted_to_num_map.clear()
    return

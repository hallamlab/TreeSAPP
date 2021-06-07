import re
import os
import sys
import logging

import pygtrie

from treesapp import logger
from treesapp import external_command_interface as eci
from treesapp import fasta
from treesapp import refpkg
from treesapp import phylo_seq

LOGGER = logging.getLogger(logger.logger_name())


def get_graftm_pkg_name(gpkg_path: str) -> str:
    pkg_name = re.sub(".gpkg", '', os.path.basename(gpkg_path))
    if not pkg_name:
        LOGGER.error("Unable to parse marker name from gpkg '{}'\n".format(gpkg_path))
        sys.exit(5)
    return pkg_name


def run_graftm_graft(input_path: str, output_dir: str, gpkg_path: str,
                     graftm_exe=None, classifier="graftm", seq_type="aminoacid", num_threads=4) -> None:
    """Wrapper for GraftM."""
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if not graftm_exe:
        graftm_exe = "graftM"

    classify_call = [graftm_exe, "graft",
                     "--forward", input_path,
                     "--graftm_package", gpkg_path,
                     "--input_sequence_type", seq_type,
                     "--threads", str(num_threads),
                     "--output_directory", output_dir,
                     "--force"]
    if classifier == "graftm":
        classify_call += ["--assignment_method", "pplacer"]
        classify_call += ["--search_method", "hmmsearch"]
    elif classifier == "diamond":
        classify_call += ["--assignment_method", "diamond"]
        classify_call += ["--search_method", "diamond"]

    LOGGER.debug("Command used:\n" + ' '.join(classify_call) + "\n")
    eci.launch_write_command(classify_call, False)

    return


def build_graftm_package(gpkg_path: str, tax_file: str, mfa_file: str, fa_file: str, threads: int):
    create_command = ["graftM", "create"]
    create_command += ["--threads", str(threads)]
    create_command += ["--alignment", mfa_file]
    create_command += ["--sequences", fa_file]
    create_command += ["--taxonomy", tax_file]
    create_command += ["--output", gpkg_path]
    create_command.append("--force")

    LOGGER.debug("Command used:\n" + ' '.join(create_command) + "\n")
    eci.launch_write_command(create_command, False)


def prep_graftm_ref_files(ref_pkg: refpkg.ReferencePackage, tmp_dir: str, target_clade: str, executables: dict) -> dict:
    """
    From the original TreeSAPP reference package files, the necessary GraftM create input files are generated
    with all reference sequences related to the target_taxon removed from the multiple sequence alignment,
    unaligned reference FASTA file and the tax_ids file.

    :param tmp_dir:  Path to write the intermediate files with target references removed
    :param target_clade: Taxonomic lineage of the clade that is being excluded from the reference package.
    :param ref_pkg: A ReferencePackage instance for the reference package being tested
    :param executables: Dictionary of paths to dependency executables indexed by their names. Must include:
         'hmmbuild', 'FastTree' and 'raxml-ng'.
    :return: A dictionary providing paths to output files
    """
    # GraftM refpkg input paths:
    output_paths = {"filtered_tax_ids": os.path.join(tmp_dir, ref_pkg.prefix + "_lineage_ids.txt"),
                    "filtered_mfa": os.path.join(tmp_dir, ref_pkg.prefix + ".mfa"),
                    "filtered_fasta": os.path.join(tmp_dir, ref_pkg.prefix + ".fa")}

    ce_refpkg = ref_pkg.clone(clone_path=tmp_dir + ref_pkg.prefix + ref_pkg.refpkg_suffix)

    ce_refpkg.exclude_clade_from_ref_files(tmp_dir=tmp_dir, executables=executables, target_clade=target_clade)

    # Write the lineage_ids file
    lineage_info = []
    for ref_leaf in ce_refpkg.generate_tree_leaf_references_from_refpkg():
        lineage_info.append("{}_{}\t{}".format(ref_leaf.number, ce_refpkg.prefix, ref_leaf.lineage))

    with open(output_paths["filtered_tax_ids"], 'w') as taxa_handler:
        taxa_handler.write("\n".join(lineage_info) + "\n")

    # Create and write the unaligned fasta file
    ce_fasta = ce_refpkg.get_fasta()  # type: fasta.FASTA
    fasta.write_new_fasta(fasta_dict=ce_fasta.fasta_dict, fasta_name=output_paths["filtered_mfa"])
    ce_fasta.unalign()
    fasta.write_new_fasta(fasta_dict=ce_fasta.fasta_dict, fasta_name=output_paths["filtered_fasta"])
    return output_paths


def grab_graftm_taxa(tax_ids_file) -> pygtrie.StringTrie:
    taxonomic_tree = pygtrie.StringTrie(separator='; ')
    with open(tax_ids_file) as tax_ids:
        header = tax_ids.readline().strip()
        last_rank = int(header[-1])
        final_index = 6 - last_rank
        if not re.search("parent_id,rank,tax_name,root,rank_0,rank_1,rank_2,rank_3,rank_4,rank_5,rank_6", header):
            LOGGER.error("Unable to handle format of " + tax_ids_file + "!")
            sys.exit(21)
        line = tax_ids.readline().strip()
        while line:
            fields = line.split(',')
            if final_index < 0:
                fields = line.split(',')[:final_index]
            try:
                _, _, _, _, _, k_, p_, c_, o_, f_, g_, s_, = fields
            except (IndexError, ValueError):
                LOGGER.error("Unexpected format of line with %d fields in " % len(line.split(',')) +
                             tax_ids_file + ":\n" + line)
                sys.exit(21)
            ranks = ["Root", k_, p_, c_, o_, f_, g_, s_]
            lineage_list = []
            # In case there are missing ranks... which is likely
            for rank in ranks:
                if rank:
                    # GraftM seems to append an 'e1' to taxa that are duplicated in the taxonomic lineage.
                    # For example: Bacteria; Aquificae; Aquificaee1; Aquificales
                    lineage_list.append(re.sub(r'_graftm_\d+$', '', rank))
                    # lineage_list.append(rank)
            lineage = re.sub('_', ' ', '; '.join(lineage_list))
            i = 0
            ranks = len(lineage)
            while i < len(lineage):
                taxonomic_tree["; ".join(lineage.split("; ")[:ranks - i])] = True
                i += 1

            line = tax_ids.readline().strip()
    return taxonomic_tree


def read_graftm_classifications(assignment_file) -> list:
    """
    Function for reading the _read_tax.tsv file generated by graftM.
    Sequences that have either multiple genes and/or subunits encoded or have homologous regions separated by
     a divergent sequence are recorded as "_split_" and only the first split is recorded and analyzed.

    :param assignment_file: Path to the _read_tax.tsv file
    :return: Dictionary indexed by taxonomic lineage whose values are headers of classified sequences
    """
    assignments = []
    try:
        assignments_handle = open(assignment_file, 'r')
    except IOError:
        LOGGER.error("Unable to open classification file '{}' for reading.\n".format(assignment_file))
        sys.exit(21)
    tax_lines = assignments_handle.readlines()
    assignments_handle.close()

    for line in tax_lines:
        fields = line.strip().split('\t')
        try:
            header, classified = fields
            if re.search("_split_", header):
                split = int(header.split('_')[-1])
                if split > 1:
                    continue
                else:
                    header = re.sub("_split_.*", '', header)
            classified = '; '.join([re.sub(r'e\d+$', '', taxon) for taxon in classified.split('; ')])
            if classified.split("; ")[0] == "Root":
                classified = "r__" + classified
            if header and classified:
                pquery = phylo_seq.PQuery()
                pquery.seq_name = header
                pquery.recommended_lineage = classified
                assignments.append(pquery)
        except ValueError:
            LOGGER.error("Unable to parse line:" + str(line))
            sys.exit(21)

    return assignments

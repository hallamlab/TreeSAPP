import re
import os
import sys
import logging

from treesapp import logger
from treesapp import external_command_interface as eci
from treesapp import fasta
from treesapp import refpkg

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

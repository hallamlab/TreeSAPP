"""
TreeSAPP command line.
"""
import sys
import argparse
import logging

import treesapp.commands as ts_commands
import treesapp.phylo_cluster as phyclust

usage = """
treesapp <command> [<args>]
** Commands include:
abundance      Calculate abundance measures for classified sequences 
assign         Classify query [protein|genomic] sequences using reference packages
colour         Colours a reference package's phylogeny based on taxonomic or phenotypic data
create         Create a reference package for a new gene, domain or orthologous group
evaluate       Evaluate the classification performance using clade exclusion analysis
layer          Layer extra annotation information on classifications with iTOL colours-style file(s)
package        Facilitate operations on reference packages
phylotu        Sort query sequences into clusters inferred from a reference package's phylogeny
purity         Characterize the sequences in a reference package using a curated database
train          Train a reference package's model used for correcting taxonomic rank classifications
update         Update an existing reference package with sequences found in a `classify` run
** Other commands:
info           Display TreeSAPP version and other information.
Use '-h' to get subcommand-specific help, e.g.
"""


def main(sys_args=None) -> int:
    commands = {
        "create": ts_commands.create,
        "evaluate": ts_commands.evaluate,
        "abundance": ts_commands.abundance,
        "assign": ts_commands.assign,
        "update": ts_commands.update,
        "info": ts_commands.info,
        "train": ts_commands.train,
        "colour": ts_commands.colour,
        "layer": ts_commands.layer,
        "purity": ts_commands.purity,
        "package": ts_commands.package,
        "phylotu": phyclust.cluster_phylogeny
    }
    parser = argparse.ArgumentParser(
        description="Phylogenetic classification of biological sequences"
    )
    parser.add_argument("command", nargs="?")
    if not sys_args:
        sys_args = sys.argv
    args = parser.parse_args(sys_args[1:2])
    if not args.command:
        sys.stderr.write(usage)
        sys.exit(1)

    if args.command not in commands:
        logging.error("Unrecognized command")
        sys.stderr.write(usage)
        sys.exit(1)

    cmd = commands.get(args.command)
    cmd(sys_args[2:])
    logging.info("TreeSAPP has finished successfully.\n")
    return 0


if __name__ == "__main__":
    main()

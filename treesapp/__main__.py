"""
TreeSAPP command line.
"""
import sys
import argparse
import logging

from .commands import (create, classify, update, info, train)
from .clade_exclusion_evaluator import main as evaluate_main

usage = """
treesapp <command> [<args>]
** Commands include:
create         Create a reference package for a new gene, domain or orthologous group
evaluate       Evaluate the classification performance using Clade-exclusion analysis
classify       Classify query [protein|genomic] sequences using reference packages 
update         Update an existing reference package with sequences found in a `classify` run
** Other commands:
info           Display TreeSAPP version and other information.
Use '-h' to get subcommand-specific help, e.g.
"""


def main():
    commands = {"create": create,
                "evaluate": evaluate_main,
                "classify": classify,
                "update": update,
                "info": info,
                "train": train}
    parser = argparse.ArgumentParser(description='Phylogenetic classification of biological sequences')
    parser.add_argument('command', nargs='?')
    args = parser.parse_args(sys.argv[1:2])

    if not args.command:
        sys.stderr.write(usage)
        sys.exit(1)

    if args.command not in commands:
        logging.error('Unrecognized command')
        sys.stderr.write(usage)
        sys.exit(1)

    cmd = commands.get(args.command)
    cmd(sys.argv[2:])


if __name__ == '__main__':
    main()

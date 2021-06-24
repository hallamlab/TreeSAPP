###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################
"""A Python package for gene-centric taxonomic and functional classification using phylogenetic placement."""
__author__ = "Connor Morgan-Lang, Ryan J. McLaughlin, Grace Zhang, Kevin Chan, Zachary Armstrong, Steven J. Hallam"
__author_email__ = "shallam@mail.ubc.ca"
__copyright__ = "Copyright 2020"
__credits__ = ["Connor Morgan-Lang"]
__description__ = "Python package for gene-centric taxonomic and functional classification using phylogenetic placement"
__license__ = "GPL-3.0"
__maintainer__ = "Connor Morgan-Lang"
__maintainer_email__ = "c.morganlang@gmail.com"
__python_requires__ = ">=3.6"
__status__ = "Production/Stable"
__title__ = "TreeSAPP"
__url__ = "https://github.com/hallamlab/TreeSAPP"
__version__ = "0.11.3"

__all__ = ['abundance', 'annotate_extra', 'assign', 'clade_annotation',
           'classy', 'commands', 'create_refpkg', 'entish', 'entrez_utils',
           'external_command_interface', 'fasta', 'file_parsers', 'hmmer_tbl_parser',
           'jplace_utils', 'lca_calculations', 'mcc_calculator',
           'phylo_cluster', 'phylo_dist', 'phylo_seq', 'phylogeny_painting',
           'placement_trainer', 'refpkg', 'rel_evo_dist', 'seq_clustering',
           'taxonomic_hierarchy', 'training_utils', 'treesapp_args', 'update_refpkg',
           'utilities', 'wrapper']

from . import (
    abundance,
    refpkg,
    fasta,
    classy,
    entrez_utils,
    taxonomic_hierarchy,
    lca_calculations
)

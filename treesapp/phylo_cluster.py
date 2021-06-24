#!/usr/bin/env python3

import os
import sys
import re
import glob
import logging

from Bio import Align
import numpy as np
from scipy import stats
from ete3 import Tree, TreeNode
from tqdm import tqdm

from treesapp.treesapp_args import TreeSAPPArgumentParser
from treesapp import file_parsers
from treesapp import jplace_utils
from treesapp import rel_evo_dist
from treesapp import refpkg
from treesapp import seq_clustering
from treesapp import phylo_seq
from treesapp import taxonomic_hierarchy as ts_taxonomy
from treesapp import classy as ts_classy
from treesapp import logger
from treesapp import fasta
from treesapp import wrapper
from treesapp import entish
from treesapp import utilities
from treesapp import training_utils
from treesapp import external_command_interface

LOGGER = logging.getLogger(logger.logger_name())


class PhylOTU:
    def __init__(self, name):
        """Initialize a PhylOTU instance."""
        self.number = name
        self.cardinality = 0
        self.edges = []
        self.taxon = None
        self.tree_node = None
        return


class PhyloClust(ts_classy.TreeSAPP):
    def __init__(self):
        """Initialize a PhyloClust instance."""
        super(PhyloClust, self).__init__("phylotu")
        # Parameters
        self.arg_parser = TreeSAPPArgumentParser(description="A tool for sorting query sequences placed on a phylogeny"
                                                             " into phylogenetically-inferred clusters.")
        self.clustering_mode = ""
        self.pre_mode = ""
        self.alpha = 0
        self.percentile = 0.99
        self.tax_rank = "species"
        self.normalize = False
        self.jplace_files = {}
        self.assign_output_dirs = {}
        self.source_paths = []
        self.source_type = ""
        self.output_dir = ""
        self.sample_re = None
        self.num_processes = 1
        self.clean = True

        # Objects for clustering
        self.clustered_pqueries = []
        self.sample_mat = {}
        self.cluster_index = {}
        self._edges_to_cluster_index = {}

        self.stages = {0: ts_classy.ModuleFunction("load_pqueries", 0),
                       1: ts_classy.ModuleFunction("phylogeny", 1),
                       2: ts_classy.ModuleFunction("cluster", 2)}

        return

    def load_args(self, args) -> None:
        # Set and create the output directory
        self.output_dir = args.output
        self.prep_log(args)
        self.set_output_dirs()
        self.check_previous_output(overwrite=False)

        self.validate_continue(args)

        self.clustering_mode = args.clustering_mode
        self.pre_mode = args.pre_mode

        if args.pkg_target:
            refpkg_dict = refpkg.gather_ref_packages(self.refpkg_dir, [args.pkg_target])
            self.ref_pkg = refpkg_dict[args.pkg_target]
        elif args.pkg_path:
            self.ref_pkg.f__pkl = args.pkg_path
            self.ref_pkg.slurp()
        else:
            self.ts_logger.error("A reference package must be provided to treesapp phylotu.\n")
            sys.exit(3)

        self.alpha = args.alpha
        self.tax_rank = args.tax_rank
        self.sample_re = re.compile(args.sample_regex)
        if args.jplace:
            self.source_type = "jplace"
            self.source_paths = args.jplace
        elif args.ts_out:
            self.source_type = "assign_dir"
            self.source_paths = args.ts_out
        self.executables = self.find_executables(args)
        self.num_processes = args.num_threads

        if args.delete is True:
            self.clean = True
        else:
            self.clean = False
        # Determine whether to normalise the evolutionary distances or not
        # if args.evo_dist == "red":
        #     self.normalize = True
        # elif args.evo_dist == "raw":
        #     self.normalize = False
        # else:
        #     self.ts_logger.error("Unexpected distance normalisation method: '{}'.\n".format(args.evo_dist))
        #     sys.exit(3)
        return

    def get_arguments(self, sys_args: list):
        """Add arguments to the TreeSAPP ArgumentParser instance, then parse and load."""
        refpkg_arg = self.arg_parser.reqs.add_mutually_exclusive_group()
        refpkg_arg.add_argument("--refpkg_path", dest="pkg_path",
                                help="Path to the reference package pickle (.pkl) file.\n")
        refpkg_arg.add_argument('--refpkg_name', dest="pkg_target",
                                help="Name of the reference package to use, referenced by its 'prefix' attribute, "
                                     "from the set packaged with TreeSAPP. "
                                     "Use `treesapp info -v` to get the available list.")
        self.arg_parser.add_output_dir()
        place_dat = self.arg_parser.reqs.add_mutually_exclusive_group()
        place_dat.add_argument("-j", "--jplace", nargs='+',
                               help="Path to one or more JPlace files generated by"
                                    " placement on a reference package's phylogeny.")
        place_dat.add_argument("--assign_output", nargs='+', dest="ts_out",
                               help="Path to one or more output directories of treesapp assign.")

        self.arg_parser.optopt.add_argument("-m", "--mode", dest="clustering_mode",
                                            choices=["de_novo", "ref_guided", "local"],
                                            default="ref_guided", required=False,
                                            help="The phylogentic clustering mode to use. [ DEFAULT = ref_guided ].")
        self.arg_parser.optopt.add_argument("-p", "--pre_cluster", dest="pre_mode",
                                            choices=["psc", "align"], default="psc", required=False,
                                            help="The method to use for pre-clustering the classified sequences, "
                                                 "based on either placement-space ('psc') or "
                                                 "pairwise alignment ('align'). [ DEFAULT = 'psc' ]")
        self.arg_parser.optopt.add_argument("-t", "--tax_rank",
                                            default="species", choices=["class", "order", "family", "genus", "species"],
                                            help="The taxonomic rank the cluster radius should approximately represent."
                                                 " [ DEFAULT = 'species' ].")
        self.arg_parser.optopt.add_argument("-a", "--alpha",
                                            default=0.0, required=False, type=float,
                                            help="The evolutionary distance threshold defining the cluster boundaries."
                                                 " [ DEFAULT = auto ].")
        self.arg_parser.optopt.add_argument("-s", "--sample_regex",
                                            default="", required=False, type=str,
                                            help="A regular expression for parsing the sample name from a query "
                                                 "sequence name. Example: '^(\d+)\.a:.*'. [ DEFAULT = None ].")
        self.arg_parser.add_compute_miscellany()
        self.arg_parser.add_delete()
        # TODO: Implement and validate these options
        # self.arg_parser.pplace_args.add_argument("-d", "--evo_dist",
        #                                          choices=["raw", "red"], default="red", required=False,
        #                                          help="The evolutionary distance normalisation method to use."
        #                                               " [ DEFAULT = red ]")

        # Parse the arguments
        return self.arg_parser.parse_args(sys_args)

    def prep_log(self, args):
        log_file = os.path.join(self.output_dir, "TreeSAPP_phyloclust_log.txt")
        logger.prep_logging(log_file, verbosity=args.verbose)
        return

    def announce_launch(self):
        banner = "\n##\t\t\t\tClustering sequences into pOTUs\t\t\t\t##\n\n"
        param_map = {"Taxonomic rank": self.tax_rank,
                     "Clustering mode": self.clustering_mode,
                     "Pre-clustering method": self.pre_mode,
                     "Alpha": self.alpha,
                     "Percentile": self.percentile,
                     "Output directory": self.output_dir}
        run_params = "PARAMETERS:\n" + "\n".join([k + "\t" + str(v) for k, v in param_map.items()])

        self.ts_logger.info(banner)
        self.ts_logger.debug("Arguments used:\n" + ' '.join(sys.argv[1:]) + "\n" +
                             run_params + "\n")
        return

    def load_sample_placement_files(self) -> None:
        # Format the JPlace files dictionary, mapping sample names to file paths
        if self.source_type == "jplace":
            if self.clustering_mode == "de_novo":
                self.ts_logger.error("TreeSAPP assign output directories must be provided for 'de_novo' clustering.\n")
                sys.exit(5)
            self.jplace_files = {os.path.basename(f_path).replace('.jplace', ''): f_path
                                 for f_path in self.source_paths}
        elif self.source_type == "assign_dir":
            for dirname in self.source_paths:  # type: str
                sample_id = os.path.basename(dirname.strip(os.sep))
                jplace_file = os.path.join(dirname, "iTOL_output", self.ref_pkg.prefix,
                                           self.ref_pkg.prefix + "_complete_profile.jplace")
                if not os.path.isfile(jplace_file):
                    self.ts_logger.warning("JPlace file was not found for {} in treesapp assign output path '{}'.\n"
                                           "".format(self.ref_pkg.prefix, dirname))
                    continue
                self.jplace_files[sample_id] = jplace_file
                self.assign_output_dirs[sample_id] = dirname
        self.source_paths.clear()
        return

    def clean_intermediate_files(self):
        """
        Removes all intermediate files, subdirectories and the 'intermediates' directory for phylotu
        """
        file_extensions = ["log", "mfa", "fasta", "nwk"]
        if self.clean:
            for order, stage in self.stages.items():  # type: (int, ts_classy.ModuleFunction)
                if os.path.isdir(stage.dir_path):
                    for ext in file_extensions:
                        for f_path in glob.glob(stage.dir_path + "*" + ext):
                            if os.path.isfile(f_path):
                                os.remove(f_path)
                    if len(os.listdir(stage.dir_path)) == 0:
                        os.rmdir(stage.dir_path)
            if len(os.listdir(self.var_output_dir)) == 0:
                os.rmdir(self.var_output_dir)
        return

    @staticmethod
    def check_monophyly(node: TreeNode, group_name: str) -> bool:
        if group_name in [t.prefix_taxon() for t in getattr(node, "taxon").lineage()]:
            return True
        return False

    @staticmethod
    def build_edge_node_index(ete_tree: Tree) -> dict:
        """Maps node.name attributes to edge integers in DFS traversal order."""
        node_edge_map = {}
        edge_name = 0
        if len(ete_tree.children) > 2:
            ete_tree.resolve_polytomy(recursive=False)
        ete_tree = ete_tree.get_tree_root()
        for node in ete_tree.traverse(strategy="postorder"):  # type: Tree
            if len(node.children) > 2:
                node.resolve_polytomy(recursive=False)
            node_edge_map[node.name] = edge_name
            edge_name += 1
        return node_edge_map

    def group_rel_dists(self, tree: Tree, hierarchy: ts_taxonomy.TaxonomicHierarchy) -> dict:
        """
        Determine the mean RED values between each pair of leaves of the same genus,
        thereby defining the expected range of values for a species.

        :param tree: A phylogenetic ete3 Tree instance
        :param hierarchy: A TaxonomicHierarchy instance with an 'accepted_ranks_depths' attribute mapping
         taxonomic ranks to their depth in the hierarchy, e.g. root -> 0, species -> 7
        :return: A dictionary mapping taxa to their group's monophyletic leaf-to-leaf distances
        """
        rel_dists = {}
        for node in tree.traverse():
            if not hasattr(node, "taxon") or not hasattr(node, "rel_dist"):
                self.ts_logger.error("Node {} is not complete with RED and taxonomy attributes.\n".format(node.name))
                sys.exit(17)

        entish.label_internal_nodes_ete(tree, attr="i_node", attr_type=int)
        group_index = hierarchy.accepted_ranks_depths[self.tax_rank]

        for r_leaf in tree.get_leaves():  # type: TreeNode
            r_tax = getattr(r_leaf, "taxon")  # type: ts_taxonomy.Taxon
            if len(r_tax.lineage()) <= group_index:
                continue
            try:
                target_group = r_tax.get_rank_in_lineage(rank_name=self.tax_rank).prefix_taxon()
            except IndexError:
                continue

            # Begin finding the query leaves related to the reference leaf
            for q_leaf in tree.get_leaves():  # type: TreeNode
                if q_leaf is r_leaf:
                    continue
                q_tax = getattr(q_leaf, "taxon")
                # Ensure the relationship between the two taxa is desired
                if ts_taxonomy.Taxon.lca(r_tax, q_tax).rank != self.tax_rank:
                    continue

                # Calculate the RED distance between the two monophyletic nodes
                ca = r_leaf.get_common_ancestor(q_leaf)
                if self.check_monophyly(ca, target_group):
                    pair_hash = utilities.elegant_pair(getattr(r_leaf, "i_node"),
                                                       getattr(q_leaf, "i_node"),
                                                       sort=True)
                    if self.normalize:
                        try:
                            rel_dists[target_group][pair_hash] = 1 - ca.rel_dist
                        except KeyError:
                            rel_dists[target_group] = {pair_hash: 1 - ca.rel_dist}
                    else:
                        try:
                            rel_dists[target_group][pair_hash] = r_leaf.get_distance(q_leaf)
                        except KeyError:
                            rel_dists[target_group] = {pair_hash: r_leaf.get_distance(q_leaf)}

        # Remove the hash index so map taxa to their branch length distances between leaves
        for target_group, pair_dists in rel_dists.items():
            rel_dists[target_group] = [round(x, 4) for x in pair_dists.values()]
        return rel_dists

    def get_phylogenetic_distances_for_rank(self, taxa_tree: Tree, taxonomy: ts_taxonomy.TaxonomicHierarchy,
                                            min_obs=3) -> np.array:
        dists = np.array([])
        # Find the 95% confidence interval for RED values between leaves of the same genus
        grouped_rel_dists = self.group_rel_dists(taxa_tree, taxonomy)
        if grouped_rel_dists:
            # Improve evenness across different taxa by augmenting branch length distances
            richest = max(len(x) for x in grouped_rel_dists.values())
            for _taxon, obs in grouped_rel_dists.items():
                if min_obs <= len(obs) < richest:
                    reps = int(richest / len(obs))
                    synth_obs = training_utils.augment_training_set(obs,
                                                                    n_reps=reps,
                                                                    feature_scale=np.std(obs)).reshape(len(obs) * reps,
                                                                                                       1)
                    obs = np.append(obs, synth_obs)
                dists = np.append(dists, obs)
        return dists

    def get_pairwise_distances_for_rank(self) -> np.array:
        dists = np.array([])
        taxon_leaf_map, _unique_taxa = self.ref_pkg.map_rank_representatives_to_leaves(self.tax_rank)
        ref_seqs = self.ref_pkg.get_fasta()  # type: fasta.FASTA
        ref_seqs.unalign()
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        for taxon, seq_names in taxon_leaf_map.items():
            r = seq_names.pop()  # type: str
            while seq_names:
                for q in seq_names:
                    alignments = aligner.align(ref_seqs.fasta_dict[r], ref_seqs.fasta_dict[q])
                    optimal = next(alignments)
                    dists = np.append(dists, [optimal.score/len(optimal.target)])
                r = seq_names.pop()
        return dists

    def calculate_distance_threshold(self, taxa_tree: Tree, taxonomy: ts_taxonomy.TaxonomicHierarchy,
                                     distribution="ci", override_rank=None) -> float:
        if override_rank is not None:
            self.stash("tax_rank")
            self.tax_rank = override_rank

        self.ts_logger.info("Estimating alpha threshold for '{}' phylogeny at rank '{}'... "
                            "".format(self.ref_pkg.prefix, self.tax_rank))

        min_obs = 3
        if self.clustering_mode in ["ref_guided", "de_novo"]:
            dists = self.get_phylogenetic_distances_for_rank(taxa_tree, taxonomy, min_obs)
        elif self.clustering_mode in ["local"]:
            dists = self.get_pairwise_distances_for_rank()
        else:
            dists = np.array([])

        self.withdraw_stash()
        if len(dists) > min_obs:
            if distribution == "ci":
                alpha = confidence_interval(dists, self.percentile)
            else:
                alpha = np.percentile(dists, 100 * self.percentile)
        else:
            alpha = -1

        self.ts_logger.info("done.\n")

        if alpha > 0:
            self.ts_logger.debug("Alpha estimated to be {}.\n".format(alpha))
        else:
            self.ts_logger.error("Unable to estimate alpha threshold for {}. "
                                 "Reference package size ({}) may be insufficient.\n"
                                 "".format(self.ref_pkg.prefix, self.ref_pkg.num_seqs))
        return alpha

    @staticmethod
    def partition_nodes(tree: Tree, alpha: float) -> dict:
        """
        Linear-time solution for the Max-diameter min-cut partitioning problem, as published in:
        Balaban, et al. (2019). TreeCluster: clustering biological sequences using phylogenetic trees.

        :param tree: The root of an ete3 Tree instance.
        :param alpha: The branch length threshold that delineates partitions
        :return: A dictionary mapping cluster names (integers) to TreeNode instances representing the root of clusters:
          d[1] -> TreeNode, d[2] -> TreeNode
        """
        node_partitions = {}
        i = 1
        for node in tree.traverse(strategy="postorder"):  # type: TreeNode
            if node.is_leaf():
                continue
            u_l, u_r = node.get_children()  # type: TreeNode
            b_ul, w_ul = node.get_distance(u_l), u_l.get_farthest_leaf()[1]
            b_ur, w_ur = node.get_distance(u_r), u_r.get_farthest_leaf()[1]
            if b_ul + w_ul + b_ur + w_ur > alpha:
                if b_ul + w_ul <= b_ur + w_ur:
                    node_partitions[i] = u_r
                    node.remove_child(u_r)
                    node.dist += b_ul + w_ul
                else:
                    node_partitions[i] = u_l
                    node.remove_child(u_l)
                    node.dist += b_ur + w_ur
                i += 1
            else:
                node.dist += max(b_ul + w_ul, b_ur + w_ur)

        # Capture the last leaf node, which should be the only one remaining in the tree
        if tree:
            node_partitions[i] = tree

        return node_partitions

    @staticmethod
    def split_node_partition_edges(alpha: float, node_partitions: dict) -> dict:
        """
        Prunes inner nodes that are part of a cluster when their accumulated distance exceeds alpha.
        Therefore, clusters may be entirely comprised of internal nodes.

        Distance is not accumulated within this function as this has been completed by self.partition_nodes()
        """
        inode_partitions = {}
        index = 1
        for subtree in node_partitions.values():  # type: TreeNode
            for node in subtree.traverse(strategy="postorder"):
                if len(subtree.children) == 0:
                    break
                if node is subtree:
                    break
                if node.dist > alpha:
                    inode_partitions[index] = node
                    node.up.dist -= node.dist
                    node.up.remove_child(node)
                    index += 1
            inode_partitions[index] = subtree
            index += 1
        return inode_partitions

    def set_pquery_leaf_node_taxa(self, tree: Tree) -> None:
        # Create a map between PQuery sequence names and taxon instances
        pquery_taxon_map = {}
        hierarchy = self.ref_pkg.taxa_trie  # type: ts_taxonomy.TaxonomicHierarchy
        for pquery in self.clustered_pqueries:  # type: phylo_seq.PQuery
            taxon_name = pquery.lct.split(hierarchy.lin_sep)[-1]
            pquery_taxon_map[pquery.place_name] = hierarchy.get_taxon(taxon_name)
            if not pquery_taxon_map[pquery.place_name]:
                pquery_taxon_map[pquery.place_name] = hierarchy.find_root_taxon()
        for leaf_node in tree.get_leaves():  # type: TreeNode
            try:
                leaf_node.taxon = pquery_taxon_map[leaf_node.name]
            except KeyError:
                continue  # this leaves the taxon element set to None

        self.ref_pkg.propagate_internal_node_taxon_labels_from_leaves(tree)
        return

    def stash(self, name: str) -> None:
        setattr(self, "_stash", {name: getattr(self, name)})
        return

    def withdraw_stash(self) -> None:
        if hasattr(self, "_stash"):
            k, v = getattr(self, "_stash").popitem()
            self.__dict__[k] = v
            delattr(self, "_stash")
        return

    def define_tree_clusters(self, tree: Tree, internal=True, override_alpha=None, verbose=False) -> None:
        if verbose:
            LOGGER.info("Defining phylogenetic clusters... ")

        if override_alpha is not None:
            self.stash("alpha")
            self.alpha = override_alpha

        cluster_nodes = self.partition_nodes(tree, self.alpha)
        if internal:
            # Split the cluster nodes into internal nodes when necessary
            cluster_nodes = self.split_node_partition_edges(self.alpha, cluster_nodes)

        chapter = len(self.cluster_index)
        for num, cluster_tree_node in cluster_nodes.items():  # type: (int, TreeNode)
            potu = PhylOTU(name=num + chapter)
            potu.tree_node = cluster_tree_node
            taxa = {getattr(leaf, "taxon") for leaf in cluster_tree_node.get_leaves() if leaf.taxon}
            if taxa:
                potu.taxon = taxa.pop()
            self.cluster_index[potu.number] = potu

        self.withdraw_stash()

        if verbose:
            LOGGER.info("done.\n")
        return

    def match_edges_to_clusters(self, tree: Tree, override_alpha=None) -> None:
        node_edge_map = self.build_edge_node_index(tree)
        if not self.cluster_index:
            self.define_tree_clusters(tree, override_alpha)
        for p_otu in self.cluster_index.values():  # type: PhylOTU
            for node in p_otu.tree_node.traverse():
                if node_edge_map[node.name] in self._edges_to_cluster_index:
                    self.ts_logger.error("Edge '{}' already associated with cluster {}.\n"
                                         "".format(node_edge_map[node.name],
                                                   self._edges_to_cluster_index[node_edge_map[node.name]]))
                    sys.exit(13)
                self._edges_to_cluster_index[node_edge_map[node.name]] = p_otu.number
        return

    def match_pqueries_to_clusters(self, pquery_precluster_map) -> dict:
        phylo_groups = {}
        centroid_pqueries = {c.representative: num for num, c in pquery_precluster_map.items()}
        for pquery in self.clustered_pqueries:  # type: phylo_seq.PQuery
            if pquery.place_name not in centroid_pqueries:
                continue
            cluster = self._edges_to_cluster_index[pquery.consensus_placement.edge_num]
            if cluster not in phylo_groups:
                phylo_groups[cluster] = {}
            phylo_groups[cluster].update({centroid_pqueries[pquery.place_name]:
                                              pquery_precluster_map[centroid_pqueries[pquery.place_name]]})
        return phylo_groups

    def add_pquery_count_to_sample_potu_mat(self, pquery: phylo_seq.PQuery) -> None:
        otu_id = getattr(pquery, "p_otu")
        if getattr(pquery, "sample_name") not in self.sample_mat:
            self.sample_mat[getattr(pquery, "sample_name")] = {n: 0 for n in self.cluster_index}
        try:
            self.sample_mat[getattr(pquery, "sample_name")][otu_id] += 1
        except KeyError:
            self.sample_mat[getattr(pquery, "sample_name")][otu_id] = 1
        return

    def assign_pqueries_to_leaf_clusters(self, pqueries: list, cluster_map: dict) -> None:
        # Build a map between the leaf node names and cluster numbers
        LOGGER.info("Indexing cluster members... ")
        leaf_cluster_map = {}
        for cluster_num, p_otu in self.cluster_index.items():  # type: (int, PhylOTU)
            for leaf in p_otu.tree_node.get_leaf_names():
                leaf_cluster_map[leaf] = cluster_num

        member_cluster_map = {}
        for _edge_num, cluster_idx in cluster_map.items():  # type: (int, dict)
            for _num, cluster in cluster_idx.items():
                member_cluster_map.update({member.place_name: cluster.representative for member in cluster.members})
        LOGGER.info("done.\n")

        p_bar = tqdm(total=len(pqueries), ncols=100, desc="Matching query sequences to clusters")
        for pquery in pqueries:  # type: phylo_seq.PQuery
            try:
                pquery.p_otu = leaf_cluster_map[pquery.place_name]
            except KeyError:
                pquery.p_otu = leaf_cluster_map[member_cluster_map[pquery.place_name]]
            self.add_pquery_count_to_sample_potu_mat(pquery)
            self.cluster_index[pquery.p_otu].cardinality += 1
            p_bar.update()
        p_bar.close()
        return

    def assign_pqueries_to_edge_clusters(self, pqueries: list) -> None:
        """
        Identifies which phylogenetic cluster each PQuery belongs to based on its edge_num attribute.
        This is saved in a new attribute called 'p_otu'.

        Resets the cardinality of each PhylOTU instance in self.cluster_index then recounts for this dataset.
        :param pqueries: A list of PQuery instances that were placed on the same reference tree the clusters are derived
        :return: None
        """
        # Reset cluster cardinality
        for p_otu in self.cluster_index.values():  # type: PhylOTU
            p_otu.cardinality = 0

        for pquery in pqueries:  # type: phylo_seq.PQuery
            cluster = self._edges_to_cluster_index[pquery.consensus_placement.edge_num]
            pquery.p_otu = cluster
            self.add_pquery_count_to_sample_potu_mat(pquery)
            self.cluster_index[cluster].cardinality += 1
        return

    def assign_pqueries_to_alignment_clusters(self, pqueries: dict, cluster_map: dict) -> None:
        pquery_cluster_map = {}
        for num, cluster in cluster_map.items():  # type: (str, seq_clustering.Cluster)
            pquery_cluster_map.update({place_name: int(num) for place_name, _id in cluster.members})
            self.cluster_index[int(num)] = PhylOTU(name=int(num))

        for sample_id in pqueries:
            if self.ref_pkg.prefix not in pqueries[sample_id]:
                continue
            for pquery in pqueries[sample_id][self.ref_pkg.prefix]:  # type: phylo_seq.PQuery
                pquery.p_otu = pquery_cluster_map[pquery.place_name]
                self.add_pquery_count_to_sample_potu_mat(pquery)
                self.cluster_index[pquery.p_otu].cardinality += 1
        return

    def set_pquery_sample_name(self, pquery: phylo_seq.PQuery, default_sample: str) -> None:
        """Uses a user-provided regular expression to parse the sample name from PQuery.seq_name attribute."""
        if self.sample_re.pattern:
            try:
                pquery.__setattr__("sample_name", self.sample_re.match(pquery.seq_name).group(1))
            except AttributeError:
                self.ts_logger.error("Regular expression {} did not match query sequence name '{}'.\n"
                                     "".format(self.sample_re.pattern, pquery.seq_name))
                sys.exit(5)
        else:
            pquery.__setattr__("sample_name", default_sample)
        return

    @staticmethod
    def open_output_file(file_path: str):
        try:
            f_handler = open(file_path, 'w')
        except IOError:
            LOGGER.error("Unable to open output file '{}' for writing.\n".format(file_path))
            sys.exit(13)
        return f_handler

    def write_cluster_taxon_labels(self) -> None:
        """Writes a two-column table mapping each cluster to its respective taxonomic lineage."""
        if self.clustering_mode == "local":
            return
        tbl_string = ""
        for cluster_num, p_otu in self.cluster_index.items():  # type: (int, PhylOTU)
            lineage = "; ".join([t.prefix_taxon() for t in p_otu.taxon.lineage()])
            tbl_string += "{}\t{}\n".format(cluster_num, lineage)

        taxa_tbl = self.open_output_file(os.path.join(self.final_output_dir, "phylotu_taxa.tsv"))
        taxa_tbl.write(tbl_string)
        taxa_tbl.close()
        return

    def write_otu_matrix(self, sep="\t") -> None:
        """Writes a typical OTU matrix with OTU IDs for rows and samples for columns."""
        otu_mat = self.open_output_file(os.path.join(self.final_output_dir, "phylotu_matrix.tsv"))

        # Write the header
        otu_mat.write(sep.join(["#OTU_ID"] + [sid for sid in sorted(self.sample_mat.keys())]) + "\n")
        buffer = []
        for cluster_num in sorted(self.cluster_index, key=int):
            buffer.append(cluster_num)
            for sid in sorted(self.sample_mat):
                buffer.append(self.sample_mat[sid][cluster_num])
            otu_mat.write(sep.join([str(x) for x in buffer]) + "\n")
            buffer.clear()

        otu_mat.close()

        return

    def write_pquery_otu_classifications(self, sep="\t") -> None:
        """Writes a tabular table mapping PQuery sequence names to their reference package phylogenetic OTUs."""
        pquery_otu_tbl = self.open_output_file(os.path.join(self.final_output_dir, "phylotu_pquery_assignments.tsv"))
        tbl_str = sep.join(
            ["PQuery", "RefPkg", "Resolution", "Mode", "pOTU", "d_Distal", "d_Pendant", "d_MeanTip"]) + "\n"
        for pquery in sorted(self.clustered_pqueries, key=lambda x: x.seq_name):  # type: phylo_seq.PQuery
            tbl_str += sep.join([pquery.seq_name,
                                 pquery.ref_name,
                                 self.tax_rank,
                                 self.clustering_mode,
                                 str(getattr(pquery, "p_otu")),
                                 str(pquery.consensus_placement.distal_length),
                                 str(pquery.consensus_placement.pendant_length),
                                 str(pquery.consensus_placement.mean_tip_length)]) + "\n"

        pquery_otu_tbl.write(tbl_str)
        pquery_otu_tbl.close()

        return

    def report_cluster_occupancy(self, sample_name, n_pqueries: int):
        occupancy = 0
        for _, p_otu in self.cluster_index.items():
            if p_otu.cardinality >= 1:
                occupancy += 1
        self.ts_logger.debug("{} classified sequences occupy {}/{} '{}' pOTUs in sample '{}'.\n"
                             "".format(n_pqueries, occupancy, len(self.cluster_index),
                                       self.ref_pkg.prefix, sample_name))

        return

    def summarise_cluster_sizes(self, clusters: dict):
        membership_sizes = []
        for _num, cluster in clusters.items():  # type: (int, seq_clustering.Cluster)
            membership_sizes.append(len(cluster.members))

        if self.clustering_mode == "local":
            cluster_state = "cluster"
        elif self.current_stage.name == "phylogeny":
            cluster_state = "pre-cluster"
        else:
            cluster_state = "cluster"

        self.ts_logger.info("Created {} {}(s) with an average size of {}.\n"
                            "".format(len(membership_sizes),
                                      cluster_state,
                                      round(sum(membership_sizes) / len(membership_sizes), 1)))
        return
    
    def gather_classified_pqueries(self) -> dict:
        classified_pqueries = {}
        for sample_id, ts_out in self.assign_output_dirs.items():
            assigner_instance = ts_classy.TreeSAPP("phylotu")
            classified_pqueries[sample_id] = file_parsers.load_classified_sequences_from_assign_output(ts_out,
                                                                                                       assigner_instance,
                                                                                                       self.ref_pkg.prefix)
            for pquery in classified_pqueries[sample_id][self.ref_pkg.prefix]:
                self.set_pquery_sample_name(pquery, sample_id)
                self.clustered_pqueries.append(pquery)
        return classified_pqueries


def confidence_interval(data, confidence=0.95) -> float:
    a = 1.0 * np.array(data)
    mu, se = np.mean(a), stats.sem(a)
    _h_b, h_u = stats.t.interval(confidence,
                                 len(a) - 1,
                                 loc=mu,
                                 scale=se)
    return h_u


def infer_cluster_phylogeny(fa_input: str, executables: dict, output_dir: str) -> tuple:
    fa_prefix, _ext = os.path.splitext(os.path.basename(fa_input))
    mfa_file = os.path.join(output_dir, fa_prefix + ".mfa")
    wrapper.run_mafft(mafft_exe=executables["mafft"],
                      fasta_in=fa_input,
                      fasta_out=mfa_file,
                      num_threads=2)
    clusters_tree = wrapper.construct_tree(tree_builder="FastTree",
                                           executables=executables,
                                           evo_model="LG",
                                           multiple_alignment_file=mfa_file,
                                           tree_output_dir=output_dir,
                                           tree_prefix=fa_prefix + "_tree",
                                           verbosity=0)

    # Ensure that the phylogeny is bifurcating
    ete_tree = Tree(clusters_tree)
    entish.label_internal_nodes_ete(ete_tree)
    return fa_input, ete_tree


def pqueries_to_fasta(pqueries, fa_name="") -> fasta.FASTA:
    pqueries_fasta = fasta.FASTA(file_name=fa_name)
    for _, refpkg_pquery_map in pqueries.items():
        for _, pqueries in refpkg_pquery_map.items():  # type: (str, list)
            for pquery in pqueries:  # type: phylo_seq.PQuery
                pqueries_fasta.fasta_dict[pquery.place_name] = pquery.seq
    pqueries_fasta.header_registry = fasta.register_headers(list(pqueries_fasta.fasta_dict.keys()))
    return pqueries_fasta


def subtree_finder(ref_tree: Tree, leaf_nodes: set, tree_size=3) -> Tree:
    """Finds the smallest clade that includes all lead node names in leaf_nodes and is at least min_size large."""
    subtree_root = None
    for subtree_root in ref_tree.traverse(strategy="postorder"):
        if leaf_nodes.issubset(subtree_root.get_leaf_names()):
            break

    while len(subtree_root) < tree_size:
        subtree_root = subtree_root.up
    new_root = subtree_root.detach().get_midpoint_outgroup()
    subtree_root.set_outgroup(new_root)

    # Prune to an even representation of the tree that is of length tree_size
    previous_dist = 0.0
    while len(subtree_root) > tree_size:
        branch_collapse_dist = sorted([n.dist for n in subtree_root])[-1 * tree_size]
        entish.collapse_ete_tree(subtree_root, branch_collapse_dist)
        if branch_collapse_dist == previous_dist:
            break
        else:
            previous_dist = branch_collapse_dist
    return subtree_root


def get_outgroup(tree: Tree, target: str) -> TreeNode:
    max_dist = 0.0
    outgroup = None
    target_leaf = tree.get_leaves_by_name(name=target).pop()
    for node in tree.get_leaves():
        if target_leaf.get_distance(node) > max_dist:
            outgroup = node
            max_dist = target_leaf.get_distance(node)
    return outgroup


def select_subtree_sequences(clusters: list, ref_pkg: refpkg.ReferencePackage, subtree_size=3) -> fasta.FASTA:
    """Find the reference sequences that are descendents of the placement edge and create a representative subset."""
    placement_edges = set()
    leaf_nodes = []
    for pquery in sum([c.members for c in clusters], []):  # type: phylo_seq.PQuery
        placement_edges.add(pquery.consensus_placement.edge_num)

    internal_node_leaves_map = ref_pkg.get_internal_node_leaf_map()
    for edge_num in placement_edges:
        leaf_nodes += internal_node_leaves_map[edge_num]

    sub_tree = subtree_finder(ref_pkg.get_ete_tree(), set(leaf_nodes), subtree_size)
    ingroup_leaf = sub_tree.get_closest_leaf()[0]  # type: TreeNode
    outgroup = get_outgroup(tree=ref_pkg.get_ete_tree(), target=ingroup_leaf.name)

    ref_leaf_sequences = ref_pkg.get_fasta()
    ref_leaf_sequences.keep_only(sub_tree.get_leaf_names() + outgroup.get_leaf_names())
    return ref_leaf_sequences


def remap_leaf_names(subtree_fasta: fasta.FASTA, de_novo_tree: Tree) -> None:
    """Adds a 'taxon' feature to all TreeNodes and replaces these Nodes' 'name' attribute with the sequence name."""
    # Propagate a 'taxon' feature - None by default - to all TreeNodes for holding Taxon instances
    for n in de_novo_tree.traverse(strategy="postorder"):  # type: Tree
        n.add_feature(pr_name="taxon", pr_value=None)
    # Add the query sequence name for each leaf
    for leaf_node in de_novo_tree.get_leaves():
        leaf_node.name = subtree_fasta.header_registry[leaf_node.name].original
    return


def generate_trees_for_preclusters(phylo_clust: PhyloClust, pre_clusters: dict, pqueries_fasta: fasta.FASTA) -> None:
    task_list = []
    fa_index = {}
    for edge_num, cluster_num in pre_clusters.items():  # type: (int, dict)
        ref_subtree_fasta = select_subtree_sequences(list(cluster_num.values()), phylo_clust.ref_pkg)
        ref_subtree_fasta.file = "{}_E{}.fasta".format(phylo_clust.ref_pkg.prefix, edge_num)
        ref_subtree_fasta.unalign()
        # Combine the query sequences from the cluster and reference sequences

        cluster_fasta = fasta.subset_fasta(pqueries_fasta, [c.representative for c in list(cluster_num.values())])
        ref_subtree_fasta.fasta_join(cluster_fasta)

        # Change header format to one without whitespace
        ref_subtree_fasta.change_dict_keys("num")
        # Write a FASTA file for multiple sequence alignment and phylogenetic inference
        fa_input = os.path.join(phylo_clust.stage_output_dir, os.path.basename(ref_subtree_fasta.file))
        fasta.write_new_fasta(ref_subtree_fasta.fasta_dict,
                              fasta_name=fa_input)
        fa_index[fa_input] = ref_subtree_fasta

        task_list.append([fa_input, phylo_clust.executables, phylo_clust.stage_output_dir])

    phylogenies = external_command_interface.run_apply_async_multiprocessing(func=infer_cluster_phylogeny,
                                                                             arguments_list=task_list,
                                                                             num_processes=phylo_clust.num_processes,
                                                                             pbar_desc="Inferring phylogenies for each partition")
    p_bar = tqdm(desc="Defining clusters for subtrees", total=len(phylogenies), ncols=100)
    for tree_map in phylogenies:  # type: tuple
        fa_file, de_novo_tree = tree_map  # type: (str, Tree)
        remap_leaf_names(fa_index[fa_file], de_novo_tree)
        # Use the taxonomy of each classified sequence as leaf node taxon attribute
        phylo_clust.set_pquery_leaf_node_taxa(tree=de_novo_tree)
        # Define tree clusters
        phylo_clust.define_tree_clusters(tree=de_novo_tree, internal=False)
        p_bar.update()
    p_bar.close()

    # Remove clusters only occupied by reference sequences from PhyloClust.cluster_index
    pquery_names = set([pq.place_name for pq in phylo_clust.clustered_pqueries])
    for num in list(phylo_clust.cluster_index):  # type: int
        if not pquery_names.intersection(phylo_clust.cluster_index[num].tree_node.get_leaf_names()):
            phylo_clust.cluster_index.pop(num)
    return


def format_precluster_map(cluster_method: str, precluster_map: dict) -> dict:
    batch_indexed_cluster_map = {}
    if cluster_method == "psc":
        # Index the preclusters by their placement edge
        for num, cluster in precluster_map.items():  # type: (str, seq_clustering.Cluster)
            rep = cluster.members.__getitem__(0)  # type: phylo_seq.PQuery
            if rep.consensus_placement.edge_num not in batch_indexed_cluster_map:
                batch_indexed_cluster_map[rep.consensus_placement.edge_num] = {}
            batch_indexed_cluster_map[rep.consensus_placement.edge_num][num] = cluster
    elif cluster_method == "align":
        # Combine so all representatives are used to infer a single phylogeny in generate_trees_for_preclusters
        batch_indexed_cluster_map[0] = precluster_map
    return batch_indexed_cluster_map


def de_novo_phylo_clusters(phylo_clust: PhyloClust, taxon_labelled_tree: Tree,
                           cluster_method=None, phylo_group="partition", psc_cluster_size=10, drep_id=1.0) -> None:
    # Gather all classified sequences for a reference package from the treesapp assign outputs
    classified_pqueries = phylo_clust.gather_classified_pqueries()
    if len(classified_pqueries) == 0:
        return

    # Load the classified PQuery sequences into a FASTA instance
    pqueries_fasta = pqueries_to_fasta(classified_pqueries,
                                       fa_name=os.path.join(phylo_clust.stage_output_dir, "classified_pqueries.fasta"))
    pqueries_fasta.summarize_fasta_sequences()
    phylo_clust.increment_stage_dir()

    # Dereplicate classified sequences into Clusters
    pquery_precluster_map = {}
    if cluster_method == "align":
        pquery_precluster_map = seq_clustering.dereplicate_by_clustering(fasta_inst=pqueries_fasta,
                                                                         prop_similarity=drep_id,
                                                                         mmseqs_exe=phylo_clust.executables["mmseqs"],
                                                                         tmp_dir=phylo_clust.stage_output_dir)
        for cluster in pquery_precluster_map.values():  # type: seq_clustering.Cluster
            member_names = [x[0] for x in cluster.members]
            cluster.members = [pq for pq in phylo_clust.clustered_pqueries if pq.place_name in member_names]
    elif cluster_method == "psc":
        pquery_precluster_map = phylo_seq.cluster_pquery_placement_space_distances(
            pqueries=phylo_clust.clustered_pqueries,
            min_cluster_size=psc_cluster_size)
        # Keep only the representative sequences in the FASTA instance
        pqueries_fasta.change_dict_keys()
        pqueries_fasta.keep_only(header_subset=[c.representative for num, c in pquery_precluster_map.items()])

    phylo_clust.summarise_cluster_sizes(pquery_precluster_map)
    if phylo_group == "partition":
        LOGGER.info("Grouping representative query sequences into clades for tree inference... ")
        broader_rank = phylo_clust.ref_pkg.taxa_trie.deeper_rank(phylo_clust.tax_rank)
        # Broadly partition reference phylogeny for phylogenetic trees
        phylo_clust.ts_logger.disabled = True
        alpha = phylo_clust.calculate_distance_threshold(taxon_labelled_tree,
                                                         phylo_clust.ref_pkg.taxa_trie,
                                                         override_rank=broader_rank)
        phylo_clust.ts_logger.disabled = False
        phylo_clust.match_edges_to_clusters(tree=phylo_clust.ref_pkg.taxonomically_label_tree(),
                                            override_alpha=alpha)
        # Find the minimal set of clusters that contain all pqueries.
        phylo_groups = phylo_clust.match_pqueries_to_clusters(pquery_precluster_map)
        LOGGER.info("done.\n")
    else:
        # Combine clusters into batches based on the pre-cluster method used
        phylo_groups = format_precluster_map(cluster_method, pquery_precluster_map)

    # Infer a phylogeny and clusters for each pre-cluster (from either pairwise alignment or placement-space)
    generate_trees_for_preclusters(phylo_clust, phylo_groups, pqueries_fasta)

    phylo_clust.increment_stage_dir()
    for sample_id in classified_pqueries:
        if phylo_clust.ref_pkg.prefix not in classified_pqueries[sample_id]:
            continue
        pqueries = classified_pqueries[sample_id][phylo_clust.ref_pkg.prefix]
        # Use pquery_precluster_map to assign pqueries to clusters that are not in the phylogeny
        phylo_clust.assign_pqueries_to_leaf_clusters(pqueries, phylo_groups)

        phylo_clust.report_cluster_occupancy(sample_name=sample_id, n_pqueries=len(pqueries))
    return


def cluster_by_local_alignment(phylo_clust: PhyloClust, proportional_identity=0.0) -> None:
    if proportional_identity == 0.0:
        proportional_identity = phylo_clust.ref_pkg.pid

    # Gather all classified sequences for a reference package from the treesapp assign outputs
    classified_pqueries = phylo_clust.gather_classified_pqueries()
    if len(classified_pqueries) == 0:
        return
    # Load the classified PQuery sequences into a FASTA instance
    pqueries_fasta = pqueries_to_fasta(classified_pqueries,
                                       fa_name=os.path.join(phylo_clust.stage_output_dir, "classified_pqueries.fasta"))
    pqueries_fasta.summarize_fasta_sequences()
    phylo_clust.increment_stage_dir()

    # Dereplicate classified sequences into Clusters
    cluster_map = seq_clustering.dereplicate_by_clustering(fasta_inst=pqueries_fasta,
                                                           prop_similarity=proportional_identity,
                                                           mmseqs_exe=phylo_clust.executables["mmseqs"],
                                                           tmp_dir=phylo_clust.stage_output_dir)

    phylo_clust.summarise_cluster_sizes(cluster_map)
    phylo_clust.increment_stage_dir()
    phylo_clust.assign_pqueries_to_alignment_clusters(classified_pqueries, cluster_map)

    for sample_id in classified_pqueries:
        phylo_clust.report_cluster_occupancy(sample_name=sample_id,
                                             n_pqueries=len(classified_pqueries[sample_id][phylo_clust.ref_pkg.prefix]))
    return


def reference_placement_phylo_clusters(phylo_clust: PhyloClust, taxon_labelled_tree: Tree,
                                       psc=False, psc_cluster_size=10) -> None:
    # Perform maximum, mean, or median distance min-cut partitioning
    phylo_clust.define_tree_clusters(phylo_clust.ref_pkg.taxonomically_label_tree(), verbose=True)

    # Find the edges belonging to each cluster and invert the dictionary to create the _edges_to_cluster_index
    phylo_clust.match_edges_to_clusters(tree=phylo_clust.ref_pkg.taxonomically_label_tree())

    if psc:
        phylo_clust.cluster_index = phylo_seq.cluster_pquery_placement_space_distances(
            pqueries=phylo_clust.clustered_pqueries,
            min_cluster_size=psc_cluster_size)

    LOGGER.info("Assigning query sequences to clusters:\n")
    p_bar = tqdm(total=len(phylo_clust.jplace_files), ncols=100)
    for sample_id, jplace_file in phylo_clust.jplace_files.items():
        p_bar.set_description(desc="'{}'".format(sample_id), refresh=True)
        # Load the JPlace file
        jplace_dat = jplace_utils.jplace_parser(jplace_file)

        pqueries = jplace_utils.demultiplex_pqueries(jplace_dat)
        for pquery in pqueries:  # type: phylo_seq.PQuery
            phylo_clust.set_pquery_sample_name(pquery, sample_id)
            pquery.process_max_weight_placement(taxon_labelled_tree)
            phylo_clust.clustered_pqueries.append(pquery)

        # Map the PQueries (max_lwr or aelw?) to clusters
        phylo_clust.assign_pqueries_to_edge_clusters(pqueries)

        # Report the number of clusters occupied
        phylo_clust.report_cluster_occupancy(sample_name=sample_id, n_pqueries=len(pqueries))
        p_bar.update()
    p_bar.close()
    return


def cluster_phylogeny(sys_args: list) -> None:
    p_clust = PhyloClust()
    p_clust.load_args(p_clust.get_arguments(sys_args))
    p_clust.announce_launch()
    p_clust.load_sample_placement_files()
    if len(p_clust.jplace_files) == 0 and len(p_clust.assign_output_dirs) == 0:
        return

    # Calculate RED distances for each node
    taxa_tree = p_clust.ref_pkg.taxonomically_label_tree()
    red_tree = rel_evo_dist.RedTree()
    red_tree.decorate_rel_dist(taxa_tree)

    if not p_clust.alpha:
        p_clust.alpha = p_clust.calculate_distance_threshold(taxa_tree, p_clust.ref_pkg.taxa_trie)
        if p_clust.alpha < 0:
            return

    if p_clust.clustering_mode == "ref_guided":
        reference_placement_phylo_clusters(phylo_clust=p_clust,
                                           taxon_labelled_tree=taxa_tree)
    elif p_clust.clustering_mode == "de_novo":
        de_novo_phylo_clusters(phylo_clust=p_clust,
                               cluster_method=p_clust.pre_mode,
                               taxon_labelled_tree=taxa_tree)
    elif p_clust.clustering_mode == "local":
        cluster_by_local_alignment(phylo_clust=p_clust, proportional_identity=p_clust.alpha)
    else:
        LOGGER.error("Unknown clustering mode specified: '{}'.\n".format(p_clust.clustering_mode))
        sys.exit(5)

    # Exit if there were no sequences classified
    if not p_clust.clustered_pqueries:
        LOGGER.warning("No sequences were classified as '{}' in any of the provided treesapp assign outputs.\n"
                       "".format(p_clust.ref_pkg.prefix))
        return

    # Write a OTU table with the abundance of each PhylOTU for each JPlace
    p_clust.write_otu_matrix()

    # Write a table mapping each OTU identifier to its taxonomy
    p_clust.write_cluster_taxon_labels()

    # Write a table mapping PQuery sequence names to their cluster
    p_clust.write_pquery_otu_classifications()

    p_clust.clean_intermediate_files()

    # TODO: Centroids? Options are discussed in TreeCluster
    return


if __name__ == "__main__":
    cluster_phylogeny(sys.argv[1:])

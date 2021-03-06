#!/usr/bin/env python3

import os
import sys
import re
import logging

import numpy as np
from ete3 import Tree, TreeNode

from treesapp.treesapp_args import TreeSAPPArgumentParser
from treesapp.refpkg import ReferencePackage
from treesapp import file_parsers
from treesapp import jplace_utils
from treesapp import rel_evo_dist
from treesapp import phylo_seq
from treesapp import taxonomic_hierarchy as ts_taxonomy
from treesapp import classy as ts_classy


class PhylOTU:
    def __init__(self, name):
        """Initialize a PhylOTU instance."""
        self.number = name
        self.cardinality = 0
        self.edges = []
        self.taxon = None
        self.tree_node = None
        return


class PhyloClust:
    def __init__(self):
        """Initialize a PhyloClust instance."""
        # Parameters
        self.arg_parser = TreeSAPPArgumentParser(description="A tool for sorting query sequences placed on a phylogeny"
                                                             " into phylogenetically-inferred clusters.")
        self.treesapp_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.sep
        self.refpkg_dir = self.treesapp_dir + 'data' + os.sep
        self.alpha = 0
        self.tax_rank = "species"
        self.normalize = False
        self.jplace_files = {}
        self.output_dir = ""
        self.sample_re = None

        # Objects for clustering
        self.refpkg = ReferencePackage()
        self.sample_mat = {}
        self.cluster_index = {}
        self._edges_to_cluster_index = {}
        return

    def load_args(self, args) -> None:
        # Set and create the output directory
        self.output_dir = args.output
        if not os.path.isdir(self.output_dir):
            try:
                os.mkdir(self.output_dir)
            except IOError:
                logging.error("Unable to make output directory '{}'.\n".format(self.output_dir))
                sys.exit(5)

        self.prep_log(args)

        if args.pkg_target:
            refpkg_dict = file_parsers.gather_ref_packages(self.refpkg_dir, [args.pkg_target])
            self.refpkg = refpkg_dict[args.pkg_target]
        elif args.pkg_path:
            self.refpkg.f__json = args.pkg_path
            self.refpkg.slurp()
        else:
            logging.error("A reference package must be provided to treesapp phylotu.\n")
            sys.exit(3)

        self.alpha = args.alpha
        self.tax_rank = args.tax_rank
        self.sample_re = re.compile(args.sample_regex)

        # Format the JPlace files dictionary, mapping sample names to file paths
        if args.jplace:
            self.jplace_files = {os.path.basename(f_path).replace('.jplace', ''): f_path for f_path in args.jplace}
        elif args.ts_out:
            for dirname in args.ts_out:  # type: str
                sample_id = os.path.basename(dirname.strip(os.sep))
                jplace_file = os.path.join(dirname, "iTOL_output", self.refpkg.prefix,
                                           self.refpkg.prefix + "_complete_profile.jplace")
                if not os.path.isfile(jplace_file):
                    logging.warning("JPlace file was not found for {} in treesapp assign output path '{}'.\n"
                                    "".format(self.refpkg.prefix, dirname))
                    continue
                self.jplace_files[sample_id] = jplace_file

        # Determine whether to normalise the evolutionary distances or not
        # if args.evo_dist == "red":
        #     self.normalize = True
        # elif args.evo_dist == "raw":
        #     self.normalize = False
        # else:
        #     logging.error("Unexpected distance normalisation method: '{}'.\n".format(args.evo_dist))
        #     sys.exit(3)
        return

    def get_arguments(self, sys_args: list) -> None:
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
        # TODO: Implement and validate these options
        # self.arg_parser.pplace_args.add_argument("-d", "--evo_dist",
        #                                          choices=["raw", "red"], default="red", required=False,
        #                                          help="The evolutionary distance normalisation method to use."
        #                                               " [ DEAULT = red ]")
        # self.arg_parser.optopt.add_argument("-m", "--partition_metric",
        #                                     choices=["mean", "median", "max"], default="max", required=False,
        #                                     help="The metric to use when deciding when to cut a cluster."
        #                                          " [ DEFAULT = max ].")

        # Parse the arguments
        self.load_args(self.arg_parser.parse_args(sys_args))
        return

    def prep_log(self, args):
        log_file = os.path.join(self.output_dir, "TreeSAPP_phyloclust_log.txt")
        ts_classy.prep_logging(log_file, verbosity=args.verbose)
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
                logging.error("Node {} is not complete with RED and taxonomy attributes.\n".format(node.name))
                sys.exit(17)

        group_index = hierarchy.accepted_ranks_depths[self.tax_rank]

        for r_leaf in tree.get_leaves():  # type: TreeNode
            r_tax = getattr(r_leaf, "taxon")  # type: ts_taxonomy.Taxon
            if len(r_tax.lineage()) <= group_index:
                continue
            try:
                target_group = r_tax.prefix_taxon()  # type: str
            except IndexError:
                continue

            # Begin finding the query leaves related to the reference leaf
            for q_leaf in tree.get_leaves():  # type: TreeNode
                if q_leaf is r_leaf:
                    continue
                q_lineage = getattr(q_leaf, "taxon").lineage()
                if len(q_lineage) <= group_index:
                    continue
                # Ensure the query leaf belongs to the same taxonomic group
                if target_group not in [t.prefix_taxon() for t in q_lineage]:
                    continue

                # Calculate the RED distance between the two monophyletic nodes
                ca = r_leaf.get_common_ancestor(q_leaf)
                if self.check_monophyly(ca, target_group):
                    if self.normalize:
                        try:
                            rel_dists[target_group].add(1-ca.rel_dist)
                        except KeyError:
                            rel_dists[target_group] = {1-ca.rel_dist}
                    else:
                        try:
                            rel_dists[target_group].add(r_leaf.get_distance(q_leaf))
                        except KeyError:
                            rel_dists[target_group] = {r_leaf.get_distance(q_leaf)}

        return rel_dists

    def calculate_distance_threshold(self, taxa_tree: Tree, taxonomy: ts_taxonomy.TaxonomicHierarchy) -> None:
        # Find the 95% confidence interval for RED values between leaves of the same genus
        grouped_rel_dists = self.group_rel_dists(taxa_tree, taxonomy)
        dists = []
        for group_dists in grouped_rel_dists.values():
            dists += group_dists
        self.alpha = np.percentile(dists, 95)
        return

    def partition_nodes(self, tree: Tree) -> dict:
        """
        Linear-time solution for the Max-diameter min-cut partitioning problem, as published in:
        Balaban, et al. (2019). TreeCluster: clustering biological sequences using phylogenetic trees.

        :param tree: The root of an ete3 Tree instance.
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
            if b_ul + w_ul + b_ur + w_ur > self.alpha:
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

    def split_node_partition_edges(self, node_partitions: dict) -> dict:
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
                if node.dist > self.alpha:
                    inode_partitions[index] = node
                    node.up.dist -= node.dist
                    node.up.remove_child(node)
                    index += 1
            inode_partitions[index] = subtree
            index += 1
        return inode_partitions

    def define_tree_clusters(self, tree: Tree) -> None:
        node_partitions = self.partition_nodes(tree)
        # Split the cluster nodes into internal nodes when necessary
        cluster_nodes = self.split_node_partition_edges(node_partitions)

        for num, cluster_tree_node in cluster_nodes.items():  # type: (int, TreeNode)
            potu = PhylOTU(name=num)
            potu.tree_node = cluster_tree_node
            potu.taxon = getattr(cluster_tree_node, "taxon")
            self.cluster_index[num] = potu
        return

    def match_edges_to_clusters(self, tree: Tree) -> None:
        node_edge_map = self.build_edge_node_index(tree)
        if not self.cluster_index:
            self.define_tree_clusters(tree)
        for p_otu in self.cluster_index.values():  # type: PhylOTU
            for node in p_otu.tree_node.traverse():
                if node_edge_map[node.name] in self._edges_to_cluster_index:
                    logging.error("Edge '{}' already associated with cluster {}.\n"
                                  "".format(node_edge_map[node.name],
                                            self._edges_to_cluster_index[node_edge_map[node.name]]))
                    sys.exit(13)
                self._edges_to_cluster_index[node_edge_map[node.name]] = p_otu.number
        return

    def assign_pqueries_to_clusters(self, pqueries: list) -> None:
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
            if getattr(pquery, "sample_name") not in self.sample_mat:
                self.sample_mat[getattr(pquery, "sample_name")] = {n: 0 for n in self.cluster_index}
            try:
                self.sample_mat[getattr(pquery, "sample_name")][cluster] += 1
            except KeyError:
                self.sample_mat[getattr(pquery, "sample_name")][cluster] = 1
            self.cluster_index[cluster].cardinality += 1
        return

    def set_pquery_sample_name(self, pquery: phylo_seq.PQuery, default_sample: str) -> None:
        """Uses a user-provided regular expression to parse the sample name from PQuery.seq_name attribute."""
        if self.sample_re.pattern:
            try:
                pquery.__setattr__("sample_name", self.sample_re.match(pquery.seq_name).group(1))
            except AttributeError:
                logging.error("Regular expression {} did not match query sequence name '{}'.\n"
                              "".format(self.sample_re.pattern, pquery.seq_name))
                sys.exit(5)
        else:
            pquery.__setattr__("sample_name", default_sample)
        return

    def write_cluster_taxon_labels(self) -> None:
        """Writes a two-column table mapping each cluster to its respective taxonomic lineage."""
        tbl_string = ""
        output_table = os.path.join(self.output_dir, "phylotu_taxa.tsv")
        for cluster_num, p_otu in self.cluster_index.items():  # type: (int, PhylOTU)
            lineage = "; ".join([t.prefix_taxon() for t in p_otu.taxon.lineage()])
            tbl_string += "{}\t{}\n".format(cluster_num, lineage)

        with open(output_table, 'w') as taxa_tbl:
            taxa_tbl.write(tbl_string)
        return

    def write_otu_matrix(self, sep="\t") -> None:
        """Writes a typical OTU matrix with OTU IDs for rows and samples for columns."""
        try:
            otu_mat = open(os.path.join(self.output_dir, "phylotu_matrix.tsv"), 'w')
        except IOError:
            logging.error("Unable to open output file '{}' for writing.\n")
            sys.exit(13)

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


def cluster_phylogeny(sys_args: list) -> None:
    p_clust = PhyloClust()
    p_clust.get_arguments(sys_args)

    # Calculate RED distances for each node
    taxa_tree = p_clust.refpkg.taxonomically_label_tree()
    red_tree = rel_evo_dist.RedTree()
    red_tree.decorate_rel_dist(taxa_tree)

    if not p_clust.alpha:
        logging.info("Alpha threshold not provided. Estimating for '{}'... ".format(p_clust.tax_rank))
        p_clust.calculate_distance_threshold(taxa_tree, p_clust.refpkg.taxa_trie)
        logging.info("done.\n")
        logging.debug("Alpha estimated to be {}.\n".format(p_clust.alpha))

    # Perform maximum, mean, or median distance min-cut partitioning
    logging.info("Defining phylogenetic clusters... ")
    p_clust.define_tree_clusters(p_clust.refpkg.taxonomically_label_tree())
    logging.info("done.\n")

    # Find the edges belonging to each cluster and invert the dictionary to create the _edges_to_cluster_index
    p_clust.match_edges_to_clusters(tree=p_clust.refpkg.taxonomically_label_tree())

    logging.info("Assigning query sequences to clusters...\n")
    for sample_id, jplace_file in p_clust.jplace_files.items():
        # Load the JPlace file
        jplace_dat = jplace_utils.jplace_parser(jplace_file)

        pqueries = jplace_utils.demultiplex_pqueries(jplace_dat)
        for pquery in pqueries:  # type: phylo_seq.PQuery
            p_clust.set_pquery_sample_name(pquery, sample_id)
            pquery.process_max_weight_placement(taxa_tree)

        # Map the PQueries (max_lwr or aelw?) to clusters
        p_clust.assign_pqueries_to_clusters(pqueries)

        # Report the number of clusters occupied
        occupancy = 0
        for num, p_otu in p_clust.cluster_index.items():
            if p_otu.cardinality >= 1:
                occupancy += 1
        logging.info("{} classified sequences occupy {}/{} '{}' pOTUs in sample '{}'.\n"
                     "".format(len(pqueries), occupancy, len(p_clust.cluster_index), p_clust.refpkg.prefix, sample_id))

    # Write a OTU table with the abundance of each PhylOTU for each JPlace
    p_clust.write_otu_matrix()

    # Write a table mapping each OTU identifier to its taxonomy
    p_clust.write_cluster_taxon_labels()

    # TODO: Return pqueries with a 'cluster' attribute

    # TODO: Decide what to do with the PQueries with pendant lengths greater than the cluster RED radius

    # TODO: Centroids? Options are discussed in TreeCluster
    return


if __name__ == "__main__":
    cluster_phylogeny(sys.argv[1:])

#!/usr/bin/env python3

import os
import sys
import re
import logging

import numpy as np
from ete3 import Tree, TreeNode

from treesapp.treesapp_args import TreeSAPPArgumentParser
from treesapp import file_parsers
from treesapp import jplace_utils
from treesapp import rel_evo_dist
from treesapp import phylo_seq
from treesapp import taxonomic_hierarchy as ts_taxonomy
from treesapp import classy as ts_classy
from treesapp import fasta
from treesapp import wrapper
from treesapp import entish


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
        self.alpha = 0
        self.tax_rank = "species"
        self.normalize = False
        self.jplace_files = {}
        self.assign_output_dirs = {}
        self.output_dir = ""
        self.sample_re = None

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

        if args.pkg_target:
            refpkg_dict = file_parsers.gather_ref_packages(self.refpkg_dir, [args.pkg_target])
            self.ref_pkg = refpkg_dict[args.pkg_target]
        elif args.pkg_path:
            self.ref_pkg.f__json = args.pkg_path
            self.ref_pkg.slurp()
        else:
            logging.error("A reference package must be provided to treesapp phylotu.\n")
            sys.exit(3)

        self.alpha = args.alpha
        self.tax_rank = args.tax_rank
        self.sample_re = re.compile(args.sample_regex)

        # Format the JPlace files dictionary, mapping sample names to file paths
        if args.jplace:
            if self.clustering_mode == "de_novo":
                logging.error("TreeSAPP assign output directories must be provided for 'de_novo' clustering.\n")
                sys.exit(5)
            self.jplace_files = {os.path.basename(f_path).replace('.jplace', ''): f_path for f_path in args.jplace}
        elif args.ts_out:
            for dirname in args.ts_out:  # type: str
                sample_id = os.path.basename(dirname.strip(os.sep))
                jplace_file = os.path.join(dirname, "iTOL_output", self.ref_pkg.prefix,
                                           self.ref_pkg.prefix + "_complete_profile.jplace")
                if not os.path.isfile(jplace_file):
                    logging.warning("JPlace file was not found for {} in treesapp assign output path '{}'.\n"
                                    "".format(self.ref_pkg.prefix, dirname))
                    continue
                self.jplace_files[sample_id] = jplace_file
                self.assign_output_dirs[sample_id] = dirname

        self.executables = self.find_executables(args)
        # Determine whether to normalise the evolutionary distances or not
        # if args.evo_dist == "red":
        #     self.normalize = True
        # elif args.evo_dist == "raw":
        #     self.normalize = False
        # else:
        #     logging.error("Unexpected distance normalisation method: '{}'.\n".format(args.evo_dist))
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
                                            choices=["de_novo", "ref_guided"], default="ref_guided", required=False,
                                            help="The phylogentic clustering mode to use. [ DEFAULT = ref_guided ].")
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

        # Parse the arguments
        return self.arg_parser.parse_args(sys_args)

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

    def set_pquery_leaf_node_taxa(self, tree: Tree) -> None:
        # Create a map between PQuery sequence names and taxon instances
        pquery_taxon_map = {}
        hierarchy = self.ref_pkg.taxa_trie  # type: ts_taxonomy.TaxonomicHierarchy
        for pquery in self.clustered_pqueries:  # type: phylo_seq.PQuery
            taxon_name = pquery.lct.split(hierarchy.lin_sep)[-1]
            pquery_taxon_map[pquery.place_name] = hierarchy.get_taxon(taxon_name)
            if not pquery_taxon_map[pquery.place_name]:
                pquery_taxon_map[pquery.place_name] = hierarchy.find_root_taxon()
        for leaf_node in tree.get_leaves():
            leaf_node.taxon = pquery_taxon_map[leaf_node.name]

        self.ref_pkg.propagate_internal_node_taxon_labels_from_leaves(tree)
        return

    def define_tree_clusters(self, tree: Tree, internal=True) -> None:
        if internal:
            # Split the cluster nodes into internal nodes when necessary
            cluster_nodes = self.split_node_partition_edges(self.partition_nodes(tree))
        else:
            cluster_nodes = self.partition_nodes(tree)

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

    def assign_pqueries_to_leaf_clusters(self, pqueries: list, pairwise_clusters: dict) -> None:
        # Build a map between the leaf node names and cluster numbers
        leaf_cluster_map = {}
        for cluster_num, p_otu in self.cluster_index.items():  # type: (int, PhylOTU)
            for leaf in p_otu.tree_node.get_leaf_names():
                leaf_cluster_map[leaf] = cluster_num

        # Match the PQuery names to clusters
        for pquery in pqueries:  # type: phylo_seq.PQuery
            try:
                pquery.p_otu = leaf_cluster_map[pquery.place_name]
            except KeyError:
                for _, cluster in pairwise_clusters.items():  # type: ts_classy.Cluster
                    if pquery.place_name in set([member[0] for member in cluster.members]):
                        pquery.p_otu = leaf_cluster_map[cluster.representative]
            if getattr(pquery, "sample_name") not in self.sample_mat:
                self.sample_mat[getattr(pquery, "sample_name")] = {n: 0 for n in self.cluster_index}
            try:
                self.sample_mat[getattr(pquery, "sample_name")][pquery.p_otu] += 1
            except KeyError:
                self.sample_mat[getattr(pquery, "sample_name")][pquery.p_otu] = 1
            self.cluster_index[pquery.p_otu].cardinality += 1
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

    @staticmethod
    def open_output_file(file_path: str):
        try:
            f_handler = open(file_path, 'w')
        except IOError:
            logging.error("Unable to open output file '{}' for writing.\n".format(file_path))
            sys.exit(13)
        return f_handler

    def write_cluster_taxon_labels(self) -> None:
        """Writes a two-column table mapping each cluster to its respective taxonomic lineage."""
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
        tbl_str = sep.join(["PQuery", "RefPkg", "Resolution", "Mode", "OTU_ID"]) + "\n"
        for pquery in sorted(self.clustered_pqueries, key=lambda x: x.seq_name):
            tbl_str += sep.join([pquery.seq_name,
                                 pquery.ref_name,
                                 self.tax_rank,
                                 self.clustering_mode,
                                 str(pquery.p_otu)]) + "\n"

        pquery_otu_tbl.write(tbl_str)
        pquery_otu_tbl.close()

        return

    def report_cluster_occupancy(self, sample_name, n_pqueries: int):
        occupancy = 0
        for _, p_otu in self.cluster_index.items():
            if p_otu.cardinality >= 1:
                occupancy += 1
        logging.info("{} classified sequences occupy {}/{} '{}' pOTUs in sample '{}'.\n"
                     "".format(n_pqueries, occupancy, len(self.cluster_index),
                               self.ref_pkg.prefix, sample_name))
        return


def dereplicate_by_clustering(fasta_inst: fasta.FASTA, prop_similarity: float, mmseqs_exe: str, tmp_dir: str,
                              subset=None, num_threads=2) -> dict:
    """
    A method for dereplicating a FASTA instance using pairwise sequence clustering with MMSeqs2.
    The FASTA instance, fasta_inst, fasta_dict attribute is modified to only store the representative (i.e. centroid)
    sequences determined from the sequence clustering.

    :param fasta_inst: A FASTA instance with the fasta_dict and header_registry loaded
    :param prop_similarity: The proportional similarity to cluster the sequences in fasta_inst
    :param mmseqs_exe: The path to a MMSeqs2 executable
    :param tmp_dir: A directory to write temporary files
    :param subset: Optionally, a list of sequences to cluster. Those not included will be removed from fasta_inst
    :param num_threads: The number of threads for MMSeqs2 to use (2 by default)
    :return: A dictionary of cluster numerical identifiers indexing Cluster instances
    """

    fasta_inst.change_dict_keys("num")
    cluster_input = os.path.join(tmp_dir, "cluster_in.fasta")
    clusters_prefix = os.path.join(tmp_dir, "linclust_out")
    clusters_table = clusters_prefix + "_cluster.tsv"
    cluster_alignments = clusters_prefix + "_cluster_aln.tsv"

    # Write a FASTA for clustering containing the formatted headers since
    # not all clustering tools + versions keep whole header - spaces are replaced with underscores
    fasta.write_new_fasta(fasta_dict=fasta_inst.fasta_dict, fasta_name=cluster_input, headers=subset)

    wrapper.cluster_sequences(software_path=mmseqs_exe,
                              fasta_input=cluster_input, output_prefix=clusters_prefix,
                              similarity=prop_similarity, num_threads=num_threads)

    cluster_map = file_parsers.create_mmseqs_clusters(clusters_table, cluster_alignments)

    # Revert headers in cluster_dict from 'formatted' back to 'original'
    fasta.rename_cluster_headers(cluster_map, fasta_inst.header_registry)
    logging.debug("\t{} sequence clusters\n".format(len(cluster_map.keys())))

    # Keep only the representative sequences in the FASTA instance
    fasta_inst.change_dict_keys()
    fasta_inst.keep_only(header_subset=[c.representative for num, c in cluster_map.items()])

    # Clean up the temporary files
    for tmp_file in [cluster_input, clusters_table, cluster_alignments]:
        os.remove(tmp_file)
    return cluster_map


def de_novo_phylogeny_from_queries(phylo_clust: PhyloClust, fasta_inst: fasta.FASTA) -> Tree:
    # TODO: Build a profile HMM from the incomplete query sequences and filter by length/coverage

    # Change header format to one without whitespace
    fasta_inst.change_dict_keys("num")
    # Write a FASTA file for multiple sequence alignment and phylogenetic inference
    fasta.write_new_fasta(fasta_inst.fasta_dict,
                          fasta_name=phylo_clust.stage_output_dir + "de_novo_clusters.fasta")
    wrapper.run_mafft(mafft_exe=phylo_clust.executables["mafft"],
                      fasta_in=phylo_clust.stage_output_dir + "de_novo_clusters.fasta",
                      fasta_out=phylo_clust.stage_output_dir + "de_novo_clusters.mfa",
                      num_threads=2)
    clusters_tree = wrapper.construct_tree(tree_builder="FastTree",
                                           executables=phylo_clust.executables,
                                           evo_model="LG",
                                           multiple_alignment_file=phylo_clust.stage_output_dir + "de_novo_clusters.mfa",
                                           tree_output_dir=phylo_clust.stage_output_dir,
                                           tree_prefix="de_novo_cluster_tree")

    # Ensure that the phylogeny is bifurcating
    ete_tree = Tree(clusters_tree)
    entish.label_internal_nodes_ete(ete_tree)
    # Propagate a 'taxon' feature - None by default - to all TreeNodes for holding Taxon instances
    for n in ete_tree.traverse(strategy="postorder"):  # type: Tree
        n.add_feature(pr_name="taxon", pr_value=None)
    # Add the query sequence name for each leaf
    for leaf_node in ete_tree.get_leaves():
        leaf_node.name = fasta_inst.header_registry[leaf_node.name].original
    return ete_tree


def pqueries_to_fasta(pqueries, fa_name="") -> fasta.FASTA:
    pqueries_fasta = fasta.FASTA(file_name=fa_name)
    for _, refpkg_pquery_map in pqueries.items():
        for _, pqueries in refpkg_pquery_map.items():  # type: (str, list)
            for pquery in pqueries:  # type: phylo_seq.PQuery
                pqueries_fasta.fasta_dict[pquery.place_name] = pquery.seq
    pqueries_fasta.header_registry = fasta.register_headers(list(pqueries_fasta.fasta_dict.keys()))
    return pqueries_fasta


def de_novo_phylo_clusters(phylo_clust: PhyloClust, drep_similarity=1.0):
    # Gather all classified sequences for a reference package from the treesapp assign outputs
    classified_pqueries = {}
    for sample_id, ts_out in phylo_clust.assign_output_dirs.items():
        classified_pqueries[sample_id] = file_parsers.load_classified_sequences_from_assign_output(ts_out,
                                                                                                   phylo_clust.ref_pkg.prefix)
        for pquery in classified_pqueries[sample_id][phylo_clust.ref_pkg.prefix]:
            phylo_clust.set_pquery_sample_name(pquery, sample_id)
            phylo_clust.clustered_pqueries.append(pquery)

    # Load the classified PQuery sequences into a FASTA instance
    pqueries_fasta = pqueries_to_fasta(classified_pqueries,
                                       fa_name=os.path.join(phylo_clust.stage_output_dir, "classified_pqueries.fasta"))
    pqueries_fasta.summarize_fasta_sequences()
    phylo_clust.increment_stage_dir()

    # Dereplicate classified sequences into Clusters
    if drep_similarity < 1.0:
        cluster_map = dereplicate_by_clustering(fasta_inst=pqueries_fasta,
                                                prop_similarity=drep_similarity,
                                                mmseqs_exe=phylo_clust.executables["mmseqs"],
                                                tmp_dir=phylo_clust.stage_output_dir)
    else:
        cluster_map = {}

    de_novo_phylogeny = de_novo_phylogeny_from_queries(phylo_clust, pqueries_fasta)
    # Use the taxonomy of each classified sequence as leaf node taxon attribute
    phylo_clust.set_pquery_leaf_node_taxa(de_novo_phylogeny)
    phylo_clust.increment_stage_dir()

    # Define tree clusters
    phylo_clust.define_tree_clusters(de_novo_phylogeny, internal=False)
    for sample_id in classified_pqueries:
        if phylo_clust.ref_pkg.prefix not in classified_pqueries[sample_id]:
            continue
        pqueries = classified_pqueries[sample_id][phylo_clust.ref_pkg.prefix]
        # Use cluster_map to assign pqueries to clusters that are not in the phylogeny
        phylo_clust.assign_pqueries_to_leaf_clusters(pqueries, cluster_map)

        phylo_clust.report_cluster_occupancy(sample_name=sample_id, n_pqueries=len(pqueries))
    return


def reference_placement_phylo_clusters(phylo_clust: PhyloClust, taxon_labelled_tree: Tree) -> None:
    # Perform maximum, mean, or median distance min-cut partitioning
    logging.info("Defining phylogenetic clusters... ")
    phylo_clust.define_tree_clusters(phylo_clust.ref_pkg.taxonomically_label_tree())
    logging.info("done.\n")

    # Find the edges belonging to each cluster and invert the dictionary to create the _edges_to_cluster_index
    phylo_clust.match_edges_to_clusters(tree=phylo_clust.ref_pkg.taxonomically_label_tree())

    logging.info("Assigning query sequences to clusters...\n")
    for sample_id, jplace_file in phylo_clust.jplace_files.items():
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
    return


def cluster_phylogeny(sys_args: list) -> None:
    p_clust = PhyloClust()
    p_clust.load_args(p_clust.get_arguments(sys_args))

    # Calculate RED distances for each node
    taxa_tree = p_clust.ref_pkg.taxonomically_label_tree()
    red_tree = rel_evo_dist.RedTree()
    red_tree.decorate_rel_dist(taxa_tree)

    if not p_clust.alpha:
        logging.info("Alpha threshold not provided. Estimating for '{}'... ".format(p_clust.tax_rank))
        p_clust.calculate_distance_threshold(taxa_tree, p_clust.ref_pkg.taxa_trie)
        logging.info("done.\n")
        logging.debug("Alpha estimated to be {}.\n".format(p_clust.alpha))

    if p_clust.clustering_mode == "ref_guided":
        reference_placement_phylo_clusters(phylo_clust=p_clust, taxon_labelled_tree=taxa_tree)
    elif p_clust.clustering_mode == "de_novo":
        de_novo_phylo_clusters(phylo_clust=p_clust, drep_similarity=0.999)
    else:
        logging.error("Unknown clustering mode specified: '{}'.\n".format(p_clust.clustering_mode))
        sys.exit(5)

    # Write a OTU table with the abundance of each PhylOTU for each JPlace
    p_clust.write_otu_matrix()

    # Write a table mapping each OTU identifier to its taxonomy
    p_clust.write_cluster_taxon_labels()

    # Write a table mapping PQuery sequence names to their cluster
    p_clust.write_pquery_otu_classifications()

    # TODO: Return pqueries with a 'cluster' attribute

    # TODO: Decide what to do with the PQueries with pendant lengths greater than the cluster RED radius

    # TODO: Centroids? Options are discussed in TreeCluster
    return


if __name__ == "__main__":
    cluster_phylogeny(sys.argv[1:])

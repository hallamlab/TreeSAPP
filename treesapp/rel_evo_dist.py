import logging

import collections

from ete3 import Tree, TreeNode

##
# The following code was adapted from https://github.com/dparks1134/PhyloRank/rel_dist.py
##


class RelativeEvolutionaryDistance:
    def __init__(self):
        """ Initialize """
        self.logger = logging.getLogger()
        return

    @staticmethod
    def _avg_descendant_rate(ete_tree) -> None:
        """
        Calculate average rate of divergence for each nodes in a tree.
        The average rate is the arithmetic mean of the branch length to all descendant taxa.

        :param ete_tree: A phylogenetic ete3 Tree instance
        :return: None
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_taxa: number of terminal taxa
        """

        # calculate the mean branch length to extant taxa
        for node in ete_tree.traverse(strategy="postorder"):  # type: TreeNode
            avg_div = 0
            if node.is_leaf():
                node.mean_dist = 0.0
                node.num_taxa = 1
            else:
                node.num_taxa = sum([1 for _ in node.get_leaves()])
                for c in node.iter_descendants():
                    avg_div += (float(c.num_taxa) / node.num_taxa) * (c.mean_dist + c.dist)

            node.mean_dist = avg_div
        return

    def decorate_rel_dist(self, ete_tree: Tree):
        """
        Calculate relative distance to each internal node.

        :param ete_tree: A phylogenetic ete3 Tree instance
        :return: None
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_taxa: number of terminal taxa
          rel_dists: relative distance of node between root and extant organisms
        """

        self._avg_descendant_rate(ete_tree)
        for node in ete_tree.traverse("preorder"):  # type: TreeNode
            if node.is_root():
                node.rel_dist = 0.0
            elif node.is_leaf():
                node.rel_dist = 1.0
            else:
                a = node.dist
                b = node.mean_dist
                x = node.up.rel_dist

                if (a + b) != 0:
                    rel_dist = x + (a / (a + b)) * (1.0 - x)
                else:
                    # Internal node has zero length to parent,
                    # so should have the same relative distance as the parent node
                    rel_dist = x

                node.rel_dist = rel_dist
        return

    def rel_dist_to_named_clades(self, tree) -> dict:
        """
        Determine relative distance to specific taxa.

        :param tree: A phylogenetic ete3 Tree instance
        :return: Dictionary of form d[rank_index][taxon] -> relative divergence
        """

        # calculate relative distance for all nodes
        self.decorate_rel_dist(tree)

        # tabulate values for internal nodes with ranks
        rel_dists = collections.defaultdict(dict)
        for node in tree.preorder_node_iter(lambda n: n != tree.seed_node):
            if not node.label or node.is_leaf():
                continue

            _support, taxon_name = parse_label(node.label)
            if not taxon_name:
                continue

            # get most-specific rank if a node represents multiple ranks
            if ';' in taxon_name:
                taxon_name = taxon_name.split(';')[-1].strip()

            most_specific_rank = taxon_name[0:3]
            rel_dists[Taxonomy.rank_index[most_specific_rank]][taxon_name] = node.rel_dist

        return rel_dists

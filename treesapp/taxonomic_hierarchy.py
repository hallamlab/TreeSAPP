__author__ = 'Connor Morgan-Lang'

import logging
import sys
import re
from pygtrie import StringTrie


class Taxon:
    def __init__(self, name, rank):
        self.name = name
        self.rank = rank
        self.parent = None
        self.prefix = ""
        self.taxid = ""
        self.coverage = 1

    def prefix_taxon(self, taxon_sep="__"):
        return self.prefix + taxon_sep + self.name

    def lineage(self) -> list:
        """
        Recursively retrieve the Taxon instances in a Taxon's lineage of parent instances

        :return: List of Taxon instances in order from Root to Species
        """
        if not self.parent:
            return [self]
        lineage = self.parent.lineage()
        lineage.append(self)
        return lineage


class TaxonomicHierarchy:
    def __init__(self, sep="; "):
        # Private variables - should not be changed after instantiation
        self.rank_name_map = {'superkingdom': 'domain', "strain": "type_strain"}
        self.accepted_ranks_depths = {"root": 0, "domain": 1, "phylum": 2, "class": 3, "order": 4,
                                      "family": 5, "genus": 6, "species": 7}
        self.no_rank_name = "no rank"  # Name of rank for non-canonical taxonomic ranks
        self.lin_sep = sep
        self.taxon_sep = "__"  # Separator between the rank prefix and the taxon name
        self.bad_taxa = []  # Optional list that can be used to ensure some taxa are removed (e.g. r__Root)
        self.canonical_prefix = re.compile(r"^[nrdpcofgs]" + re.escape(self.taxon_sep))
        self.proper_species_re = re.compile("^s__[A-Z][a-z]+ [a-z]+$")
        self.no_rank = re.compile(r"^" + re.escape(self.no_rank_name[0] + self.taxon_sep) + r".*")
        # Main data structures
        self.trie = StringTrie(separator=sep)  # Trie used for more efficient searches of taxa and whole lineages
        self.hierarchy = dict()  # Dict of prefix_taxon (e.g. p_Proteobacteria) name to Taxon instances
        self.rank_prefix_map = {self.no_rank_name[0]: {self.no_rank_name}}  # Dict to track prefixes representing ranks
        # The following are used for tracking the state of the instance's data structures
        self.include_prefix = True  # Keeps track of the trie's prefix for automated updates
        self.clean_trie = False  # To track whether taxa of "no rank" are included in the trie
        self.rank_prefix_map_values = set
        self.lineages_fed = 0
        self.lineages_into_trie = 0

    def get_state(self) -> dict:
        instance_state = {"clean_trie": self.clean_trie,
                          "include_prefix": self.include_prefix,
                          "rank_prefix_map_values": str(self.rank_prefix_map_values),
                          "accepted_ranks_depths": self.accepted_ranks_depths,
                          "lineages_fed": str(self.lineages_fed),
                          "lineages_into_tree": str(self.lineages_into_trie),
                          "rank_prefix_map": str(self.rank_prefix_map),
                          "taxon_sep": "'" + self.taxon_sep + "'",
                          "lin_sep": "'" + self.lin_sep + "'"}
        return instance_state

    def so_long_and_thanks_for_all_the_fish(self):
        summary_string = "Summary of TaxonomicHierarchy instance state:\n"
        for key, value in self.get_state().items():
            summary_string += "\t" + key + " = " + str(value) + "\n"
        logging.error(summary_string)
        sys.exit(17)

    def get_taxon_names(self, with_prefix=False):
        if with_prefix:
            return set([self.hierarchy[taxon].prefix_taxon() for taxon in self.hierarchy])
        else:
            return set([self.hierarchy[taxon].name for taxon in self.hierarchy])

    def get_taxon(self, prefix_taxon):
        """
        Retrieve the Taxon instance mapped to prefix_taxon in self.hierarchy.
        If prefix_taxon isn't in self.hierarchy, a warning is issued and None returned

        :param prefix_taxon: A string formatted such that the first letter of the taxonomic rank is a prefix, separated
         from the taxon's name by two underscores, like 'd__Bacteria'
        :return: Taxon instance mapped to the prefix_taxon in self.hierarchy
        """
        try:
            return self.hierarchy[prefix_taxon]  # type: Taxon
        except KeyError:
            logging.warning("Taxon name '{}' not present in taxonomic hierarchy.\n".format(prefix_taxon))
            return

    def resolved_as(self, lineage, rank_name="species") -> bool:
        """
        This determines whether a lineage is resolved to at least the rank provided,
        within the existing TaxonomicHierarchy instance.

        :param lineage: A taxonomic lineage *with prefixed ranks* separated by self.lin_sep
        :param rank_name: Name of the rank being compared to
        :return: Boolean indicating whether a lineage is as resolved as a rank
        """
        if self.rank_prefix_map_values is set:
            self.validate_rank_prefixes()

        if rank_name not in self.accepted_ranks_depths:
            logging.error("Rank name '{0}' is not valid in taxonomic hierarchy.\n")
            self.so_long_and_thanks_for_all_the_fish()

        taxon = lineage.split(self.lin_sep)[-1]
        taxon_rank = self.rank_prefix_map[taxon[0]]
        if self.accepted_ranks_depths[taxon_rank] >= self.accepted_ranks_depths[rank_name]:
            return True
        else:
            return False

    def validate_rank_prefixes(self):
        rank_prefix_map = {}
        if self.rank_prefix_map_values is set:  # Work here is already done
            for prefix in self.rank_prefix_map:
                rank_names = self.rank_prefix_map[prefix]  # type: set
                if len(rank_names) == 1:
                    rank_prefix_map[prefix] = rank_names.pop()
                elif len(rank_names) > 1:
                    logging.error("Conflicting rank names detected for rank prefix '{0}': "
                                  "'{1}'\n".format(prefix, ','.join(rank_names)))
                    sys.exit(5)
                else:
                    logging.warning("Prefix exists for missing rank name...? This isn't going to be good.\n")
            self.rank_prefix_map = rank_prefix_map
            self.rank_prefix_map_values = str
        return

    def trie_check(self):
        """
        Ensure the trie has contemporary with all lineages that have been fed into the hierarchy.
        This is accomplished by comparing the lineages_fed and lineages_into_trie properties.
        If they differ, trie rebuilds.

        :return: None
        """
        if self.lineages_fed != self.lineages_into_trie:
            self.build_multifurcating_trie(with_prefix=self.include_prefix)
        return

    def feed(self, lineage: str, lineage_ex: list) -> Taxon:
        """
        Takes lineage string, where each taxon is separated by self.lin_sep and, guided by the lineage_ex
        which is a list of extra lineage information, it adds each taxon to the dictionary self.hierarchy.
        Each taxon's name and rank is used to either find an existing Taxon object in self.hierarchy or create new ones.
        The lineage string is used to ensure the order of the taxon-parent relationships (linked-list) are correct.

        :param lineage: A taxonomic lineage string. Example:
         "cellular organisms; Bacteria"
        :param lineage_ex: A dictionary with taxon names as keys and taxonomic rank as values. Example:
         [{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'},
                          {'TaxId': '2', 'ScientificName': 'Bacteria', 'Rank': 'superkingdom'}]
        :return: None
        """
        previous = None

        # Ensure the rank_prefix_map values are sets that can be added to, in case self.validate_rank_prefixes was run
        if self.rank_prefix_map_values is not set:
            self.rank_prefix_map = {k: {v} for k, v in self.rank_prefix_map.items()}
            self.rank_prefix_map_values = set

        lineage = self.clean_lineage_string(lineage)
        taxa = lineage.split(self.lin_sep)
        while taxa and lineage_ex:
            taxon_info = lineage_ex.pop(0)
            taxon = taxa.pop(0)
            if taxon != taxon_info["ScientificName"]:
                logging.error("Lineage and Entrez lineage extra information list don't match. Current state:\n"
                              "taxon = {0}\n"
                              "taxon_info = {1}\n".format(taxon, str(taxon_info)))
                sys.exit(11)
            rank = taxon_info["Rank"]
            # Decide what the prefix will be
            if rank in self.accepted_ranks_depths:
                pass
            elif rank in self.rank_name_map:
                rank = self.rank_name_map[rank]
            else:
                rank = self.no_rank_name
            rank_prefix = rank[0]
            # Record the rank's prefix in rank_prefix_map
            try:
                self.rank_prefix_map[rank_prefix].add(rank)
            except KeyError:
                self.rank_prefix_map[rank_prefix] = {rank}

            prefix_name = rank_prefix + self.taxon_sep + taxon
            if prefix_name not in self.hierarchy:
                ti = Taxon(taxon, rank)
                ti.parent = previous
                ti.prefix = rank_prefix
                self.hierarchy[prefix_name] = ti
                previous = ti
            else:
                self.hierarchy[prefix_name].coverage += 1
                previous = self.hierarchy[prefix_name]

        if len(taxa) > 0 or len(lineage_ex) > 0:
            logging.error("Not all elements popped from paired lineage and lineage_ex information.\n"
                          "lineage list = {0}\n"
                          "lineage information dictionary = {1}\n".format(taxa, lineage_ex))

        # Update the number of lineages provided to TaxonomicHierarchy
        self.lineages_fed += 1

        return previous

    def emit(self, prefix_taxon: str, with_prefix=False) -> str:
        """
        Taking a prefixed-taxon (e.g. g__Actinomyces) as input it
        retrieves the corresponding Taxon instance from self.hierarchy.
        The parents of the taxon are traversed while a parent exists to recreate the original taxonomic lineage.

        :param prefix_taxon: A string formatted such that the first letter of the taxonomic rank is a prefix, separated
         from the taxon's name by two underscores, like 'd__Bacteria'
        :param with_prefix: Boolean controlling whether each taxon's prefix is included in the lineage string
        :return: A string of the taxon's lineage where each rank is separated by self.lin_sep
        """
        taxon = self.get_taxon(prefix_taxon)
        if not taxon:
            return ""

        anno_lineage = taxon.lineage()
        if with_prefix:
            return self.lin_sep.join([taxon.prefix_taxon() for taxon in anno_lineage])
        else:
            return self.lin_sep.join([taxon.name for taxon in anno_lineage])

    def build_multifurcating_trie(self, with_prefix=None) -> None:
        """
        The initial TaxonomicHierarchy is just made from a dictionary of String (e.g. d__Bacteria): Taxon pairs.
        After finding the correct key, a lineage can be formed by traversing upwards through the linked-list of 'parent's.

        However, finding the Keys across a set of taxa names is inefficient coming in at O(n^2).
        This function builds a prefix trie that could be leveraged for more efficiently determining whether
        all taxa in a lineage are present.
        Using self.trie (pygtrie.StringTrie object), it is populated by the lineage of each taxon in self.hierarchy.
        Nodes in the tree follow the format (lineage: taxon) where taxon lacks its rank-prefix. An example:

        StringTrie(n__cellular organisms: cellular organisms,
        n__cellular organisms; d__Bacteria: Bacteria,
        n__cellular organisms; d__Bacteria; n__Terrabacteria group: Terrabacteria group,
        n__cellular organisms; d__Bacteria; n__Terrabacteria group; p__Actinobacteria: Actinobacteria, separator=; )

        :type with_prefix: bool
        :param with_prefix: Flag indicating whether the lineages (node keys) have their rank-prefix
        :return: None
        """
        lineages = set()

        if with_prefix is not None:
            # Set the include_prefix boolean for future automated updates
            self.include_prefix = with_prefix

        for taxon_name in self.hierarchy:  # type: str
            lineages.add(self.emit(taxon_name, True))
        for lin in lineages:
            if self.clean_trie:
                lin = self.clean_lineage_string(lin, self.include_prefix)
            self.trie[lin] = self.canonical_prefix.sub('', lin.split("; ")[-1])

        self.lineages_into_trie = self.lineages_fed
        return

    def project_lineage(self, lineage_str):
        """
        Function for checking whether all taxa in a lineage exist within TaxonomicHierarchy

        :param lineage_str: A lineage string where taxa are separated by self.lin_sep. Can be prefixed (e.g. d__Archaea)
         or not, but it must match the state of self.trie. This is not checked!
        :return: Boolean indicating whether the entire lineage is present in the multifurcating taxonomic trie
        """
        self.trie_check()
        return lineage_str in self.trie

    def clean_lineage_string(self, lineage: str, with_prefix=None) -> str:
        """
        Removes superfluous taxonomic ranks and characters that make lineage comparisons difficult

        :type with_prefix: bool
        :param with_prefix: Flag indicating whether the lineages have their rank-prefix
        :param lineage: A taxonomic lineage string where each rank is separated by a semi-colon
        :return: String with the purified taxonomic lineage adhering to the NCBI hierarchy
        """
        reconstructed_lineage = []
        if with_prefix is not None:
            # Set the include_prefix boolean for future automated updates
            self.include_prefix = with_prefix

        if self.bad_taxa:
            # Potential list of bad strings: "cellular organisms; ", "delta/epsilon subdivisions; ", "\(miscellaneous\)"
            for bs in self.bad_taxa:
                lineage = re.sub(bs, '', lineage)

        ranks = lineage.split(self.lin_sep)
        for rank in ranks:
            if not self.no_rank.match(rank):
                if not self.include_prefix:
                    rank = self.canonical_prefix.sub('', rank)
                reconstructed_lineage.append(str(rank))
        return self.lin_sep.join(reconstructed_lineage)

    def check_lineage(self, lineage: str, organism_name: str, verbosity=0) -> str:
        """
        Sometimes the NCBI lineage is incomplete or the rank prefixes are out of order.
        This function checks (and fixes) the following things:
        1. Uses organism_name to add Species to the lineage if it is missing
        2. Ensure the progression of rank (i.e. from domain to species) is ordered properly

        :param lineage: A taxonomic lineage *with prefixed ranks* separated by self.lin_sep
        :param organism_name: Name of the organism. Parsed from the sequence header (usually at the end in square brackets)
        :param verbosity: 1 prints debugging messages
        :return: A list of elements for each taxonomic rank representing the taxonomic lineage
        """
        if not self.include_prefix or not self.clean_trie:
            self.clean_trie = True
            logging.debug("Switching multifurcating trie to include rank prefixes.\n")
            self.build_multifurcating_trie(True)

        if verbosity:
            logging.debug("check_lineage():\n\t"
                          "lineage = '{0}'\n\t"
                          "organism = '{1}'\n\t"
                          "include_prefix = {2}\n\t"
                          "clean_trie = {3}\n".format(lineage, organism_name, self.include_prefix, self.clean_trie))

        if not lineage:
            return ""

        lineage = self.clean_lineage_string(lineage)

        if not self.project_lineage(lineage):
            logging.error("Lineage '{0}' not present in taxonomic hierarchy.\n".format(lineage))
            self.so_long_and_thanks_for_all_the_fish()

        # Handle prefix discrepancy between lineage and organism_name if organism_name doesn't have rank prefix
        if not self.canonical_prefix.search(organism_name):
            for child_lineage, taxon in self.trie.items(prefix=lineage):  # type: (str, str)
                if taxon == organism_name:
                    organism_name = child_lineage.split(self.lin_sep)[-1]

        lineage_list = lineage.split(self.lin_sep)
        if self.proper_species_re.match(lineage_list[-1]):
            if verbosity:
                logging.debug("check_lineage(): Perfect lineage.\n")
        elif 6 <= len(lineage_list) < 7 and self.proper_species_re.match(organism_name):
            if verbosity:
                logging.debug("check_lineage(): Organism name added to complete the lineage.\n")
            lineage_list.append(organism_name)
        else:
            if verbosity:
                logging.debug("check_lineage(): Truncated lineage.\n")

        self.validate_rank_prefixes()  # Ensure the rank-prefix map is formatted correctly
        # Ensure the order and progression of ranks is correct (no domain -> phylum -> species for example)
        i = 0
        for taxon in lineage_list:  # type: str
            if self.accepted_ranks_depths[self.rank_prefix_map[taxon[0]]] != i+1:
                logging.warning("Order of taxonomic ranks in cleaned lineage '{0}' is unexpected.\n"
                                "Lineage will be truncated to '{1}'.\n".format(self.lin_sep.join(lineage_list),
                                                                               self.lin_sep.join(lineage_list[0:i])))
                lineage_list = lineage_list[0:i]
                break
            i += 1

        return self.lin_sep.join(lineage_list)

    def remove_leaf_nodes(self, taxa):
        if type(taxa) is str:
            taxa = [taxa]
        for prefix_taxon_name in taxa:
            leaf_taxon = self.get_taxon(prefix_taxon=prefix_taxon_name)  # type: Taxon
            if leaf_taxon:
                for taxon in leaf_taxon.lineage():
                    # Decrease Taxon.coverage for every Taxon instance in the lineage
                    taxon.coverage -= 1
                    # Remove a Taxon from self.hierarchy dictionary if its coverage equals zero
                    if taxon.coverage == 0:
                        self.hierarchy.pop(taxon.prefix_taxon())
                        # Update self.lin_fed
                        self.lineages_fed -= 1

        # Update self.trie is the number of lineages fed is not equal to number of lineages in the trie
        self.trie_check()
        return

    def summarize_taxa(self):
        """
        Function for enumerating the representation of each taxonomic rank within a TaxonomicHierarchy instance

        :return: A formatted, human-readable string stating the number of unique taxa at each rank
        """
        taxonomic_summary_string = ""
        taxa_counts = dict()

        for taxon_name in self.hierarchy:
            taxon = self.hierarchy[taxon_name]  # type: Taxon
            try:
                taxa_counts[taxon.rank] += 1
            except KeyError:
                taxa_counts[taxon.rank] = 1

        taxonomic_summary_string += "Number of unique lineages:\n"
        for rank in sorted(self.accepted_ranks_depths, key=lambda x: self.accepted_ranks_depths[x]):
            try:
                counts = taxa_counts[rank]
            except KeyError:  # Pass if a rank wasn't found in the lineages
                continue

            buffer = " "
            while len(rank) + len(str(counts)) + len(buffer) < 12:
                buffer += ' '
            taxonomic_summary_string += "\t" + rank + buffer + str(counts) + "\n"

        return taxonomic_summary_string

    # TODO: Move these to tests
    def test_th_feed(self):
        print("Testing feed()")
        dup_lineage = "cellular organisms; Bacteria; Terrabacteria group; Actinobacteria; Actinobacteria; Actinomycetales; Actinomycetaceae; Actinomyces; Actinomyces nasicola"
        trunc_lineage = "cellular organisms; Bacteria; Terrabacteria group; Actinobacteria; Actinobacteria"
        diff_lineage = "cellular organisms; Bacteria; Terrabacteria group; Actinobacteria; Actinobacteria; Bifidobacteriales; Bifidobacteriaceae; Bifidobacterium"
        dup_lineage_ex = [{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'},
                          {'TaxId': '2', 'ScientificName': 'Bacteria', 'Rank': 'superkingdom'},
                          {'TaxId': '1783272', 'ScientificName': 'Terrabacteria group', 'Rank': 'no rank'},
                          {'TaxId': '201174', 'ScientificName': 'Actinobacteria', 'Rank': 'phylum'},
                          {'TaxId': '1760', 'ScientificName': 'Actinobacteria', 'Rank': 'class'},
                          {'TaxId': '2037', 'ScientificName': 'Actinomycetales', 'Rank': 'order'},
                          {'TaxId': '2049', 'ScientificName': 'Actinomycetaceae', 'Rank': 'family'},
                          {'TaxId': '1654', 'ScientificName': 'Actinomyces', 'Rank': 'genus'},
                          {'ScientificName': "Actinomyces nasicola", 'Rank': "species"}]
        trunc_lineage_ex = [{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'},
                            {'TaxId': '2', 'ScientificName': 'Bacteria', 'Rank': 'superkingdom'},
                            {'TaxId': '1783272', 'ScientificName': 'Terrabacteria group', 'Rank': 'no rank'},
                            {'TaxId': '201174', 'ScientificName': 'Actinobacteria', 'Rank': 'phylum'},
                            {'TaxId': '1760', 'ScientificName': 'Actinobacteria', 'Rank': 'class'}]
        diff_lineage_ex = [{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'},
                           {'TaxId': '2', 'ScientificName': 'Bacteria', 'Rank': 'superkingdom'},
                           {'TaxId': '1783272', 'ScientificName': 'Terrabacteria group', 'Rank': 'no rank'},
                           {'TaxId': '201174', 'ScientificName': 'Actinobacteria', 'Rank': 'phylum'},
                           {'TaxId': '1760', 'ScientificName': 'Actinobacteria', 'Rank': 'class'},
                           {'TaxId': '85004', 'ScientificName': 'Bifidobacteriales', 'Rank': 'order'},
                           {'TaxId': '31953', 'ScientificName': 'Bifidobacteriaceae', 'Rank': 'family'},
                           {'TaxId': '1678', 'ScientificName': 'Bifidobacterium', 'Rank': 'genus'}]
        self.feed(lineage=dup_lineage, lineage_ex=dup_lineage_ex)
        self.feed(lineage=trunc_lineage, lineage_ex=trunc_lineage_ex)
        self.feed(lineage=diff_lineage, lineage_ex=diff_lineage_ex)
        print(self.hierarchy)
        return

    def test_th_emit(self):
        print("Testing emit()")
        taxon = "g__Actinomyces"
        print(self.emit(taxon))
        print(self.emit("g__Bifidobacterium"))
        return

    def test_check_lineage(self):
        print("Testing check_lineage()")
        t1_lineage = "n__cellular organisms; d__Bacteria; n__Terrabacteria group; p__Actinobacteria; c__Actinobacteria"
        t1_organism = "Actinobacteria"
        t2_lineage = "d__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces"
        t2_organism = "s__Actinomyces nasicola"
        t3_lineage = "d__Archaea"
        t3_organism = "Archaea"
        t4_lineage = "d__Bacteria; p__Cyanobacteria; o__Synechococcales; f__Prochloraceae; g__Prochlorococcus"
        t4_organism = "s__Prochlorococcus marinus"
        r = self.check_lineage(lineage=t1_lineage, organism_name=t1_organism)
        print("Test one complete:", r)
        r = self.check_lineage(lineage=t2_lineage, organism_name=t2_organism)
        print("Test two complete:", r)
        self.feed(lineage="cellular organisms; Archaea; Euryarchaeota",
                  lineage_ex=[{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'},
                              {'ScientificName': 'Archaea', 'Rank': 'domain'},
                              {'ScientificName': 'Euryarchaeota', 'Rank': 'phylum'}])
        self.feed("Bacteria; Cyanobacteria; Synechococcales; Prochloraceae; Prochlorococcus",
                  [{'ScientificName': 'Bacteria', 'Rank': 'domain'},
                   {'ScientificName': 'Cyanobacteria', 'Rank': 'phylum'},
                   {'ScientificName': 'Synechococcales', 'Rank': 'order'},
                   {'ScientificName': 'Prochloraceae', 'Rank': 'family'},
                   {'ScientificName': 'Prochlorococcus', 'Rank': 'genus'}])
        r = self.check_lineage(lineage=t3_lineage, organism_name=t3_organism)
        print("Test three complete:", r)
        r = self.check_lineage(lineage=t4_lineage, organism_name=t4_organism)
        print("Test four complete:", r)

        return

    def test_clean_lineage_string(self):
        print("Testing clean_lineage_string()")
        l1 = "n__cellular organisms; d__Bacteria; n__Terrabacteria group"
        l2 = "d__Bacteria"
        self.include_prefix = False
        print(self.clean_lineage_string(l1) == self.clean_lineage_string(l2))
        self.include_prefix = True
        print(self.clean_lineage_string(l1) == self.clean_lineage_string(l2))
        return

    def test_remove_leaf_nodes(self):
        print("Testing remove_leaf_nodes()")
        self.remove_leaf_nodes("p__Euryarchaeota")


def main():
    th = TaxonomicHierarchy()
    th.test_clean_lineage_string()
    th.test_th_feed()
    th.test_th_emit()
    th.build_multifurcating_trie()
    th.test_check_lineage()
    print(th.summarize_taxa())
    th.test_remove_leaf_nodes()
    print(th.summarize_taxa())


if __name__ == "__main__":
    main()

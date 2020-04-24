__author__ = 'Connor Morgan-Lang'

import logging
import sys


class Taxon:
    def __init__(self, name, rank):
        self.name = name
        self.rank = rank
        self.parent = None
        self.taxid = ""


class TaxonomicHierarchy:
    def __init__(self, sep="; "):
        self.hierarchy = dict()  # A dictionary of prefix_taxon (e.g. p_Proteobacteria) name to Taxon instances
        self.lin_sep = sep
        self.taxon_sep = "__"
        self.no_rank_name = "no rank"
        self.rank_prefixes = dict()
        self.prefix_map = {'superkingdom': 'domain', "strain": "type_strain"}
        self.accepted_ranks = {"domain", "phylum", "class", "order", "family", "genus", "species"}

    def feed(self, lineage: str, lineage_ex: list) -> None:
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
            if rank in self.accepted_ranks:
                pass
            elif rank in self.prefix_map:
                rank = self.prefix_map[rank]
            else:
                rank = self.no_rank_name
            prefix_name = rank[0] + self.taxon_sep + taxon
            if prefix_name not in self.hierarchy:
                ti = Taxon(taxon, rank)
                ti.parent = previous
                self.hierarchy[prefix_name] = ti
                previous = ti
            else:
                previous = self.hierarchy[prefix_name]

        if len(taxa) > 0 or len(lineage_ex) > 0:
            logging.error("Not all elements popped from paired lineage and lineage_ex information.\n"
                          "lineage list = {0}\n"
                          "lineage information dictionary = {1}\n".format(taxa, lineage_ex))
        return

    def emit(self, prefix_taxon: str) -> str:
        """
        Taking a prefixed-taxon (e.g. g__Actinomyces) as input it
        retrieves the corresponding Taxon instance from self.hierarchy.
        The parents of the taxon are traversed while a parent exists to recreate the original taxonomic lineage.

        :param prefix_taxon: A string formatted such that the first letter of the taxonomic rank is a prefix, separated
         from the taxon's name by two underscores, like 'd__Bacteria'
        :return: A string of the taxon's lineage where each rank is separated by self.lin_sep
        """
        anno_lineage = []
        try:
            taxon = self.hierarchy[prefix_taxon]  # type: Taxon
        except KeyError:
            logging.warning("Taxon name '{}' not present in taxonomic hierarchy.\n".format(prefix_taxon))
            return ""

        anno_lineage.append(taxon.name)
        while taxon.parent:  # type: Taxon
            anno_lineage.append(taxon.parent.name)
            taxon = taxon.parent

        return self.lin_sep.join(reversed(anno_lineage))

    # TODO: Move these to tests
    def test_th_feed(self):
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
        taxon = "g__Actinomyces"
        print(self.emit(taxon))
        print(self.emit("g__Bifidobacterium"))
        return

import logging
import os
import sys

from treesapp import fasta
from treesapp import wrapper
from treesapp import logger

LOGGER = logging.getLogger(logger.logger_name())


class Cluster:
    def __init__(self, rep_name):
        self.representative = rep_name
        self.members = list()
        self.lca = ''

    def __str__(self):
        return "Cluster instance with {} members represented by '{}'".format(len(self.members), self.representative)

    def get_info(self):
        info_string = "Representative: " + str(self.representative) + "\n" + \
                      "LCA: " + self.lca + "\n" + \
                      "Members:\n\t" + "\n\t".join([', '.join(member) for member in self.members]) + "\n"
        return info_string


class BlastAln:
    def __init__(self):
        self.subject = ""
        self.query = ""
        self.start = 0
        self.end = 0
        self.alnlen = 0
        self.gapopen = 0
        self.mismatch = 0
        self.qstart = 0
        self.qend = 0
        self.tstart = 0
        self.tend = 0
        self.evalue = 0.0
        self.pident = 0.0
        self.bits = 0.0

    def load_blast_tab(self, line: str, sep="\t") -> None:
        try:
            fields = line.strip().split(sep)
        except ValueError:
            LOGGER.error("Unable to parse line in from alignment table:\n{}\n".format(line))
            sys.exit(7)

        try:
            self.subject, self.query = fields[0], fields[1]
            self.pident, self.alnlen, self.mismatch, self.gapopen = fields[2:6]
            self.qstart, self.qend, self.tstart, self.tend = fields[6:10]
            self.evalue, self.bits = fields[10], fields[11]
        except ValueError:
            LOGGER.error("Incorrect format for line in alignment file. Twelve were expected, found {}.\n{}\n"
                          "".format(len(fields), line))
            sys.exit(7)

        return


def dereplicate_by_clustering(fasta_inst: fasta.FASTA, prop_similarity: float, mmseqs_exe: str, tmp_dir: str,
                              subset=None, num_threads=2) -> dict:
    """
    A method for dereplicating a FASTA instance using pairwise sequence clustering with MMSeqs2.

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

    cluster_map = create_mmseqs_clusters(clusters_table, cluster_alignments)

    # Revert headers in cluster_dict from 'formatted' back to 'original'
    fasta.rename_cluster_headers(cluster_map, fasta_inst.header_registry)
    LOGGER.debug("\t{} sequence clusters\n".format(len(cluster_map.keys())))

    # Keep only the representative sequences in the FASTA instance
    fasta_inst.change_dict_keys()
    fasta_inst.keep_only(header_subset=[c.representative for num, c in cluster_map.items()])

    # Clean up the temporary files
    for tmp_file in [cluster_input, clusters_table, cluster_alignments]:
        os.remove(tmp_file)
    return cluster_map


def read_uc(uc_file: str) -> dict:
    """
    Function to read a VSEARCH cluster (.uc) file

    :param uc_file: Path to a .uc file produced by VSEARCH
    :return: Dictionary where keys are numerical identifiers and values are Cluster objects
        The Cluster object
    """
    cluster_dict = dict()
    try:
        uc = open(uc_file, 'r')
    except IOError:
        LOGGER.error("Unable to open VSEARCH cluster file " + uc_file + " for reading!\n")
        sys.exit(13)

    LOGGER.debug("Reading VSEARCH cluster file... ")

    # Find all clusters with multiple identical sequences
    for line in uc:
        cluster_type, num_id, _length, identity, _, _, _, _cigar, header, _representative = line.strip().split("\t")
        if cluster_type == "S":
            cluster_dict[num_id] = Cluster(header)
        elif cluster_type == "H":
            cluster_dict[num_id].members.append([header, identity])
        elif cluster_type == "C":
            pass
        else:
            LOGGER.error("Unexpected cluster type '" + str(cluster_type) + "' in " + uc_file + "\n")
            sys.exit(13)

    uc.close()
    LOGGER.debug("done.\n")
    return cluster_dict


def read_linclust_clusters(clusters_tsv_file: str) -> dict:
    """
    Reads the two-column TSV file representing MMSeqs clusters, where the first column is the name of the
    representative sequence and the second is the name of the member sequence.

    :param clusters_tsv_file: Path to the two-column TSV file written by `mmseqs createtsv`.
    :return: A dictionary mapping the cluster number (in order of iteration) mapped to its Cluster instance.
    """
    try:
        cluster_tbl = open(clusters_tsv_file)
    except IOError:
        LOGGER.error("Unable to open MMSeqs clusters table '{}' for reading!\n".format(clusters_tsv_file))
        sys.exit(13)

    previous = ""
    cluster_acc = 0
    clusters = dict()
    cl_inst = None
    for line in cluster_tbl:
        if not line:
            continue
        try:
            rep_name, mem_name = line.strip().split("\t")
        except ValueError:
            LOGGER.error("Unacceptable line format in MMSeqs cluster table:\n{}\n".format(line))
            sys.exit(17)
        if rep_name != previous:
            cl_inst = Cluster(rep_name)
            clusters[str(cluster_acc)] = cl_inst
            cluster_acc += 1
        cl_inst.members.append([mem_name, 0.0])
        previous = rep_name

    cluster_tbl.close()

    return clusters


def create_mmseqs_clusters(clusters_tbl: str, aln_tbl: str) -> dict:
    """
    Joins an cluster-specific alignment table created by MMSeqs.

    :param clusters_tbl: Path to a tab-separated values file with 2 columns, one each for the sequence names of
     the cluster representative and for the cluster member.
    :param aln_tbl: Path to a tab-separated values file with 12 columns:
     (1,2) identifiers for query and target sequences/profiles,
     (3) sequence identity,
     (4) alignment length,
     (5) number of mismatches,
     (6) number of gap openings,
     (7-8, 9-10) domain start and end-position in query and in target,
     (11) E-value, and (12) bit score.
    :return: Dictionary where keys are numerical identifiers representing the cluster ID and values are Cluster objects
    """
    # Define the clusters
    clusters = read_linclust_clusters(clusters_tbl)

    try:
        aln = open(aln_tbl, 'r')
    except IOError:
        LOGGER.error("Unable to open MMSeqs alignment table '{}' for reading!\n".format(aln_tbl))
        sys.exit(13)

    # Read the alignments, storing them in a dictionary indexed by subject/target/representative names
    aln_map = {}
    for line in aln:
        if not line:
            continue
        alignment = BlastAln()
        alignment.load_blast_tab(line)
        try:
            aln_map[alignment.subject].append(alignment)
        except KeyError:
            aln_map[alignment.subject] = [alignment]
    aln.close()

    for _, cluster in clusters.items():  # type: (int, Cluster)
        rep_alignments = aln_map[cluster.representative]
        for member in cluster.members:  # type: [str, float]
            x = 0
            while x < len(rep_alignments):
                alignment = rep_alignments[x]  # type: BlastAln
                if alignment.query == member[0]:
                    member[1] = alignment.pident
                    rep_alignments.pop(x)
                else:
                    x += 1
            if member[1] == 0:
                LOGGER.error("Unable to find an alignment between representative '{}' and query '{}' from MMSeqs.\n")
                sys.exit(11)

    return clusters

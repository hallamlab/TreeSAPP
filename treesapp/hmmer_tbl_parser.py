import sys
import re
import logging
from collections import namedtuple

from treesapp import logger

LOGGER = logging.getLogger(logger.logger_name())


def format_hmmer_domtbl_line(line):
    stats = []
    stat = ""
    for c in line:
        if c == ' ':
            if len(stat) > 0:
                stats.append(str(stat))
                stat = ""
            else:
                pass
        else:
            stat += c
    stats.append(str(stat))
    return stats


class HmmSearchStats:
    def __init__(self):
        self.raw_alignments = 0
        self.seqs_identified = 0
        self.fragmented = 0
        self.inverted = 0
        self.glued = 0
        self.dropped = 0
        self.bad = 0
        self.short = 0
        self.multi_alignments = 0  # matches of the same query to a different HMM (>1 lines)

    def num_dropped(self):
        self.dropped = self.inverted + self.bad + self.short
        return self.dropped

    def summarize(self):
        alignment_stat_string = ""
        alignment_stat_string += "\tInitial alignments:\t" + str(self.raw_alignments) + "\n"
        alignment_stat_string += "\tAlignments discarded:\t" + str(self.dropped) + "\n"
        alignment_stat_string += "\tFragmented alignments:\t" + str(self.fragmented) + "\n"
        alignment_stat_string += "\tInversions detected:\t" + str(self.inverted) + "\n"
        alignment_stat_string += "\tAlignments scaffolded:\t" + str(self.glued) + "\n"
        alignment_stat_string += "\tMulti-alignments:\t" + str(self.multi_alignments) + "\n"
        alignment_stat_string += "\tSequences identified:\t" + str(self.seqs_identified) + "\n"
        return alignment_stat_string


class HmmMatch:
    def __init__(self):
        self.genome = ""  # Name of the input file (Metagenome, SAG, MAG, or isolate genome)
        self.target_hmm = ""  # Name of the HMM aligned to
        self.orf = ""  # Name of the ORF, or more generally contig sequence
        self.hmm_len = 0  # Length of the hidden Markov model
        self.start = 0  # Alignment start position on the contig
        self.end = 0  # Alignment end position on the contig
        self.pstart = 0  # Alignment start position on the hmm profile
        self.pend = 0  # Alignment end position on the hmm profile
        self.seq_len = 0  # Length of the query sequence
        self.num = 0
        self.of = 0
        self.desc = ""
        self.acc = 0.0
        self.ieval = 0.0
        self.eval = 0.0  # Full-sequence E-value (in the case a sequence alignment is split)
        self.full_score = 0
        self.next_domain = None  # The next domain aligned by hmmsearch

    def get_info(self):
        info_string = "Info for query " + str(self.orf) + ":\n"
        info_string += "\tHMM = " + self.target_hmm + ", length = " + str(self.hmm_len) + "\n"
        info_string += "\tSequence length = " + str(self.seq_len) + "\n"
        info_string += "\tAligned length = " + str(self.end - self.start) + "\n"
        info_string += "\tAlignment start = " + str(self.start) + ", alignment stop = " + str(self.end) + "\n"
        info_string += "\tProfile start = " + str(self.pstart) + ", profile stop = " + str(self.pend) + "\n"
        info_string += "\tNumber " + str(self.num) + " of " + str(self.of) + "\n"
        info_string += "\tcE-value = " + str(self.ieval) + "\n"
        info_string += "\tacc = " + str(self.acc) + "\n"
        info_string += "\tfull score = " + str(self.full_score) + "\n"
        return info_string

    def subsequent_matches(self):
        if not self.next_domain:
            return [self]
        return [self] + self.next_domain.subsequent_matches()

    def collinear(self) -> bool:
        """
        This checks the profile-start and profile-stop positions against the other alignment fragments to ensure
        they are the smallest (derived from the left-most position of the HMM profile).

        :return: Boolean representing whether all domains/alignment fragments are collinear (True) or not (False)
        """
        if not self.next_domain:
            return True
        if self.next_domain.collinear():
            if self.pend < self.next_domain.pend and self.pstart < self.next_domain.pstart:
                return True
        return False

    def contains_duplicate_loci(self, overlap_threshold=0.33) -> bool:
        """
        Meant to identify whether the query sequence contains duplicate loci of the same HMM profile. If at least two
        HMM alignments overlap on the profile by more than the threshold (default 33%) they are deemed duplicate loci.

        :param overlap_threshold: A float representing the maximum HMM profile length proportion the two alignments
         can overlap on the HMM profile to be considered discrete loci; they are considered duplicates if exceeded.
        :return: Boolean, True indicates there are duplicate loci
        """
        if not self.next_domain:
            return False
        if self.next_domain.contains_duplicate_loci(overlap_threshold):
            return True
        for hmm_match in sorted(self.subsequent_matches(), key=lambda x: x.num):  # type: HmmMatch
            # Duplicate loci should not overlap on the query
            # Duplicate loci should overlap significantly on the HMM
            if self.num != hmm_match.num:
                # Check for redundant profile HMM coverage
                p_overlap_len = overlap_length(self.pstart, self.pend,
                                               hmm_match.pstart, hmm_match.pend)
                q_overlap_len = overlap_length(self.start, self.end,
                                               hmm_match.start, hmm_match.end)
                if q_overlap_len > 0 or p_overlap_len/self.hmm_len < overlap_threshold:
                    return False
                else:
                    return True
            else:
                pass
        return False

    def enumerate(self) -> None:
        i = 1
        next_matches = sorted(self.subsequent_matches(), key=lambda x: x.num)
        of = len(next_matches)
        for hmm_match in next_matches:  # type: HmmMatch
            hmm_match.num = i
            hmm_match.of = of
            i += 1
        return

    def copy(self, new_match):
        for i in self.__dict__:
            self.__dict__[i] = new_match.__dict__[i]
        return

    def drop_match_at(self, index: int) -> None:
        """
        Function for removing an element from the linked list at a position

        :param index: The index in the linked list that is to be dropped
        :return: None
        """
        if index < 0 or index >= len(self.subsequent_matches()):
            LOGGER.error("Unable to drop HmmMatch from next_domain linked list at index {}.\n"
                          "HmmMatch instance status:\n{}\n".format(index, self.get_info()))
            sys.exit(9)
        elif index == 0:
            self.copy(self.next_domain)
        else:
            try:
                self.subsequent_matches()[index-1].next_domain = self.subsequent_matches()[index+1]
            except IndexError:
                self.subsequent_matches()[index-1].next_domain = None

        return

    def merge_alignment_fragments(self, index: int) -> None:
        # Update the alignments that are being merged into the base/reference alignment
        merging_aln = self.subsequent_matches()[index]
        self.start = min([self.start, merging_aln.start])
        self.end = max([self.end, merging_aln.end])
        self.pstart = min([self.pstart, merging_aln.pstart])
        self.pend = max([self.pend, merging_aln.pend])
        self.ieval = min([self.eval, self.ieval, merging_aln.ieval])
        self.full_score = max([self.full_score, merging_aln.full_score])

        self.drop_match_at(index)
        return

    def scaffold_domain_alignments(self, seq_length_wobble=1.2) -> None:
        """
        If one or more alignments do not completely redundantly cover the HMM profile,
        overlap or are within a few BPs of each other of the query sequence,
        and do not generate an alignment 120% longer than the HMM profile, then
        merge the alignment co-ordinates, average the acc, Eval, cEval and make 'num' and 'of' reflect number of alignments
        Takes this:
        -------------
                    --------------
                                    ----------------            --------------------------------------
        and converts it to:
        ---------------------------------------------           --------------------------------------

        :param seq_length_wobble: The scalar threshold controlling the maximum size of a scaffold.
        Since we're dealing with HMMs, its difficult to estimate the sequence length variance so we're allowing for
        some 'wobble' in how long or short paralogs could be.
        :return: None
        """
        if not self.next_domain:
            return
        self.next_domain.scaffold_domain_alignments()

        if self.contains_inversion():
            return
        i = 1
        while i < len(self.subsequent_matches()):
            next_match = sorted(self.subsequent_matches(), key=lambda x: x.num)[i]  # type: HmmMatch
            if self.num == next_match.num:
                LOGGER.error("Iteration of HmmMatch.subsequent_matches() begins at current instance (%s).\n" % self.num)
                sys.exit(9)
            # Check for sub- or super-sequence orientation on the query sequence
            q_overlap_len = overlap_length(self.start, self.end, next_match.start, next_match.end)
            # Check for redundant profile HMM coverage
            p_overlap_len = overlap_length(self.pstart, self.pend, next_match.pstart, next_match.pend)

            min_profile_covered = min([self.pend - self.pstart,
                                       next_match.pend - next_match.pstart])
            try:
                aln_overlap_proportion = p_overlap_len / min_profile_covered
            except ZeroDivisionError:
                if self.pend - self.pstart < next_match.pend - next_match.pstart:
                    self.drop_match_at(i-1)
                else:
                    self.drop_match_at(i)
                continue
            if aln_overlap_proportion > 0.5:
                # They overlap significant regions - are they overlapping sequence or are they repeats?
                if q_overlap_len == next_match.seq_len:
                    # The projected alignment is a subsequence of the base alignment - DROP
                    self.drop_match_at(i)
                    i -= 1
                elif q_overlap_len < p_overlap_len:
                    # The two alignments represent repeats of a single profile - KEEP
                    pass
                else:
                    # The two alignments represent something...
                    self.merge_alignment_fragments(i)
                    i -= 1
            else:
                # They do not significantly overlap in the profile, likely part of the same alignment
                a_new_start = min([self.start, next_match.start])
                a_new_end = max([self.end, next_match.end])
                if float(a_new_end - a_new_start) < float(seq_length_wobble * int(self.hmm_len)):
                    # The proximal alignments overlap so they would probably form a homologous sequence - MERGE
                    self.merge_alignment_fragments(i)
                    i -= 1
                else:
                    # The alignments are very far apart from each other and would form a monster sequence - KEEP
                    pass
            i += 1
        return

    def contains_inversion(self) -> bool:
        """
        Asserting this HmmMatch instance is the first alignment (and therefore the left-most alignment on the query)
        this checks for collinearity along the HMM profile and query sequences.
        If they are not collinear an additional operation to look for duplicated HMM profile regions in the query,
        indicating multiple domains, is performed.
        Otherwise the positions have been inverted, or in other words, a section of the gene has been rearranged.

        :return: Failing both tests for collinearity and duplication return True, else False
        """
        if self.of == 1 or not self.next_domain:
            return False
        if self.collinear():
            return False
        elif self.contains_duplicate_loci(0.1):
            return False
        else:
            return True

    def orient_alignments(self) -> dict:
        """
        Creates a dictionary of alignment relation classes ('supersequence', 'subsequence', 'overlap' and 'satellite')
        indexed by the HmmMatch.num of a pair of alignments.

        :return: Dictionary
        """
        alignment_relations = dict()
        if not self.next_domain:
            return alignment_relations
        alignment_relations.update(self.next_domain.orient_alignments())

        for next_match in self.subsequent_matches():
            if self.num != next_match.num:
                alignment_relations[(self.num, next_match.num)] = detect_orientation(self.start, self.end,
                                                                                     next_match.start, next_match.end)
        return alignment_relations

    def consolidate_subalignments(self, alignment_relations: dict) -> dict:
        """
        What is the point of this function?
        Identify alignments that can be

        :param alignment_relations:
        :return: Dictionary of HmmMatch instances indexed by their respective query names
        """
        alignments_for_disposal = set()
        scaffolded_alignments = dict()
        for pair in alignment_relations:
            # If there are multiple alignments that span the whole hmm profile, report them both
            base, projected = pair
            if alignment_relations[pair] == "satellite" or alignment_relations[pair] == "overlap":
                pass
            elif alignment_relations[pair] == "supersequence":
                alignments_for_disposal.add(projected)
            elif alignment_relations[pair] == "subsequence":
                alignments_for_disposal.add(base)
            else:
                LOGGER.error("Unexpected alignment comparison state: '" + alignment_relations[pair] + "'\n")
                sys.exit(31)

        for hmm_match in self.subsequent_matches():
            # Filter out the worst of the overlapping alignments that couldn't be scaffolded
            if hmm_match.num not in alignments_for_disposal:
                query_header = ' '.join([hmm_match.orf, hmm_match.desc]) + \
                               '_' + str(hmm_match.num) + '_' + str(hmm_match.of)
                scaffolded_alignments[query_header] = hmm_match

        return scaffolded_alignments


class DomainTableParser(object):

    def __init__(self, dom_tbl):
        self.alignments = {}
        self.i = 0
        self.lines = []
        self.size = 0
        try:
            self.commentPattern = re.compile(r'^#')
            self.src = open(dom_tbl)
        except IOError:
            LOGGER.error("Could not open " + dom_tbl + " or file is not available for reading.\n")
            sys.exit(0)

    def read_domtbl_lines(self):
        """
        Function to read the lines in the domain table file,
        skipping those matching the comment pattern

        :return: self.lines is a list populated with the lines
        """
        line = self.src.readline()
        while line:
            comment = self.commentPattern.match(line)
            if not comment:
                self.lines.append(line.strip())
            if not line:
                break
            line = self.src.readline()
        self.size = len(self.lines)

    def next(self):
        """
        Reformat the raw lines of the domain table into
        an easily accessible hmm_domainTable format and perform
        QC to validate the significance of the alignments
        """
        if self.i < self.size:
            hit = format_hmmer_domtbl_line(self.lines[self.i])
            self.prepare_data(hit)
            self.i += 1

            try:
                return self.alignments
            except ValueError:
                return None
        else:
            self.src.close()
            return None

    def prepare_data(self, hit):
        self.alignments['query'] = str(hit[0])
        self.alignments['query_len'] = int(hit[2])
        self.alignments['hmm_name'] = str(hit[3])
        self.alignments['hmm_len'] = int(hit[5])
        self.alignments['Eval'] = float(hit[6])  # Full-sequence E-value (in the case a sequence alignment is split)
        self.alignments['full_score'] = float(hit[7])  # Full-sequence score
        self.alignments['num'] = int(hit[9])  # HMMER is able to detect whether there are multi-hits
        self.alignments['of'] = int(hit[10])  # This is the number of multi-hits for a query
        self.alignments['cEval'] = float(hit[11])  # conditional E-value
        self.alignments['iEval'] = float(hit[12])  # conditional E-value
        self.alignments['pstart'] = int(hit[15])  # First position on HMM profile
        self.alignments['pend'] = int(hit[16])  # Last position on HMM profile
        self.alignments['qstart'] = int(hit[19])  # env coord from
        self.alignments['qend'] = int(hit[20])  # env coord to
        self.alignments['acc'] = float(hit[21])
        self.alignments['desc'] = ' '.join(hit[22:])


def prep_args_for_parsing(args) -> namedtuple:
    """
    Check whether specific attributes used for filtering alignments exist in
    the args object created by Argparse.parse_args(), and add them if they do not.
    Create a namedtuple object with max_e, max_ie, min_acc, min_score and perc_aligned attributes and
    populate it with the filtering attributes in args.

    :param args: An object created by Argparse.parse_args()
    :return: A namedtuple with max_e, max_ie, min_acc, min_score and perc_aligned attributes
    """
    thresholds = namedtuple("thresholds", "max_e max_ie min_acc min_score perc_aligned profile_match")
    if not hasattr(args, "max_e"):
        args.max_e = 1E-5
    thresholds.max_e = args.max_e
    if not hasattr(args, "max_ie"):
        args.max_ie = 1E-3
    thresholds.max_ie = args.max_ie
    if not hasattr(args, "min_acc"):
        args.min_acc = 0.7
    thresholds.min_acc = args.min_acc
    if not hasattr(args, "min_score"):
        args.min_score = 20
    thresholds.min_score = args.min_score
    if not hasattr(args, "perc_aligned"):
        args.perc_aligned = 60
    thresholds.perc_aligned = args.perc_aligned
    if not hasattr(args, "query_aligned"):
        args.query_aligned = args.perc_aligned
    thresholds.query_aligned = args.query_aligned
    thresholds.profile_match = True

    # Print some stuff to inform the user what they're running and what thresholds are being used.
    info_string = "Filtering HMM alignments using the following thresholds:\n"
    info_string += "\tMaximum E-value = " + str(thresholds.max_e) + "\n"
    info_string += "\tMaximum i-Evalue = " + str(thresholds.max_ie) + "\n"
    info_string += "\tMinimum acc = " + str(thresholds.min_acc) + "\n"
    info_string += "\tMinimum score = " + str(thresholds.min_score) + "\n"
    info_string += "\tMinimum percentage of the HMM covered = " + str(thresholds.perc_aligned) + "%\n"
    info_string += "\tMinimum percentage of the query aligned = " + str(thresholds.query_aligned) + "%\n"
    LOGGER.debug(info_string)

    return thresholds


def detect_orientation(q_i: int, q_j: int, r_i: int, r_j: int) -> str:
    """
    Returns the class of orientation ('supersequence', 'subsequence', 'overlap' and 'satellite') based on
    the number of positions the query (base) alignment overlaps with the reference (projected) alignment.

    :param q_i: query start position
    :param q_j: query end position
    :param r_i: reference start position
    :param r_j: reference end position
    :return: String representing the class of orientation in relation to the query (q)
    """
    if q_i <= r_i <= q_j:
        if q_i <= r_j <= q_j:
            return "supersequence"
        else:
            return "overlap"
    elif r_i <= q_i <= r_j:
        if r_i <= q_j <= r_j:
            return "subsequence"
        else:
            return "overlap"
    else:
        return "satellite"


def overlap_length(r_i: int, r_j: int, q_i: int, q_j: int) -> int:
    """
    Returns the number of positions the query (base) alignment overlaps with the reference (projected) alignment

    :param q_i: query start position
    :param q_j: query end position
    :param r_i: reference start position
    :param r_j: reference end position
    :return: Number of positions the two alignments overlap
    """
    if r_j < q_i or q_j < r_i:
        # Satellite alignments
        return 0
    else:
        return min(r_j, q_j) - max(r_i, q_i)


def assemble_domain_alignments(first_match: HmmMatch, search_stats: HmmSearchStats):
    distinct_alignments = dict()
    if first_match.next_domain:
        # STEP 1: Scaffold the alignments covering overlapping regions on the query sequence
        frags_pre_scaffolding = len(first_match.subsequent_matches())
        first_match = remove_redundant_alignments(first_match)
        first_match.scaffold_domain_alignments()
        if first_match.contains_inversion():
            search_stats.inverted += 1
            return distinct_alignments
        first_match.enumerate()
        search_stats.glued += (frags_pre_scaffolding - len(first_match.subsequent_matches()))

        # STEP 2: Determine the order and orientation of the alignments
        # alignment_relations = orient_alignments(first_match.subsequent_matches())
        alignment_relations = first_match.orient_alignments()

        # STEP 3: Decide what to do with the fragmented alignments: join or split?
        distinct_alignments = first_match.consolidate_subalignments(alignment_relations)
        first_match.enumerate()
    return distinct_alignments


def format_split_alignments(domain_table: DomainTableParser, search_stats: HmmSearchStats) -> dict:
    """
    Handles the alignments where 'of' > 1
    If the alignment covers the whole target HMM or if the distance between the two parts of the alignment
    are very far apart, then the alignment will be divided into two unrelated alignments
    If the alignment parts are near together and/or each part covers a portion of the HMM, then they will be joined

    :param domain_table: DomainTableParser() object
    :param search_stats: HmmSearchStats() instance containing accumulators to track alignment parsing
    :return:
    """
    # Dictionary of single sequence alignments to return
    distinct_alignments = dict()

    # Query-relevant parsing variables
    first_match = previous_match = HmmMatch()
    previous_query_header = ""
    while domain_table.next():
        data = domain_table.alignments
        hmm_match = HmmMatch()
        hmm_match.target_hmm = data['hmm_name']
        hmm_match.hmm_len = data['hmm_len']
        hmm_match.seq_len = data['query_len']
        hmm_match.orf = data['query']
        hmm_match.desc = data['desc']
        hmm_match.start = data['qstart']
        hmm_match.end = data['qend']
        hmm_match.pstart = data['pstart']
        hmm_match.pend = data['pend']
        hmm_match.num = data['num']
        hmm_match.of = data['of']
        hmm_match.acc = data['acc']  # Used for filtering
        hmm_match.ieval = data['iEval']  # Used for filtering
        hmm_match.eval = data['Eval']  # Used for filtering
        hmm_match.full_score = data['full_score']  # Used for filtering

        search_stats.raw_alignments += 1

        query_header = ' '.join([hmm_match.orf, hmm_match.desc])
        # Finish off "old business" (sub-alignments)
        if previous_match.orf != hmm_match.orf and first_match.orf == previous_match.orf:
            distinct_alignments.update(assemble_domain_alignments(first_match, search_stats))
        if hmm_match.target_hmm != previous_match.target_hmm and query_header == previous_query_header:
            # New HMM (target), same ORF (query)
            search_stats.multi_alignments += 1

        # Carry on with this new alignment
        query_header_desc_aln = ' '.join([hmm_match.orf, hmm_match.desc]) + \
                                '_' + str(hmm_match.num) + '_' + str(hmm_match.of)
        if not hmm_match.orf:
            LOGGER.error("Double-line parsing encountered: hmm_match.orf is empty!\n")
            sys.exit(9)

        if data["of"] == 1:
            distinct_alignments[query_header_desc_aln] = hmm_match
        elif hmm_match.num == 1:
            search_stats.fragmented += 1
            first_match = hmm_match
        else:
            search_stats.fragmented += 1
            previous_match.next_domain = hmm_match

        previous_query_header = query_header
        previous_match = hmm_match

    # Check to see if the last alignment was part of multiple alignments, just like before
    if first_match.next_domain and first_match.orf == previous_match.orf:
        distinct_alignments.update(assemble_domain_alignments(first_match, search_stats))

    return distinct_alignments


def filter_poor_hits(thresholds: namedtuple, distinct_alignments: dict, search_stats: HmmSearchStats) -> dict:
    """
    Filters the homology matches based on their E-values and mean posterior probability of aligned residues from
    the maximum expected accuracy (MEA) calculation.
    Takes into account multiple homology matches of an ORF to a single gene and determines the total length of the
    alignment instead of treating them as individual alignments. This information is used in the next filtering step.

    :param thresholds: A namedtuple with max_e, max_ie, min_acc, min_score and perc_aligned attributes
     that must be exceeded for alignments to be included.
    :param distinct_alignments: A dictionary of HmmMatch instances indexed by their respective header names
    :param search_stats: An HmmSearchStats instance used for tracking various alignment parsing stats
    :return: A dictionary of HmmMatch instances that pass thresholds indexed by their respective header names
    """

    purified_matches = dict()

    for query_header_desc_aln in sorted(distinct_alignments):
        hmm_match = distinct_alignments[query_header_desc_aln]

        query_header_desc = (hmm_match.orf, hmm_match.desc)
        if query_header_desc not in purified_matches:
            purified_matches[query_header_desc] = list()

        if hmm_match.eval <= float(thresholds.max_e) and hmm_match.ieval <= float(thresholds.max_ie):
            if hmm_match.acc >= float(thresholds.min_acc) and hmm_match.full_score >= float(thresholds.min_score):
                purified_matches[query_header_desc].append(hmm_match)
                continue
        search_stats.dropped += 1
        search_stats.bad += 1

    return purified_matches


def filter_incomplete_hits(thresholds: namedtuple, purified_matches: dict, search_stats: HmmSearchStats) -> list:
    """
    Removes all alignments of each query-HMM pair that do not meet the thresholds.perc_alignment cut-off.
    The alignment length is based on the start and end positions on the HMM profile, not the query sequence.

    :param thresholds:  A namedtuple with max_e, max_ie, min_acc, min_score and perc_aligned attributes
     that must be exceeded for alignments to be included.
    :param purified_matches: A dictionary of HmmMatch instances indexed by their respective header names
    :param search_stats: An HmmSearchStats instance for tracking the number of sequences filtered out
    :return: List of HmmMatch instances that meet or exceed the minimum percentage aligned threshold to be retained
    """
    complete_gene_hits = list()

    for query in purified_matches:
        for hmm_match in purified_matches[query]:  # type: HmmMatch
            ali_len = hmm_match.pend - hmm_match.pstart  # length of alignment on profile
            query_aligned = (hmm_match.end-hmm_match.start)*100/hmm_match.seq_len  # query seq % aligned to profile
            perc_aligned = (float((int(ali_len)*100)/int(hmm_match.hmm_len)))  # % of profile covered by alignment
            if query_aligned >= thresholds.query_aligned and thresholds.profile_match is False:
                complete_gene_hits.append(hmm_match)
            elif perc_aligned >= thresholds.perc_aligned:
                complete_gene_hits.append(hmm_match)
            else:
                search_stats.dropped += 1
                search_stats.short += 1

    return complete_gene_hits


def renumber_multi_matches(complete_gene_hits: list):
    orf_name_index = dict()
    # Find which names have duplicates
    for match in sorted(complete_gene_hits, key=lambda x: x.orf):  # type: HmmMatch
        if match.orf not in orf_name_index:
            orf_name_index[match.orf] = []
        orf_name_index[match.orf].append(match)
    # Rename the ORFs with more than one alignment
    for orf_name in orf_name_index:
        n = 1
        of = len(orf_name_index[orf_name])
        for match in sorted(orf_name_index[orf_name], key=lambda x: x.num):
            match.num = n
            match.of = of
            n += 1
    return


def drop_current_match(match: HmmMatch) -> HmmMatch:
    i = 0
    while i < len(match.subsequent_matches()):
        if match.subsequent_matches()[i] == match:
            match.subsequent_matches().pop(i)
            i = len(match.subsequent_matches())
        i += 1
    match = match.next_domain  # type: HmmMatch
    match.of -= 1
    match.num -= 1
    return match


def drop_next_match(match: HmmMatch) -> None:
    """
    Function for removing an element from the linked list formed by HmmMatch.next_domain

    :param match: The current HmmMatch of which the next is to be removed from the linked list subsequent_matches
    :return: None
    """
    i = 0
    while i < len(match.subsequent_matches()):
        if match.subsequent_matches()[i] == match.next_domain:
            match.subsequent_matches().pop(i)
            i = len(match.subsequent_matches())
            try:
                match.next_domain = match.subsequent_matches()[i+1]
            except IndexError:
                match.next_domain = None
            match.of -= 1
            return
        i += 1
    LOGGER.warning("Next HmmMatch was not dropped from linked list.\n")
    return


def remove_redundant_alignments(match: HmmMatch, index=0) -> HmmMatch:
    """
    If multiple alignments exist, HMMER will not report alignments that overlap in both profile and query sequence,
    however, there is a chance some redundantly cover a high proportion of the HMM profile.
    It is the purpose of this function to identify those alignments and eliminate the worse of the two.

    Algorithm:
        For each alignment in the linked-list of alignments (subsequent_matches):
         calculate HMM profile alignment overlaps
          if non-zero:
           Worse of the two alignments is determined (based on profile alignment length and score).
           Worst alignment is dropped from the linked-list

    This function is currently used as a pre-filter for alignment scaffolding and inversion detection.
    It clears out the riff raff that cannot be scaffolded (since profile alignment co-ordinates overlap entirely)
    and shouldn't be maintained because the alignment is garbage. Consequences of keeping these spurious alignments
    include, but are not limited to, false alarms during inversion detection (with HmmMatch.contains_inversion).

    :return: HmmMatch pointing to the head of the linked list
    """
    if not match.next_domain:
        return match
    match.next_domain = remove_redundant_alignments(match.next_domain, index+1)  # Patched in 0.8.9
    query_orientation = detect_orientation(match.start, match.end,
                                           match.next_domain.start, match.next_domain.end)
    profile_orientation = detect_orientation(match.pstart, match.pend,
                                             match.next_domain.pstart, match.next_domain.pend)

    if query_orientation == "satellite":
        if profile_orientation == "subsequence":
            if ((match.pend - match.pstart)/match.hmm_len) < 0.1:
                match = drop_current_match(match)
        elif profile_orientation == "supersequence":
            if ((match.next_domain.pend - match.next_domain.pstart)/match.hmm_len) < 0.1:
                drop_next_match(match)

    return match

#include "helper.h"
using namespace std;


#define _MAX 1000000000

unsigned long create_contigs_dictionary(std::string contigs_file, std::map<std::string, CONTIG> &contigs_dictionary) {
    
     FastaReader fastareader(contigs_file);
     map<string, unsigned long> contig_lengths;
     map<string, unsigned long>::iterator it_contig_lens;
     
     fastareader.get_fasta_sequence_info(contig_lengths);
     unsigned long genome_length = 0;
     CONTIG contig;
     for(it_contig_lens = contig_lengths.begin(); it_contig_lens != contig_lengths.end(); it_contig_lens++ ) {
        genome_length += it_contig_lens->second;
        contig.L = it_contig_lens->second;
        contigs_dictionary[it_contig_lens->first] = contig; 
     }

     return genome_length;
}


RUN_STATS consume_sam(const std::string &SAM_file, const std::string &format,\
     vector<MATCH> &all_reads, map<std::string, float > &multireads, bool show_status) {

    MatchOutputParser *parser = ParserFactory::createParser(SAM_file, format);
    if ( parser == 0 ) {
        std::cout << "ERROR : Cannot open a parser to parse the file " << SAM_file << std::endl;
    }

    MATCH match;

    map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > reads_dict;
    vector<std::string> holder;
   
    if ( show_status )
        std::cout << "Number of lines processed: " << std::endl;

    RUN_STATS stats;

    int i;
    struct QUADRUPLE <bool, bool, unsigned int, unsigned int> p;
    for ( i =0; ; i++ ) {
        if (show_status && i % 10000 == 0)
            std::cout << "\n\033[F\033[J" << i;

        if (!parser->nextline(match))
            break;

        if (i >= _MAX) break;

        if (match.mapped)
            stats.num_mapped_reads++;
        else
            stats.num_unmapped_reads++;

        if (match.parity)
            stats.num_reads_2++;
        else stats.num_reads_1++;

        if (reads_dict.find(match.query) == reads_dict.end()) {
            p.first = false;
            p.second = false;
            p.third = 0;
            p.fourth = 0;
            reads_dict[match.query] = p;
        }
        stats.num_alignments++;

        // if it is not mapped then ignore it
        if (!match.mapped)
            continue;

        if (match.parity) {
            reads_dict[match.query].first = true;
            reads_dict[match.query].third++;
        }
        else {
            reads_dict[match.query].second = true;
            reads_dict[match.query].fourth++;
        }

        // store it to process later by looking up the dictionary
        try {
            all_reads.push_back(match);
        }
        catch (...) {
            cout << "failing " << match.query << "   " << all_reads.size() << endl;
        }

    }

    for( map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> >::iterator it = reads_dict.begin();
         it != reads_dict.end();
         it++) {
        if( !(it->second.first && it->second.second))
            stats.num_singleton_reads++;
        if( it->second.third > 1) {
            stats.num_multireads++;
            multireads[it->first] = 0.0;
            stats.num_secondary_hits += it->second.third-1;
        }
        if( it->second.fourth  > 1) {
            stats.num_multireads++;
            multireads[it->first] = 0.0;
            stats.num_secondary_hits += it->second.fourth-1;
        }
    }

    stats.num_distinct_reads_unmapped = stats.num_unmapped_reads;
    stats.num_distinct_reads_mapped = stats.num_mapped_reads - stats.num_secondary_hits;

    for ( vector<MATCH>::iterator it = all_reads.begin(); it != all_reads.end(); it++)  {

        if ( it->parity == 0  ) {
            if( reads_dict[it->query].first && reads_dict[it->query].second )
                it->w = 0.5/static_cast<float>(reads_dict[it->query].third);
            else
                it->w = 1/static_cast<float>(reads_dict[it->query].third);
        }
        else  { //parity 1
            if( reads_dict[it->query].first && reads_dict[it->query].second )
                it->w = 0.5/static_cast<float>(reads_dict[it->query].fourth);
            else
                it->w = 1/static_cast<float>(reads_dict[it->query].fourth);
        }
    }

    delete parser;
    return stats;
}


void process_SAM(const std::string & reads_map_file, std::map<string, CONTIG> &contigs_dictionary,
                 const std::string &reads_map_file_format,
                 std::vector<MATCH> &all_reads,
                 map<std::string, float > &multireads,
                 bool show_status) {
    /*
     * At this point, multireads is a map with only the keys (read names) entered and
     * has empty vector<string> objects as values - they are populated here.
     */
    MATCH match;

    READ_DATA read_data;
    READ_DATUM read_datum;
    
    int i =0;
    
    if( show_status ) std::cout << "Number of hits processed : " ;
    // iterate through individual hits/alignments
    for(vector<MATCH>::iterator it=all_reads.begin();  it!= all_reads.end(); it++ )  {

        if( i >=_MAX ) break;

        if (show_status && i%10000==0) {
           std::cout << "\n\033[F\033[J";
           std::cout << i ;
        }
        i++;

        if ( contigs_dictionary.find(it->subject)==contigs_dictionary.end() ) {
            std::cout << " Missing contig " << it->subject << std::endl;
            std::cerr << "ERROR : Could not find the matched contig in the contig file " << std::endl;
            exit(1);
        }
        read_datum.name = it->query;
        if ( multireads.find(it->query) != multireads.end() ) {
            multireads[it->query] += 1;
            read_datum.multi = true;
        }
        else
            read_datum.multi = false;

        if (  it->start < it->end ) {
            read_datum.start = it->start;
            read_datum.end = it->end;
        }
        else {
            read_datum.start = it->end;
            read_datum.end = it->start;
        }

        contigs_dictionary[it->subject].M.push_back(read_datum);
    }
}

unsigned int getMaxReadSize( std::map<string, vector<MATCH> > &orf_dictionary,
                             std::map<string, CONTIG> &contigs_dictionary) {
    unsigned int size = 0;
    std::map<string, vector<MATCH> >::iterator itcont;

    for(itcont= orf_dictionary.begin(); itcont != orf_dictionary.end(); itcont++)  {
       for(std::vector<READ_DATUM>::iterator it = contigs_dictionary[itcont->first].M.begin(); it != contigs_dictionary[itcont->first].M.end(); it++) {
          if( size < it->end  - it->start ) size = it->end  - it->start ;
       }
    }

    return size;
} 


std::vector<READ_DATUM>::iterator  binary_search(std::vector<READ_DATUM> &A, int seekValue) {
  // continually narrow search until just one element remains
  unsigned long imin, imax;
  imin = 0;
  imax = A.size();

  while (imin < imax)
    {
      unsigned int imid = (imin+imax)/2;
 
      // code must guarantee the interval is reduced at each iteration
      assert(imid < imax);

      // note: 0 <= imin < imax implies imid will always be less than imax
 
      // reduce the search
      if (A[imid].start < static_cast<unsigned int>(seekValue) )
        imin = imid + 1;
      else
        imax = imid;
    }

    std::vector<READ_DATUM>::iterator it = A.begin() + imin;
    return it ;
}

void substring_coverage(std::map<string, CONTIG> &contigs_dictionary, const std::string &contig,
                        unsigned long start, unsigned long end,
                        COVERAGE &coverage, unsigned int maxReadLength,
                        map<std::string, float > &multireads, bool multi_reads) {

    if ( contigs_dictionary.find(contig) == contigs_dictionary.end() || contigs_dictionary[contig].L == 0 ) {
        coverage.coverage = 0 ;
        coverage.numreads = 0 ;
        coverage.sequence_length = end - start ;
        coverage.uncovered_length = 0;
    }

    float numreads = 0;
    float _coverage = 0;
    unsigned long uncovered_length = 0;
    unsigned long p_end = start;

    int _seekValue =  maxReadLength == 0 ? 0 :  (start < maxReadLength ? 0 : start-maxReadLength );

    std::vector<READ_DATUM>::iterator it= contigs_dictionary[contig].M.begin();

    if ( _seekValue >0 )
        it =  binary_search(contigs_dictionary[contig].M, _seekValue);

    // iterate through every read that aligned to that contig
    for ( ; it != contigs_dictionary[contig].M.end(); it++) {

        uncovered_length  +=  ( p_end > it->start  || it->start > end) ? 0 : it->start - p_end;
        //make sure the read start and end are not going beyond the contig
        if( it->end > p_end )
            p_end = it->end;

        if( (start <= it->start && it->start <= end) ||  (start <= it->end && it->end <= end)  ) {
            numreads += 1;
        }

        // the subsequent reads are going off the end of the orf
        if( it->start > end )
            break;

        if (multi_reads && it->multi) {
            float read_multiplicity = multireads.find(it->name)->second;
            numreads += 1.0/read_multiplicity;
        }
    }


    uncovered_length += (p_end > end ) ? 0 : (end - p_end);

    unsigned long sequence_length = end - start;
    if( sequence_length > 0 )
        _coverage = ((float)(sequence_length - uncovered_length )/(float)sequence_length)*100;

    coverage.numreads = numreads;
    coverage.coverage = _coverage;
    coverage.sequence_length  = sequence_length;
    coverage.uncovered_length =  uncovered_length;
}
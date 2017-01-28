#include <map>
#include <vector>
#include <algorithm>
#include "types.h"
#include "utilities.h"
#include "helper.h"
#include <assert.h>

using namespace std;
bool compare_triplets(const READ_DATUM &a, const READ_DATUM &b) {
   return a.start < b.start ? true : false; 
}
 
int main( int argc, char **argv ){
    // parse options
    Options options;

    if (options.SetOptions(argc, argv)==false) {
        options.print_usage(argv[0]);
        exit(0);
    }
    map<string, CONTIG> contigs_dictionary;
    unsigned long total_contig_length = create_contigs_dictionary(options.contigs_file, contigs_dictionary);
    
    bool print_stats_file = false;

    std::ostream *output;
    std::ofstream rpkm_output;
    rpkm_output.open(options.output_file.c_str(), std::ifstream::out);
    output = &rpkm_output;

    std::ostream *stats_out;
    std::ofstream Sout;
    if ( options.stats_file.size() > 0) {
        print_stats_file = true; 
        Sout.open(options.stats_file.c_str(), std::ifstream::out);
        stats_out = &Sout;
    }   
 
    vector<MATCH> all_reads;
    all_reads.reserve(80000000);
    
    // creating the read multiplicity counts if there is a single end read file 
    map<std::string, float > multireads;

    RUN_STATS stats;
    RUN_STATS _stats;
    for( vector<string>::iterator it = options.read_map_files.begin() ; it != options.read_map_files.end(); it++)  { 
        if( it->size() == 0 ) 
            continue;
        all_reads.clear();
        if( options.show_status ) 
            std::cout << "\n\n" << "Parsing alignments from SAM file " << *it << std::endl;
        _stats = detect_multireads_samfile(*it, options.reads_map_file_format, all_reads, multireads, options.show_status);
        
        if( options.show_status )
            _stats.printStats(&std::cout);

        if(print_stats_file) {
            *stats_out << "\nStats for file :  " << *it << std::endl;
            _stats.printStats(stats_out);
        }

        stats = stats + _stats;
        process_SAM(*it, contigs_dictionary, options.reads_map_file_format, all_reads, multireads, options.show_status);
        stats.num_multireads = multireads.size();
    }
   
    if ( options.show_status )
       std::cout << "Composite stats for all files " << std::endl;

    if ( print_stats_file)
       *stats_out << "\nComposite stats for all files " << std::endl;
 
    if ( options.show_status )
        stats.printStats(&std::cout);

    if (print_stats_file)
        stats.printStats(stats_out);
 

    if ( options.show_status )
        std::cout << "\n\nSorting the read matches...";
    for (map<string, CONTIG>::iterator it = contigs_dictionary.begin(); it != contigs_dictionary.end(); it++) {
        std::sort(it->second.M.begin(), it->second.M.end(), compare_triplets);

    }
    if (options.show_status)
        std::cout << "done" << endl;

//    float total_distinct_reads;
    if (options.num_reads == 0) {
//        total_distinct_reads = stats.num_distinct_reads_mapped + stats.num_distinct_reads_unmapped;
    }
    else {
//        total_distinct_reads = static_cast<double>(stats.num_distinct_reads_mapped);
        stats.num_total_reads = options.num_reads;
    }

    unsigned long total_covered_length = 0;
    float rpkm_sum = 0.0;
    COVERAGE coverage;
    if (options.show_status)
        std::cout << "Computing coverage and RPKM values for all contigs...";

    for (map<string, CONTIG>::iterator it = contigs_dictionary.begin(); it != contigs_dictionary.end(); it++) {
        substring_coverage(contigs_dictionary, it->first, 1, it->second.L, coverage, 0, multireads, options.multi_reads);
        total_covered_length += coverage.coverage*contigs_dictionary[it->first].L;
        it->second.rpkm = (1E9/ static_cast<float>(stats.num_total_reads))*
                (static_cast<float>(coverage.numreads)/ static_cast<float>(coverage.sequence_length));
        rpkm_sum += it->second.rpkm;

        *output << it->first << ',' << it->second.rpkm << endl;
    }
    rpkm_output.close();

    if (options.show_status)
        std::cout << "done\n";

    char buf[1000];

    std::cout << std::endl;
    sprintf(buf, "Number of Contigs                         : %ld ", (long int)contigs_dictionary.size() );
    std::cout << buf << std::endl;

    sprintf(buf,"Total Contig Length                       : %ld ", total_contig_length);
    std::cout << buf  << std::endl;

    sprintf(buf,"Contig length covered                     : %ld ", total_covered_length/100);
    std::cout << buf  << std::endl;

    sprintf(buf,"Percentage of contig coverage             : %-5.2f%%",
            (float)total_covered_length/(float)total_contig_length);
    std::cout << buf  << std::endl;

    sprintf(buf,"Percentage of uniquely-aligned reads      : %-5.2f%%",
            100*(float)stats.num_distinct_reads_mapped/stats.num_total_reads);
    std::cout << buf  << std::endl;

    sprintf(buf,"Percentage of multireads                  : %-5.2f%%",
            100*(float)stats.num_multireads/stats.num_total_reads);
    std::cout << buf  << std::endl;

    sprintf(buf,"Percentage of unaligned reads             : %-5.2f%%",
            100*(float)stats.num_distinct_reads_unmapped/stats.num_total_reads) ;
    std::cout << buf  << std::endl;

    sprintf(buf,"Percentage of mapped reads                : %-5.2f%%",
            100*(static_cast<double>(stats.num_mapped_reads)/stats.num_total_reads )) ;
    std::cout << buf  << std::endl;

    sprintf(buf,"Total number of reads                     : %ld ",stats.num_total_reads);
    std::cout << buf  << std::endl;

    sprintf(buf,"Number of multireads                      : %d ", stats.num_multireads);
    std::cout << buf  << std::endl;

    sprintf(buf,"Avg rpkm across contigs in sample         : %.2f ", rpkm_sum/contigs_dictionary.size());
    std::cout << buf  << std::endl;

   if ( print_stats_file ) {
       *stats_out << std::endl;
       *stats_out << "Mapping stats for the sample " << std::endl;
       sprintf(buf, "Number of Contigs                         : %ld ", (long int)contigs_dictionary.size() );
       *stats_out << buf<< std::endl;

       sprintf(buf,"Total Contig Length                       : %ld ", total_contig_length);
       *stats_out << buf  << std::endl;

       sprintf(buf,"Contig length covered                     : %ld ", total_covered_length/100);
       *stats_out << buf  << std::endl;

       sprintf(buf,"Percentage of contig coverage             : %-5.2f%%",
               (float)total_covered_length/(float)total_contig_length);
       *stats_out << buf  << std::endl;

       sprintf(buf,"Percentage of uniquely-aligned reads      : %-5.2f%%",
               100*(float)stats.num_distinct_reads_mapped/stats.num_total_reads);
       *stats_out << buf  << std::endl;

       sprintf(buf,"Percentage of multireads                  : %-5.2f%%",
               100*(float)stats.num_multireads/stats.num_total_reads);
       *stats_out << buf  << std::endl;

       sprintf(buf,"Percentage of unaligned reads             : %-5.2f%%",
               100*(float)stats.num_distinct_reads_unmapped/stats.num_total_reads) ;
       *stats_out << buf  << std::endl;

       sprintf(buf,"Percentage of mapped reads                : %-5.2f%%",
               100*(static_cast<double>(stats.num_mapped_reads)/stats.num_total_reads )) ;
       *stats_out << buf  << std::endl;

       sprintf(buf,"Total number of reads                     : %ld ",stats.num_total_reads);
       *stats_out << buf  << std::endl;

       sprintf(buf,"Number of multireads                      : %d ", stats.num_multireads);
       *stats_out << buf  << std::endl;

       sprintf(buf,"Avg rpkm across contigs in sample         : %.2f ", rpkm_sum/contigs_dictionary.size());
       *stats_out << buf  << std::endl;

       Sout.close();
   }

}
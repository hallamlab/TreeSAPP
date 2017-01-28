#ifndef __RPKM_TYPE
#define __RPKM_TYPE
#include <map>
#include <vector>
#include <iostream>

using namespace std;

typedef struct _READ_DATUM {
    unsigned int start,
            end;
    bool multi;
    std::string name;
} READ_DATUM;

typedef vector<READ_DATUM> READ_DATA;

typedef struct _CONTIG {
    /*
     * L is the length of the contig
     * rpkm is the contig's RPKM value - calculated in main
     * M is the list of all reads aligned to the contig with start, end, boolean for multiread and its name
     */
    unsigned int L;
    float rpkm;
    READ_DATA M;
} CONTIG;


typedef struct _MATCH {
    /*
     * parity is 1 if it is second in the pair
     * mapped is 1 if the read was not unmapped
     * orphan is 1 if the mate was not successfully aligned
     * multi is 1 if the read is not a primary alignment
     * chimeric is 1 if parts of the read aligned to different loci
     * singleton is 1 if the mate was not successfully aligned
     */
    std::string query, subject;
    unsigned int start, end;
    bool parity; 
    bool mapped;
    bool orphan;
//    bool multi;
    bool chimeric;
    bool singleton;
    float  w; //no idea what this does...
    _MATCH(): w(0) { } 
} MATCH;


template< typename A, typename B, typename C, typename D>
struct QUADRUPLE {
     A first;
     B second;
     C third;
     D fourth;
};


typedef struct _RUN_STATS {
    int num_unmapped_reads ; 
    int num_mapped_reads ; 
    int num_singleton_reads ; 
    int num_reads_1 ; 
    int num_reads_2 ; 
    int num_multireads ;
    long int num_total_reads ;
    int num_secondary_hits ;
    int num_distinct_reads_unmapped;
    int num_distinct_reads_mapped;
    
    _RUN_STATS() {
        num_unmapped_reads =0; 
        num_mapped_reads =0; 
        num_singleton_reads =0; 
        num_reads_1 = 0;
        num_reads_2 = 0;
        num_multireads = 0;
        num_total_reads = 0;
        num_secondary_hits = 0;
        num_distinct_reads_unmapped =0 ;
        num_distinct_reads_mapped = 0;
    }

    struct _RUN_STATS&  operator+( const struct _RUN_STATS & stats) {
        this->num_unmapped_reads += stats.num_unmapped_reads;  
        this->num_mapped_reads  += stats.num_mapped_reads  ;
        num_singleton_reads += stats.num_singleton_reads ;
        num_reads_1 += stats.num_reads_1;  
        num_reads_2 += stats.num_reads_2; 
        num_multireads += stats.num_multireads; 
        num_total_reads  += stats.num_total_reads  ;
        num_secondary_hits  += stats.num_secondary_hits;
        num_distinct_reads_unmapped  += stats.num_distinct_reads_unmapped  ;
        num_distinct_reads_mapped  += stats.num_distinct_reads_mapped  ;
        return *this;
    }


    void printStats( std::ostream *output ) {
        *output << std::endl;
        *output << "Number of lines in SAM file      : " << num_total_reads << std::endl; 
        *output << "Number of unmapped lines         : " << num_unmapped_reads << std::endl; 
        *output << "Number of mapped lines           : " << num_mapped_reads << std::endl; 
        *output << "Number of singletons             : " << num_singleton_reads << std::endl; 
        *output << "Number of read 1                 : " << num_reads_1  << std::endl; 
        *output << "Number of read 2                 : " << num_reads_2 << std::endl; 
        *output << "Number of multireads             : " << num_multireads << std::endl; 
        *output << "Number of secondary hits         : " << num_secondary_hits << std::endl; 
        *output << "Number distinct reads mapped     : " << num_distinct_reads_mapped << std::endl; 
        *output << "Number distinct reads unmapped   : " << num_distinct_reads_unmapped << std::endl; 
    }



} RUN_STATS;

typedef struct _COVERAGE {
    float coverage;
    float numreads;
    unsigned int sequence_length,
            uncovered_length;
} COVERAGE;
#endif //__RPKM_TYPE




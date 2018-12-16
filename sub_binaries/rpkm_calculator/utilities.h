#ifndef _UTILITIES
#define _UTILITIES
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <string.h>
#include <stdlib.h>

using namespace std;

// Structure for RPKM input options
struct Options {
   
    /* Input files for RPKM */
    string contigs_file; // the reference sequence file used for indexing and alignment
    string output_file; // location to write CSV file
    string stats_file; // location to write summary stats file
    string reads_map_file_format; // aligner type BWA or BLAST, two SAM files or one
    vector<string> read_map_files;
    long int num_reads;
    
    /* Flags and settings */
    bool multi_reads; // flag for detecting multiple mapping of reads
    bool show_status; // shows the counter that counts the number of reads processed, and other info on screen
    
    // Constructor with default settings 
    Options() {
        contigs_file = "";
        stats_file = "";
        output_file = "";
        read_map_files.clear();
        num_reads = 0;

        multi_reads = false;
        show_status = false;
        reads_map_file_format = "sam-1";
    };
    
    void print_usage( char *arg);
 
    bool SetOptions(int argc, char *argv[]);
};

void split(const std::string  &strn, std::vector<char *> &v, char *buf, char d='\t');

bool matchString(const string &str, const string & stringtomatch, bool fromstart=false);

void get_fasta_sequence_info(const std::string &fasta_file_name);

std::string extract_sequence_name(const std::string &name);

#endif //_UTILITIES


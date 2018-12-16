#ifndef _HELPER
#define _HELPER
#include <map>
#include <string>
#include <ostream>
#include <iterator>
#include <assert.h>
#include "types.h"
#include "fastareader.h"
#include "matchoutputparser.h"

using namespace std;

unsigned long create_contigs_dictionary(std::string contigs_file, std::map<std::string, CONTIG> &contigs_dictionary);

RUN_STATS consume_sam(const std::string &SAM_file,
        const std::string &format,
        vector <MATCH> &all_reads,
        map<std::string, float > &multireads,
        bool show_progress_counter = false);
/*
 * multireads is a map object containing the header of the read in the first position and
 * the names of all contigs in a string vector
 */

void process_SAM(const std::string & reads_map_file, std::map<string, CONTIG> &contigs_dictionary,
                 const std::string &reads_map_file_format,
                 vector<MATCH> &all_reads,
                 map<std::string, float > &multireads,
                 bool show_status= false);

void substring_coverage(std::map<string, CONTIG> &contigs_dictionary, const std::string &contig,
                        unsigned long start, unsigned long end,
                        COVERAGE &coverage, unsigned int maxReadLength,
                        map<std::string, float > &multireads, bool multi_reads);

#endif //_HELPER

#ifndef __FASTAREADER
#define __FASTAREADER
#include "utilities.h"
#include <string.h>


class FastaReader {

private:
     string contigs_file;
public:
     FastaReader(const string & contigs_file);
     void get_fasta_sequence_info(map<string, unsigned long> &contigs_dictionary) ;
     std::string extract_sequence_name(const std::string &name);
     string getContigsFileName();
};

#endif // __FASTAREADER

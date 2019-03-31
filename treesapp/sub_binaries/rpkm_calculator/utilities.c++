#include "utilities.h"
#include <stdlib.h>


void Options::print_usage(char *arg) {
    std::cout << "\nUSAGE : " << arg << " -c contigs.fasta -a alignments.sam [-a alignments_N.sam] -o RPKMs.csv\n\n"\
    << "Required arguments:   \n"
    << "-c  <contigs.fasta>   \n"\
    << "-a  <alignments.sam>  \n"\
    << "-o  <output_file.csv> \n\n"\
    << "Optional arguments:   \n"
    << "-n  <num_reads>       [default = number of reads in SAM]\n"\
    << "-s  <stats_file>      \n\n"\

    << "--m                   \n"\
    << "\tUse this flag to use multireads in RPKM calculations "\
    << "(i.e., turn off masking of low-complexity regions).\n"\
    << "\tReads aligning to multiple loci will evenly distributed amongst the contigs.\n"
    << "--verbose             [shows run status]\n"\
    << "--h                   [ for help ]\n"\
    << std::endl;
}

bool Options::SetOptions(int argc, char *argv[]) {
   for(int i = 1; i < argc ; i++) {
       if( strncmp(argv[i], "-c", strlen("-c")) == 0 ) {
          this->contigs_file = argv[++i];
       }   
       else if( strncmp(argv[i], "--h", strlen("--h")) == 0 ) {
           print_usage(argv[0]);
           exit(0);
       }   
       else if( strncmp(argv[i], "-s", strlen("-s")) == 0 ) {
          this->stats_file = argv[++i];
       }
       else if(strncmp(argv[i], "-a", strlen("-a")) == 0 ) {
          this->read_map_files.push_back(string(argv[++i]));
       }
       else if(strncmp(argv[i], "--verbose", strlen("--verbose")) == 0 ) {
          this->show_status = true;
       }
       else if( strncmp(argv[i], "-o", strlen("-o")) == 0 ) {
           this->output_file = argv[++i];
       }
       else if( strncmp(argv[i], "--m", strlen("--m")) == 0  ) {
          this->multi_reads = true;
       }
       else if( strncmp(argv[i], "-n", strlen("-n")) == 0  ) {
           this->num_reads = atoi(argv[++i]);
       }
    } //for loop for arguments processing

   if ( this->contigs_file.empty()) {
       std::cout << "ERROR: There must be a contigs file" << std::endl;
       return false;
   }

   if ( this->read_map_files.empty()) {
       std::cout << "ERROR: There must be at least one SAM file" << std::endl;
       return false;
   }

   if ( this->output_file.empty()) {
       std::cout << "ERROR: There must be a output file path" << std::endl;
       return false;
   }
   return true;
};

void get_fasta_sequence_info(const std::string & fasta_file_name) {
    /*
     * This function prints a map of each sequence's name and its length
     */
    std::ifstream input(fasta_file_name.c_str());
    if(!input.good()){
        std::cerr << "Error opening '" << fasta_file_name << "'. Bailing out." << std::endl;
        return ;
    }
 
    std::string line, name, content;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Sequence record name
            if( !name.empty() ){ // Print out what we read from the last entry
                std::cout << extract_sequence_name(name) << " : " << content.size() << std::endl;
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){ // Print out what we read from the last entry
        name = extract_sequence_name(name);
        std::cout << name << " : " << content.size() << std::endl;
    }
    input.close();
 
}

std::string extract_sequence_name(const std::string &name) {
     char  cstr[1000000];
     strcpy(cstr, name.c_str());
     
     char * cptr = cstr;

     while( *cptr != '\t' && *cptr !=  ' ' && *cptr != '\0' )  cptr++; 
     (*cptr) ='\0';

     std::string sname(cstr);
     return sname;
}

/*
CML -- this function has been made robust to lines that are >1000 characters
by allocating more space for buf if required
*/
void split(const string  &strn, std::vector<char *> &v, char *buf, char d) {
  if (strn.length() > 1000 )
      buf = (char *) malloc ((strn.length() + 1) * sizeof(char));
  strcpy(buf, strn.c_str());
  char *s1 = buf;
  v.clear();
  v.reserve(15);
  v.push_back(s1);
  while(*s1 != '\0') {
     if(*s1==d) { 
       *s1 = '\0';
       v.push_back(s1+1);
     }
     s1++;
  }
}

bool matchString(const string &str, const string & stringtomatch, bool fromstart) {

    unsigned long pos = str.find(stringtomatch);
    if(fromstart && pos ==0 ) return true;

    return !fromstart && pos >= 0;

}

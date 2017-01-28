#include "utilities.h"
#include "fastareader.h"
//#include <iostream>
//#include <fstream>

FastaReader::FastaReader( const string & contigs_file) {
    this->contigs_file = contigs_file;
}



string FastaReader::getContigsFileName() {
    return this->contigs_file;
}

void FastaReader::get_fasta_sequence_info(map< string, unsigned long> &contigs_dictionary) {

    std::ifstream input(this->contigs_file.c_str());
    if(!input.good()){
        std::cerr << "Error opening '"<<this->contigs_file << "'. Bailing out." << std::endl;
        return ;
    }
 
    std::string line, name, content;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                name = extract_sequence_name(name);
                contigs_dictionary[name] = content.size() ;
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
        //std::cout << name << " : " << content.size() << std::endl;
        contigs_dictionary[name] = content.size() ;
    }
    input.close();
 
}

std::string FastaReader::extract_sequence_name(const std::string &name) {
     char  cstr[1000000];
     strcpy(cstr, name.c_str());
     
     char * cptr = cstr;

     while( *cptr != '\t' && *cptr !=  ' ' && *cptr != '\0' )  cptr++; 
     (*cptr) ='\0';

     std::string sname(cstr);
     return sname;
}

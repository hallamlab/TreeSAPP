#include "matchoutputparser.h"

using namespace std;

MatchOutputParser::MatchOutputParser(const std::string &filename, const std::string &format) {
     this->filename = filename;
     this->format = format;
     this->num_unmapped_reads =0;
};

unsigned long MatchOutputParser::get_Num_Unmapped_Reads() {
    return  this->num_unmapped_reads;
}

MatchOutputParser::~MatchOutputParser() {
}


SamFileParser::SamFileParser(const std::string &filename, const std::string &format):MatchOutputParser(filename, format) {
     this->input.open(filename.c_str(), std::ifstream::in);

     if(!this->input.good()){
         std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
         return ;
     }  
     this->count = 0;

}

SamFileParser::~SamFileParser() {
   this->input.close();
}

bool SamFileParser::getMateInfo(unsigned int i, MATCH &match)  {
    /*
     * Details for the match object are outlined in types.h
     */

    unsigned int a = i;
    bool singleton = 0;
    a = a >> 2;
    match.mapped = !(a&1); 
    singleton = a&1;

    a = a >> 1;
    match.orphan = a&1; 
    singleton = singleton^(a&1);

    a = a >> 3;
    if ( a&1 )  {
         match.parity = 0; 
         a = a >> 1;
    }
    else {
         a = a >> 1;
         if ( a&1 )
             match.parity  = 1;
         else
             return false;
    }

    a = a >> 4;
    match.chimeric = a&1; 
    match.singleton = singleton;
    return true;

}
bool SamFileParser::nextline(MATCH &match) {
     string line;
     std::string skipPattern("@");
     std::string skipStar("*");

     bool _success = false;
     while ( std::getline(this->input, line ).good()) {
         if ( matchString(line, skipPattern, true) )
             continue;

         fields.clear();
         split(line, fields, this->buf,'\t');

         if(fields.size() < 9)  continue;

         _success = true;
         break;
     }
    
     if ( _success )  {
         match.query =  fields[0];
         match.subject = std::string(fields[2]);
         match.start = atoi(fields[3]);
         match.end =  match.start + std::string(fields[9]).size();
         getMateInfo(static_cast<unsigned int>(atoi(fields[1])), match);

         return true;
     }
     
    return false;
}
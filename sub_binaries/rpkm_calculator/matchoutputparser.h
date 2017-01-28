#ifndef _MATHOUTPUTPARSER
#define _MATHOUTPUTPARSER
#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "utilities.h"
#include "types.h"


using namespace std;


class MatchOutputParser {

public:
      std::string filename;
      std::string format;
      std::ifstream input;
      char buf[1000];
      vector<char *> tempv;
      virtual ~MatchOutputParser() = 0;
      vector<char *> fields;
      MatchOutputParser(const std::string &filename, const std::string &format);
      unsigned long get_Num_Unmapped_Reads();
      virtual bool nextline(MATCH &match)=0;
protected:
      unsigned long num_unmapped_reads;
};


//subclass of the MathOutputParser
class SamFileParser: virtual public MatchOutputParser {
private:
      unsigned long count;
public:
      SamFileParser(const std::string &filename, const std::string &format);
      virtual bool nextline(MATCH &match);
      bool getMateInfo(unsigned int i, MATCH &match);
      ~SamFileParser();
};


class ParserFactory {
public:
     static MatchOutputParser * createParser(const std::string &filename, const std::string &format) {
            if(format.find("sam-1") != string::npos  || format.find("sam-2") !=string::npos )  {
                return new SamFileParser( filename, format);
            }

            return 0;
     }
};



#endif //_MATHOUTPUTPARSER

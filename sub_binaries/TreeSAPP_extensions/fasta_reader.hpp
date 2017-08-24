#ifndef FASTA_H
#define FASTA_H

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>

#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <set>
#include <map>

#include <Python.h>

using namespace std;

static PyObject *read_format_fasta(PyObject *self, PyObject *args);

static char module_docstring[] =
        "This module rapidly loads and formats sequences from a FASTA file using C from within TreeSAPP";
static char read_format_fasta_docstring[] =
        "Reads the FASTA file and formats it (checking duplicate headers, ambiguity characters, etc.) for TreeSAPP";

static PyMethodDef module_methods[] = {
        {"_read_format_fasta",
        read_format_fasta,
        METH_VARARGS,
        read_format_fasta_docstring},
        {NULL, NULL, 0, NULL}
};

class Fasta {
protected:
    int record_header( std::string );
    int record_sequence();

public:
    // Initialization functions
    Fasta( char * , char *, char *);
    Fasta( void );
    ~Fasta();
    Fasta( const Fasta& other);
    void clear( void );

    // Class objects
    std::set<std::string> header_base;
    PyObject *fasta_list;
    std::string sequence_buffer;
    std::string molecule;
    long int N_contigs;
    long int count_ambiguity;
    long int count_total;
    long int substitutions;
    char* log_file;
    char* write_buffer;
    ifstream *fasta_file;
    ofstream *parse_log;

    // Class functions
    int parse_fasta(int min_length);
//    int find_longest_contig();
//    int writeNx(std::string output, bool verbose);
};

#endif

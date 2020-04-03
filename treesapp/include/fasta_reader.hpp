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

static char read_format_fasta_docstring[] =
        "Reads the FASTA file and formats it (checking duplicate headers, ambiguity characters, etc.) for TreeSAPP";

static PyMethodDef module_methods[] = {
        {"_read_format_fasta",
        read_format_fasta,
        METH_VARARGS,
        read_format_fasta_docstring},
        {NULL, NULL, 0, NULL}
};

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

#if PY_MAJOR_VERSION >= 3

    static int module_traverse(PyObject *m, visitproc visit, void *arg) {
        Py_VISIT(GETSTATE(m)->error);
        return 0;
    }

    static int module_clear(PyObject *m) {
        Py_CLEAR(GETSTATE(m)->error);
        return 0;
    }

    static struct PyModuleDef module_def = {
        PyModuleDef_HEAD_INIT,
        "_fasta_reader",
        "This module rapidly loads and formats sequences from a FASTA file using C from within TreeSAPP",
        sizeof(struct module_state),
        module_methods,    /* m_methods */
        NULL,                /* m_reload */
        module_traverse,                /* m_traverse */
        module_clear,                /* m_clear */
        NULL,                /* m_free */
    };

#define INITERROR return NULL

PyMODINIT_FUNC PyInit__fasta_reader(void)

#else
#define INITERROR return

PyMODINIT_FUNC init_fasta_reader(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *m = PyModule_Create(&module_def);
#else
    static char module_docstring[] =
        "This module rapidly loads and formats sequences from a FASTA file using C from within TreeSAPP";
    PyObject *m = Py_InitModule3("_fasta_reader", module_methods, module_docstring);
#endif

    if (m == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(m);

    st->error = PyErr_NewException("_fasta_reader.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(m);
        INITERROR;
    }

#if PY_MAJOR_VERSION >= 3
    return m;
#endif
}

class Fasta {
protected:
    int record_header(std::string);
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
    int parse_fasta(int min_length, std::size_t max_header_length);
//    int find_longest_contig();
//    int writeNx(std::string output, bool verbose);
};

#endif

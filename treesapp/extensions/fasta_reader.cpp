/*
 * Fasta.cpp Connor Morgan-Lang <c.morganlang@gmail.com>
 * Hallam Lab, Department of Microbiology and Immunology, UBC
 * ------------------------------------------
 * Last modified: 23 August 2017 (CML)
 * ------------------------------------------
 */

#include "fasta_reader.hpp"

Fasta::Fasta(char * input, char * output_dir, char * molecule_type) {
    fasta_list = PyList_New(0);
    sequence_buffer.clear();
    header_base.clear();
    N_contigs = 0;
    count_ambiguity = 0;
    count_total = 0;
    substitutions = 0;
    molecule.assign(molecule_type);
    log_file = new char[1000];
    write_buffer = new char[1000];
    sprintf(log_file, "%s/fasta_reader_log.txt", output_dir);

    fasta_file = new ifstream(input);
    parse_log = new ofstream(log_file);

    if ( !fasta_file->is_open() ) {
        cerr << "Unable to open '" << input << "' for reading. Exiting now!" << endl;
        exit(0);
    }

    if ( !parse_log->is_open() ) {
        cerr << "Unable to open '" << log_file << "' for writing. Exiting now!" << endl;
        exit(0);
    }
}

Fasta::Fasta(void) { clear(); }

void Fasta::clear(void) {
    sequence_buffer.clear();
    count_total = 0;
    N_contigs = 0;
}

Fasta::~Fasta() {
    fasta_file->close();
    parse_log->close();
    sequence_buffer.clear();
    header_base.clear();
    delete[] write_buffer;
    delete[] log_file;
}

char replace_operators(char it) {
    if (it == ' ' || it == '\t')
        it = '_';
    if (it == ',' || it == '|')
        it = '_';
    if (it == ';')
        it = '_';
    return it;
}

string erase_characters(string str, const char* targets) {
    string purified;
    bool match;
    int i;
    std::size_t k;
    int j = 0;
    // Get the size of targets
    while (targets[j])
        j++;
    for (k = 0; k<str.length(); ++k) {
        i = 0;
        match = false;
        while ( i < j ) {
            if (str[k] == targets[i])
                match = true;
            i++;
        }
        if (!match)
            purified += str[k];
    }

    return purified;
}

bool ascii_range(char c) {
    return !(c>=0 && c <128);
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

string remove_non_ascii_char(std::string str) {
    str.erase(std::remove_if(str.begin(),
                             str.end(),
                             ascii_range),
              str.end());
    return str;
}

int Fasta::record_header(string line) {
    /*
    * Function for inserting a header into the fasta_list
    * record_header additionally checks for duplicate headers and appends "_" + occurrence
    */
    std::pair<std::set<string>::iterator,bool> ret;
    string new_header;
    new_header = line;

    // Replace whitespace characters and other ASCII operators with underscores
    rtrim(new_header);
    new_header = remove_non_ascii_char(new_header);
    new_header = erase_characters(new_header, "()[]/\\\\'<>");
    transform(new_header.begin(), new_header.end(), new_header.begin(), replace_operators);

    // Must remove period if it is the last character - BLAST ignores it causing problems downstream
    if (new_header[new_header.length()-1] == '.') {
        string tmp;
        while (new_header[new_header.length()-1] == '.') {
            tmp = new_header.substr(0, new_header.length()-1);
            new_header = tmp;
        }
    }

    ret = header_base.insert(new_header);
    if (!ret.second) {
        string redundant_header;
        int dup_num = 1;
        char buffer [10];

        string num_buffer;
        sprintf(buffer, "%d", dup_num);
        redundant_header = new_header + "_" + num_buffer.assign(buffer);
        ret = header_base.insert(redundant_header);
        while (!ret.second) {
            dup_num++;
            sprintf(buffer, "%d", dup_num);
            redundant_header = new_header + "_" + num_buffer.assign(buffer);
            ret = header_base.insert(redundant_header);
        }
        new_header = redundant_header;
        sprintf(write_buffer, "WARNING: Duplicate header (%s) replaced with %s\n", line.c_str(), new_header.c_str());
        parse_log->write(write_buffer, 44+line.length()+new_header.length());
    }
    PyList_Append(fasta_list, Py_BuildValue("s", new_header.c_str()));

    return 0;
}

int Fasta::record_sequence() {
    /*
    * Convert all lowercase characters to uppercase (done)
    * Find all `X`s and `N`s in the sequence and add to count_ambiguity (done)
    * Warn about the number of ambiguity character substitutions made
    * Raise error when unknown character is encountered
    */
    int seq_ambiguity = 0;
    int replacements = 0;
    int pos = 0;
    char *purified_seq = (char*) malloc (sequence_buffer.length() + 10);
    string acceptable_characters;
    std::map<char, char> iupac_map;
        iupac_map['R'] = 'A';
        iupac_map['Y'] = 'C';
        iupac_map['S'] = 'G';
        iupac_map['W'] = 'A';
        iupac_map['K'] = 'G';
        iupac_map['M'] = 'A';
        iupac_map['B'] = 'C';
        iupac_map['D'] = 'A';
        iupac_map['H'] = 'A';
        iupac_map['V'] = 'A';
    std::size_t acceptable;

    transform(sequence_buffer.begin(), sequence_buffer.end(), sequence_buffer.begin(), ::toupper);

    if (molecule == "prot") {
        acceptable_characters.assign("ACDEFGHIKLMNPQRSTUVWYBXZ.-*");
    }
    else if (molecule == "dna" || molecule == "rrna") {
        acceptable_characters.assign("-ACGTU.");
    }
    else {
        cerr << "ERROR: Unrecognized molecule. Unable to determine which characters to use." << endl;
        free(purified_seq);
        return 5;
    }

    for (string::iterator it=sequence_buffer.begin(); it!=sequence_buffer.end(); ++it){
        acceptable = acceptable_characters.find(*it);
        if (acceptable!=std::string::npos) {
            purified_seq[pos] = *it;
        }
        // Count the number of ambiguity characters
        else if (*it == 'X' || *it == 'N') {
            purified_seq[pos] = *it;
            seq_ambiguity++;
        }
        // Replace ambiguity characters
        else if (molecule == "dna" && iupac_map.find(*it) != iupac_map.end()) {
            purified_seq[pos] = iupac_map.find(*it)->second;
            replacements++;
        }
        // Skip sequence if the character is not a standard amino acid or nucleotide character
        else if (acceptable==std::string::npos) {
            cerr << "\nWARNING: '" << *it << "' at position " << pos << " is not an accepted character! ";
            free(purified_seq);
            return 4;
        }
        pos++;
    }
    purified_seq[pos] = '\0';
    sequence_buffer.assign(purified_seq);
    free(purified_seq);

    // Remove from analysis if greater than 75% is ambiguous
    if ( seq_ambiguity >= (sequence_buffer.length()*0.75) )
        return 3;
    count_total += sequence_buffer.length();
    count_ambiguity += seq_ambiguity;
    N_contigs++;
    if (replacements > (sequence_buffer.length()*0.1)) {
        parse_log->write("WARNING: Ambiguity character replacements accounted for >10% of sequence", 72);
        return 2;
    }
    if ( (seq_ambiguity*2) > (signed)sequence_buffer.length() )
        return 1;
    else
        return 0;
}

int Fasta::parse_fasta(int min_length, std::size_t max_header_length) {
    string line;
    std::string header;
    int status;

    if (fasta_file->good()) {
        getline( *fasta_file, line);
    }
    else {
        fprintf(stderr, "The fasta file cannot be parsed!\n");
        exit(0);
    }
    while ( fasta_file->good() ) {
        if ( line.empty() )
            ;
        else {
            if (line[0] == '>') {
                // If the sequence meets the length threshold, evaluate it's composition
                if ( (signed)sequence_buffer.length() >= min_length ) {
                    status = record_sequence();
                    if (status == 4) {
                        cerr << "Skipping header:\n'" << header << "'" << endl;
                        sprintf(write_buffer, "WARNING: Skipping %s due to bad character.\n", header.c_str());
                        parse_log->write(write_buffer, 53+header.length());
                    }
                    else if (status == 3) {
                        sprintf(write_buffer, "WARNING: Skipping %s due to >75%% ambiguity characters.\n", header.c_str());
                        parse_log->write(write_buffer, 53+header.length());
                    }
                    else if (status == 2) {
                        sprintf(write_buffer, " %s\n", header.c_str());
                        parse_log->write(write_buffer, 2+header.length());
                        record_header(header);
                        PyList_Append(fasta_list, Py_BuildValue("s", sequence_buffer.c_str()));
                    }
                    else if (status < 2) {
                        record_header(header);
                        PyList_Append(fasta_list, Py_BuildValue("s", sequence_buffer.c_str()));
                    }
                    else
                        return 5;
                }
                header = line.substr(1, max_header_length - 1);
                sequence_buffer.clear();
            }
            else {
                rtrim(line);
                sequence_buffer.append(line);
            }
        }
        getline( *fasta_file, line );
    }
    rtrim(line);
    sequence_buffer.append(line);
    // Ensure the last sequence is appended to fasta_list
    if ( (signed)sequence_buffer.length() >= min_length ) {
        status = record_sequence();
        if (status == 4) {
            cerr << "Skipping " << header << endl;
            sprintf(write_buffer, "WARNING: Skipping %s due to bad character.\n", header.c_str());
            parse_log->write(write_buffer, 53+header.length());
        }
        else if (status == 3) {
            sprintf(write_buffer, "WARNING: Skipping %s due to >75%% ambiguity characters.\n", header.c_str());
            parse_log->write(write_buffer, 53+header.length());
        }
        else if (status == 2) {
            sprintf(write_buffer, " %s\n", header.c_str());
            parse_log->write(write_buffer, 2+header.length());
            record_header(header);
            PyList_Append(fasta_list, Py_BuildValue("s", sequence_buffer.c_str()));
        }
        else if (status < 2) {
            record_header(header);
            PyList_Append(fasta_list, Py_BuildValue("s", sequence_buffer.c_str()));
        }
        else
            return 5;
    }
    // header_base.size() is assumed to be equal to N_contigs is everything was parsed properly...
    if ( (signed)header_base.size() != N_contigs )
        return 1;
    // Exit if no sequences were found
    if (count_total == 0)
        return 2;
    // Exit the program if less than half of the characters are standard nucleotides or amino acids
    if (count_ambiguity == count_total)
        return 3;
    else
        return 0;
}

static PyObject *read_format_fasta(PyObject *self, PyObject *args) {
    char * fasta_file;
    char * output_dir;
    int min_length;
    std::size_t max_header_length;
    char * molecule;
    if (!PyArg_ParseTuple(args, "sissn", &fasta_file, &min_length, &output_dir, &molecule, &max_header_length)) {
        return NULL;
    }
    /*
    * Create a Fasta class object to store the PyList
    * Begin reading in the fasta file, replacing ambiguity characters on the fly
    * Check for duplicate headers by looking through a sorted list. If duplicate found, append _N and insert
    */
    Fasta fasta_object(fasta_file, output_dir, molecule);
    int return_status = fasta_object.parse_fasta(min_length, max_header_length);
    if (return_status > 0)
        fasta_object.fasta_list = PyList_New(0);
    if (return_status == 1)
        fprintf(stderr, "ERROR: The input was not parsed correctly (N_contigs differs from vector size)!\n");
    if (return_status == 2)
        fprintf(stderr, "ERROR: Your input file appears to be corrupted. No sequences were found!\n");
    if (return_status == 3)
        fprintf(stderr, "ERROR: The majority of sequence(s) are completely ambiguous (only X or N)!\n");

    return fasta_object.fasta_list;
}

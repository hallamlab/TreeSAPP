#ifndef DYNAMITEsequencedbHEADERFILE
#define DYNAMITEsequencedbHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "sequence.h"
#include "database.h"
#include "hscore.h"


#define SequenceDBLISTLENGTH 128

enum SequenceDBFormat {
  SEQ_DB_UNKNOWN = 32,
  SEQ_DB_FASTA };

/* Object FileSource
 *
 * Descrip: This object represents a single
 *        file source for a database. At
 *        the moment only multiple fasta
 *        files are catered for
 *
 *
 */
struct Wise2_FileSource {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * filename;     
    FILE * input;   /*  could be stdin!  */ 
    int format;  
    int type;    
    } ;  
/* FileSource defined */ 
#ifndef DYNAMITE_DEFINED_FileSource
typedef struct Wise2_FileSource Wise2_FileSource;
#define FileSource Wise2_FileSource
#define DYNAMITE_DEFINED_FileSource
#endif


/* Object SequenceDB
 *
 * Descrip: This is the basic Sequence database
 *        wrapper - it handles all the formats
 *        and the on-the-fly indexing.
 *
 *        Generally it wont be directly used by
 *        an algorithm, which will be using something
 *        specific to the sequence type produce complex 
 *        sequence type objects: ie you will be using
 *        proteindb or cdnadb which internally
 *        will be using this object
 *
 *
 */
struct Wise2_SequenceDB {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    FileSource ** fs;    
    int len;/* len for above fs  */ 
    int maxlen; /* maxlen for above fs */ 
    int current_source;  
    FILE * current_file;     
    int sequence_no;     
    int byte_position;   
    int has_warned_single;   
    int seq_start;   
    int seq_end;     
    } ;  
/* SequenceDB defined */ 
#ifndef DYNAMITE_DEFINED_SequenceDB
typedef struct Wise2_SequenceDB Wise2_SequenceDB;
#define SequenceDB Wise2_SequenceDB
#define DYNAMITE_DEFINED_SequenceDB
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  get_Sequence_from_SequenceDB(sdb,de)
 *
 * Descrip:    Quite a mindless function which retrieves sequences
 *             via indexes
 *
 *             Going to spend too much time in fopen if this is used
 *             too much
 *
 *
 * Arg:        sdb [UNKN ] Undocumented argument [SequenceDB *]
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_get_Sequence_from_SequenceDB(SequenceDB * sdb,DataEntry * de);
#define get_Sequence_from_SequenceDB Wise2_get_Sequence_from_SequenceDB


/* Function:  add_SequenceDB_info_DataEntry(sdb,de)
 *
 * Descrip:    A function which places data into dataentry so we can
 *             be guarenteed to retrieve it sometime.
 *
 *             It uses 0 and 1 points in the Data array.
 *
 *
 * Arg:        sdb [UNKN ] Undocumented argument [SequenceDB *]
 * Arg:         de [UNKN ] Undocumented argument [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SequenceDB_info_DataEntry(SequenceDB * sdb,DataEntry * de);
#define add_SequenceDB_info_DataEntry Wise2_add_SequenceDB_info_DataEntry


/* Function:  close_SequenceDB(last,sdb)
 *
 * Descrip:    top level function that closes the SequenceDB
 *             after the last sequence is read.
 *
 *
 * Arg:        last [WRITE] Sequence object to be freed  [Sequence *]
 * Arg:         sdb [READ ] database to be closed [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_close_SequenceDB(Sequence * last,SequenceDB * sdb);
#define close_SequenceDB Wise2_close_SequenceDB


/* Function:  init_SequenceDB(sdb,return_status)
 *
 * Descrip:    top level function that starts a database read on
 *             SequenceDB
 *
 *
 *
 * Arg:                  sdb [READ ] sequence database [SequenceDB *]
 * Arg:        return_status [WRITE] returns the database status as found in database.h [int *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_init_SequenceDB(SequenceDB * sdb,int * return_status);
#define init_SequenceDB Wise2_init_SequenceDB


/* Function:  reload_SequenceDB(last,sdb,return_status)
 *
 * Descrip:    top level function that reloads a sequence database
 *
 *
 *
 * Arg:                 last [WRITE] previous sequence to be used: will simply be freed at the moment [Sequence *]
 * Arg:                  sdb [UNKN ] sequence database [SequenceDB *]
 * Arg:        return_status [WRITE] returns the database status as found in database.h [int *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
Sequence * Wise2_reload_SequenceDB(Sequence * last,SequenceDB * sdb,int * return_status);
#define reload_SequenceDB Wise2_reload_SequenceDB


/* Function:  SequenceDB_from_FILE_and_format(input,format)
 *
 * Descrip:    makes a SequencDB from a straight file stream.
 *
 *             This means SequenceDB will *not* close it when
 *             the SequenceDB is closed.
 *
 *
 * Arg:         input [READ ] filestream [FILE *]
 * Arg:        format [UNKN ] format as defined by /word_to_format [int]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * Wise2_SequenceDB_from_FILE_and_format(FILE * input,int format);
#define SequenceDB_from_FILE_and_format Wise2_SequenceDB_from_FILE_and_format


/* Function:  single_fasta_SequenceDB(filename)
 *
 * Descrip:    pre-packed single fasta file db
 *
 *
 *
 * Arg:        filename [UNKN ] name of fastadb [char *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * Wise2_single_fasta_SequenceDB(char * filename);
#define single_fasta_SequenceDB Wise2_single_fasta_SequenceDB


/* Function:  read_SequenceDB_line(line,ifp)
 *
 * Descrip:    Reads a SequenceDB definition from
 *
 *             seqdb <name>
 *             <filename> <format> <type>
 *             ...
 *             endseqdb
 *
 *
 *
 * Arg:        line [UNKN ] starting line (seqdb line) [char *]
 * Arg:         ifp [UNKN ] file input [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * Wise2_read_SequenceDB_line(char * line,FILE * ifp);
#define read_SequenceDB_line Wise2_read_SequenceDB_line


/* Function:  word_to_format(word)
 *
 * Descrip:    converts char * to format for SequenceDB FileSources
 *
 *
 * Arg:        word [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_word_to_format(char * word);
#define word_to_format Wise2_word_to_format


/* Function:  hard_link_FileSource(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FileSource *]
 *
 * Return [UNKN ]  Undocumented return value [FileSource *]
 *
 */
FileSource * Wise2_hard_link_FileSource(FileSource * obj);
#define hard_link_FileSource Wise2_hard_link_FileSource


/* Function:  FileSource_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FileSource *]
 *
 */
FileSource * Wise2_FileSource_alloc(void);
#define FileSource_alloc Wise2_FileSource_alloc


/* Function:  free_FileSource(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FileSource *]
 *
 * Return [UNKN ]  Undocumented return value [FileSource *]
 *
 */
FileSource * Wise2_free_FileSource(FileSource * obj);
#define free_FileSource Wise2_free_FileSource


/* Function:  add_SequenceDB(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SequenceDB *]
 * Arg:        add [OWNER] Object to add to the list [FileSource *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SequenceDB(SequenceDB * obj,FileSource * add);
#define add_SequenceDB Wise2_add_SequenceDB


/* Function:  flush_SequenceDB(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SequenceDB(SequenceDB * obj);
#define flush_SequenceDB Wise2_flush_SequenceDB


/* Function:  SequenceDB_alloc_std(void)
 *
 * Descrip:    Equivalent to SequenceDB_alloc_len(SequenceDBLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * Wise2_SequenceDB_alloc_std(void);
#define SequenceDB_alloc_std Wise2_SequenceDB_alloc_std


/* Function:  SequenceDB_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * Wise2_SequenceDB_alloc_len(int len);
#define SequenceDB_alloc_len Wise2_SequenceDB_alloc_len


/* Function:  hard_link_SequenceDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * Wise2_hard_link_SequenceDB(SequenceDB * obj);
#define hard_link_SequenceDB Wise2_hard_link_SequenceDB


/* Function:  SequenceDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * Wise2_SequenceDB_alloc(void);
#define SequenceDB_alloc Wise2_SequenceDB_alloc


/* Function:  free_SequenceDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SequenceDB *]
 *
 * Return [UNKN ]  Undocumented return value [SequenceDB *]
 *
 */
SequenceDB * Wise2_free_SequenceDB(SequenceDB * obj);
#define free_SequenceDB Wise2_free_SequenceDB


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
FileSource * Wise2_access_fs_SequenceDB(SequenceDB * obj,int i);
#define access_fs_SequenceDB Wise2_access_fs_SequenceDB
int Wise2_access_current_source_SequenceDB(SequenceDB * obj);
#define access_current_source_SequenceDB Wise2_access_current_source_SequenceDB
boolean Wise2_replace_current_source_SequenceDB(SequenceDB * obj,int current_source);
#define replace_current_source_SequenceDB Wise2_replace_current_source_SequenceDB
boolean Wise2_replace_current_file_SequenceDB(SequenceDB * obj,FILE * current_file);
#define replace_current_file_SequenceDB Wise2_replace_current_file_SequenceDB
boolean Wise2_replace_sequence_no_SequenceDB(SequenceDB * obj,int sequence_no);
#define replace_sequence_no_SequenceDB Wise2_replace_sequence_no_SequenceDB
int Wise2_access_sequence_no_SequenceDB(SequenceDB * obj);
#define access_sequence_no_SequenceDB Wise2_access_sequence_no_SequenceDB
boolean Wise2_replace_name_SequenceDB(SequenceDB * obj,char * name);
#define replace_name_SequenceDB Wise2_replace_name_SequenceDB
int Wise2_access_format_FileSource(FileSource * obj);
#define access_format_FileSource Wise2_access_format_FileSource
boolean Wise2_replace_byte_position_SequenceDB(SequenceDB * obj,int byte_position);
#define replace_byte_position_SequenceDB Wise2_replace_byte_position_SequenceDB
int Wise2_length_fs_SequenceDB(SequenceDB * obj);
#define length_fs_SequenceDB Wise2_length_fs_SequenceDB
boolean Wise2_replace_format_FileSource(FileSource * obj,int format);
#define replace_format_FileSource Wise2_replace_format_FileSource
int Wise2_access_byte_position_SequenceDB(SequenceDB * obj);
#define access_byte_position_SequenceDB Wise2_access_byte_position_SequenceDB
FILE * Wise2_access_current_file_SequenceDB(SequenceDB * obj);
#define access_current_file_SequenceDB Wise2_access_current_file_SequenceDB
boolean Wise2_replace_filename_FileSource(FileSource * obj,char * filename);
#define replace_filename_FileSource Wise2_replace_filename_FileSource
char * Wise2_access_name_SequenceDB(SequenceDB * obj);
#define access_name_SequenceDB Wise2_access_name_SequenceDB
char * Wise2_access_filename_FileSource(FileSource * obj);
#define access_filename_FileSource Wise2_access_filename_FileSource
boolean Wise2_replace_type_FileSource(FileSource * obj,int type);
#define replace_type_FileSource Wise2_replace_type_FileSource
boolean Wise2_replace_input_FileSource(FileSource * obj,FILE * input);
#define replace_input_FileSource Wise2_replace_input_FileSource
int Wise2_access_type_FileSource(FileSource * obj);
#define access_type_FileSource Wise2_access_type_FileSource
FILE * Wise2_access_input_FileSource(FileSource * obj);
#define access_input_FileSource Wise2_access_input_FileSource
Sequence * Wise2_get_next_SequenceDB(SequenceDB * sdb);
#define get_next_SequenceDB Wise2_get_next_SequenceDB
boolean Wise2_SequenceDB_at_end(SequenceDB * sdb);
#define SequenceDB_at_end Wise2_SequenceDB_at_end
boolean Wise2_load_next_fs_SequenceDB(SequenceDB * sdb);
#define load_next_fs_SequenceDB Wise2_load_next_fs_SequenceDB
boolean Wise2_close_last_fs_SequenceDB(SequenceDB * sdb);
#define close_last_fs_SequenceDB Wise2_close_last_fs_SequenceDB
FileSource * Wise2_FileSource_from_FILE_and_format(FILE * input,int format);
#define FileSource_from_FILE_and_format Wise2_FileSource_from_FILE_and_format
FileSource * Wise2_FileSource_from_line(char * line);
#define FileSource_from_line Wise2_FileSource_from_line
void Wise2_swap_SequenceDB(FileSource ** list,int i,int j) ;
#define swap_SequenceDB Wise2_swap_SequenceDB
void Wise2_qsort_SequenceDB(FileSource ** list,int left,int right,int (*comp)(FileSource * ,FileSource * ));
#define qsort_SequenceDB Wise2_qsort_SequenceDB
void Wise2_sort_SequenceDB(SequenceDB * obj,int (*comp)(FileSource *, FileSource *));
#define sort_SequenceDB Wise2_sort_SequenceDB
boolean Wise2_expand_SequenceDB(SequenceDB * obj,int len);
#define expand_SequenceDB Wise2_expand_SequenceDB

#ifdef _cplusplus
}
#endif

#endif

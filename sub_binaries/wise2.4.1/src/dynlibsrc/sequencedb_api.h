

/* Functions that create, manipulate or act on SequenceDB
 *
 * Wise2_close_SequenceDB
 * Wise2_hard_link_SequenceDB
 * Wise2_SequenceDB_alloc_std
 * Wise2_replace_name_SequenceDB
 * Wise2_access_name_SequenceDB
 * Wise2_access_fs_SequenceDB
 * Wise2_length_fs_SequenceDB
 * Wise2_flush_SequenceDB
 * Wise2_add_SequenceDB
 * Wise2_replace_current_source_SequenceDB
 * Wise2_access_current_source_SequenceDB
 * Wise2_replace_current_file_SequenceDB
 * Wise2_access_current_file_SequenceDB
 * Wise2_replace_sequence_no_SequenceDB
 * Wise2_access_sequence_no_SequenceDB
 * Wise2_replace_byte_position_SequenceDB
 * Wise2_access_byte_position_SequenceDB
 * Wise2_free_SequenceDB [destructor]
 *
 */



/* Functions that create, manipulate or act on FileSource
 *
 * Wise2_hard_link_FileSource
 * Wise2_FileSource_alloc
 * Wise2_replace_filename_FileSource
 * Wise2_access_filename_FileSource
 * Wise2_replace_input_FileSource
 * Wise2_access_input_FileSource
 * Wise2_replace_format_FileSource
 * Wise2_access_format_FileSource
 * Wise2_replace_type_FileSource
 * Wise2_access_type_FileSource
 * Wise2_free_FileSource [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_single_fasta_SequenceDB
 *

/* API for object SequenceDB */
/* Function:  Wise2_close_SequenceDB(last,sdb)
 *
 * Descrip:    top level function that closes the SequenceDB
 *             after the last sequence is read.
 *
 *
 * Arg:        last         Sequence object to be freed  [Wise2_Sequence *]
 * Arg:        sdb          database to be closed [Wise2_SequenceDB *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_close_SequenceDB( Wise2_Sequence * last,Wise2_SequenceDB * sdb);

/* Function:  Wise2_hard_link_SequenceDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_SequenceDB *]
 *
 * Returns Undocumented return value [Wise2_SequenceDB *]
 *
 */
Wise2_SequenceDB * Wise2_hard_link_SequenceDB( Wise2_SequenceDB * obj);

/* Function:  Wise2_SequenceDB_alloc_std(void)
 *
 * Descrip:    Equivalent to SequenceDB_alloc_len(SequenceDBLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_SequenceDB *]
 *
 */
Wise2_SequenceDB * Wise2_SequenceDB_alloc_std();

/* Function:  Wise2_replace_name_SequenceDB(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SequenceDB *]
 * Arg:        name         New value of the variable [char *]
 *
 * Returns member variable name [boolean]
 *
 */
boolean Wise2_replace_name_SequenceDB( Wise2_SequenceDB * obj,char * name);

/* Function:  Wise2_access_name_SequenceDB(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SequenceDB *]
 *
 * Returns member variable name [char *]
 *
 */
char * Wise2_access_name_SequenceDB( Wise2_SequenceDB * obj);

/* Function:  Wise2_access_fs_SequenceDB(obj,i)
 *
 * Descrip:    Access members stored in the fs list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_SequenceDB *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_FileSource *]
 *
 */
Wise2_FileSource * Wise2_access_fs_SequenceDB( Wise2_SequenceDB * obj,int i);

/* Function:  Wise2_length_fs_SequenceDB(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_SequenceDB *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_fs_SequenceDB( Wise2_SequenceDB * obj);

/* Function:  Wise2_flush_SequenceDB(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_SequenceDB *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_SequenceDB( Wise2_SequenceDB * obj);

/* Function:  Wise2_add_SequenceDB(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_SequenceDB *]
 * Arg:        add          Object to add to the list [Wise2_FileSource *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SequenceDB( Wise2_SequenceDB * obj,Wise2_FileSource * add);

/* Function:  Wise2_replace_current_source_SequenceDB(obj,current_source)
 *
 * Descrip:    Replace member variable current_source
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SequenceDB *]
 * Arg:        current_source New value of the variable [int]
 *
 * Returns member variable current_source [boolean]
 *
 */
boolean Wise2_replace_current_source_SequenceDB( Wise2_SequenceDB * obj,int current_source);

/* Function:  Wise2_access_current_source_SequenceDB(obj)
 *
 * Descrip:    Access member variable current_source
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SequenceDB *]
 *
 * Returns member variable current_source [int]
 *
 */
int Wise2_access_current_source_SequenceDB( Wise2_SequenceDB * obj);

/* Function:  Wise2_replace_current_file_SequenceDB(obj,current_file)
 *
 * Descrip:    Replace member variable current_file
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SequenceDB *]
 * Arg:        current_file New value of the variable [FILE *]
 *
 * Returns member variable current_file [boolean]
 *
 */
boolean Wise2_replace_current_file_SequenceDB( Wise2_SequenceDB * obj,FILE * current_file);

/* Function:  Wise2_access_current_file_SequenceDB(obj)
 *
 * Descrip:    Access member variable current_file
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SequenceDB *]
 *
 * Returns member variable current_file [FILE *]
 *
 */
FILE * Wise2_access_current_file_SequenceDB( Wise2_SequenceDB * obj);

/* Function:  Wise2_replace_sequence_no_SequenceDB(obj,sequence_no)
 *
 * Descrip:    Replace member variable sequence_no
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SequenceDB *]
 * Arg:        sequence_no  New value of the variable [int]
 *
 * Returns member variable sequence_no [boolean]
 *
 */
boolean Wise2_replace_sequence_no_SequenceDB( Wise2_SequenceDB * obj,int sequence_no);

/* Function:  Wise2_access_sequence_no_SequenceDB(obj)
 *
 * Descrip:    Access member variable sequence_no
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SequenceDB *]
 *
 * Returns member variable sequence_no [int]
 *
 */
int Wise2_access_sequence_no_SequenceDB( Wise2_SequenceDB * obj);

/* Function:  Wise2_replace_byte_position_SequenceDB(obj,byte_position)
 *
 * Descrip:    Replace member variable byte_position
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SequenceDB *]
 * Arg:        byte_position New value of the variable [int]
 *
 * Returns member variable byte_position [boolean]
 *
 */
boolean Wise2_replace_byte_position_SequenceDB( Wise2_SequenceDB * obj,int byte_position);

/* Function:  Wise2_access_byte_position_SequenceDB(obj)
 *
 * Descrip:    Access member variable byte_position
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SequenceDB *]
 *
 * Returns member variable byte_position [int]
 *
 */
int Wise2_access_byte_position_SequenceDB( Wise2_SequenceDB * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_SequenceDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_SequenceDB *]
 *
 * Returns Undocumented return value [Wise2_SequenceDB *]
 *
 */
Wise2_SequenceDB * Wise2_free_SequenceDB( Wise2_SequenceDB * obj);

/* API for object FileSource */
/* Function:  Wise2_hard_link_FileSource(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_FileSource *]
 *
 * Returns Undocumented return value [Wise2_FileSource *]
 *
 */
Wise2_FileSource * Wise2_hard_link_FileSource( Wise2_FileSource * obj);

/* Function:  Wise2_FileSource_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_FileSource *]
 *
 */
Wise2_FileSource * Wise2_FileSource_alloc();

/* Function:  Wise2_replace_filename_FileSource(obj,filename)
 *
 * Descrip:    Replace member variable filename
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_FileSource *]
 * Arg:        filename     New value of the variable [char *]
 *
 * Returns member variable filename [boolean]
 *
 */
boolean Wise2_replace_filename_FileSource( Wise2_FileSource * obj,char * filename);

/* Function:  Wise2_access_filename_FileSource(obj)
 *
 * Descrip:    Access member variable filename
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_FileSource *]
 *
 * Returns member variable filename [char *]
 *
 */
char * Wise2_access_filename_FileSource( Wise2_FileSource * obj);

/* Function:  Wise2_replace_input_FileSource(obj,input)
 *
 * Descrip:    Replace member variable input
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_FileSource *]
 * Arg:        input        New value of the variable [FILE *]
 *
 * Returns member variable input [boolean]
 *
 */
boolean Wise2_replace_input_FileSource( Wise2_FileSource * obj,FILE * input);

/* Function:  Wise2_access_input_FileSource(obj)
 *
 * Descrip:    Access member variable input
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_FileSource *]
 *
 * Returns member variable input [FILE *]
 *
 */
FILE * Wise2_access_input_FileSource( Wise2_FileSource * obj);

/* Function:  Wise2_replace_format_FileSource(obj,format)
 *
 * Descrip:    Replace member variable format
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_FileSource *]
 * Arg:        format       New value of the variable [int]
 *
 * Returns member variable format [boolean]
 *
 */
boolean Wise2_replace_format_FileSource( Wise2_FileSource * obj,int format);

/* Function:  Wise2_access_format_FileSource(obj)
 *
 * Descrip:    Access member variable format
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_FileSource *]
 *
 * Returns member variable format [int]
 *
 */
int Wise2_access_format_FileSource( Wise2_FileSource * obj);

/* Function:  Wise2_replace_type_FileSource(obj,type)
 *
 * Descrip:    Replace member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_FileSource *]
 * Arg:        type         New value of the variable [int]
 *
 * Returns member variable type [boolean]
 *
 */
boolean Wise2_replace_type_FileSource( Wise2_FileSource * obj,int type);

/* Function:  Wise2_access_type_FileSource(obj)
 *
 * Descrip:    Access member variable type
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_FileSource *]
 *
 * Returns member variable type [int]
 *
 */
int Wise2_access_type_FileSource( Wise2_FileSource * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_FileSource(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_FileSource *]
 *
 * Returns Undocumented return value [Wise2_FileSource *]
 *
 */
Wise2_FileSource * Wise2_free_FileSource( Wise2_FileSource * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_single_fasta_SequenceDB(filename)
 *
 * Descrip:    pre-packed single fasta file db
 *
 *
 *
 * Arg:        filename     name of fastadb [char *]
 *
 * Returns Undocumented return value [Wise2_SequenceDB *]
 *
 */
Wise2_SequenceDB * Wise2_single_fasta_SequenceDB( char * filename);


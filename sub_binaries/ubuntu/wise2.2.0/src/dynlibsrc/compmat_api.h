

/* Functions that create, manipulate or act on CompMat
 *
 * Wise2_fail_safe_CompMat_access
 * Wise2_write_Blast_CompMat
 * Wise2_read_Blast_file_CompMat
 * Wise2_read_Blast_CompMat
 * Wise2_hard_link_CompMat
 * Wise2_CompMat_alloc
 * Wise2_replace_name_CompMat
 * Wise2_access_name_CompMat
 * Wise2_free_CompMat [destructor]
 *
 */

/* API for object CompMat */
/* Function:  Wise2_fail_safe_CompMat_access(cm,aa1,aa2)
 *
 * Descrip:    gives the fail form of the macro CompMat_AAMATCH which 
 *             checks that aa1 and a2 are sensible and that cm is not NULL.
 *
 *
 * Arg:        cm           compmat object [Wise2_CompMat *]
 * Arg:        aa1          first amino acid [int]
 * Arg:        aa2          second amino acid [int]
 *
 * Returns Undocumented return value [Score]
 *
 */
Score Wise2_fail_safe_CompMat_access( Wise2_CompMat * cm,int aa1,int aa2);

/* Function:  Wise2_write_Blast_CompMat(cm,ofp)
 *
 * Descrip:    writes a protien CompMat with a standard
 *             alphabet.
 *
 *
 * Arg:        cm           CompMat object [Wise2_CompMat *]
 * Arg:        ofp          file to output [FILE *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_write_Blast_CompMat( Wise2_CompMat * cm,FILE * ofp);

/* Function:  Wise2_read_Blast_file_CompMat(filename)
 *
 * Descrip:    Opens file, reads matrix, closes file.
 *             calls /read_Blast_CompMat for the actual format
 *             reading. Uses /openfile to open the file,
 *             so will open from config files.
 *
 *
 * Arg:        filename     Undocumented argument [char *]
 *
 * Returns Undocumented return value [Wise2_CompMat *]
 *
 */
Wise2_CompMat * Wise2_read_Blast_file_CompMat( char * filename);

/* Function:  Wise2_read_Blast_CompMat(ifp)
 *
 * Descrip:    reads a BLAST format matrix and
 *             allocates a new ComMat structure.
 *
 *
 * Arg:        ifp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [Wise2_CompMat *]
 *
 */
Wise2_CompMat * Wise2_read_Blast_CompMat( FILE * ifp);

/* Function:  Wise2_hard_link_CompMat(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_CompMat *]
 *
 * Returns Undocumented return value [Wise2_CompMat *]
 *
 */
Wise2_CompMat * Wise2_hard_link_CompMat( Wise2_CompMat * obj);

/* Function:  Wise2_CompMat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_CompMat *]
 *
 */
Wise2_CompMat * Wise2_CompMat_alloc();

/* Function:  Wise2_replace_name_CompMat(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_CompMat *]
 * Arg:        name         New value of the variable [char *]
 *
 * Returns member variable name [boolean]
 *
 */
boolean Wise2_replace_name_CompMat( Wise2_CompMat * obj,char * name);

/* Function:  Wise2_access_name_CompMat(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_CompMat *]
 *
 * Returns member variable name [char *]
 *
 */
char * Wise2_access_name_CompMat( Wise2_CompMat * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_CompMat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_CompMat *]
 *
 * Returns Undocumented return value [Wise2_CompMat *]
 *
 */
Wise2_CompMat * Wise2_free_CompMat( Wise2_CompMat * obj);


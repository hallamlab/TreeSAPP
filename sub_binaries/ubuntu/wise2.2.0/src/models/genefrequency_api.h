

/* Functions that create, manipulate or act on GeneFrequency21
 *
 * Wise2_hard_link_GeneFrequency21
 * Wise2_GeneFrequency21_alloc
 * Wise2_replace_ss5_GeneFrequency21
 * Wise2_access_ss5_GeneFrequency21
 * Wise2_replace_ss3_GeneFrequency21
 * Wise2_access_ss3_GeneFrequency21
 * Wise2_free_GeneFrequency21 [destructor]
 *
 */



/* Functions that create, manipulate or act on GeneConsensus
 *
 * Wise2_hard_link_GeneConsensus
 * Wise2_GeneConsensus_alloc_std
 * Wise2_replace_center_GeneConsensus
 * Wise2_access_center_GeneConsensus
 * Wise2_access_gsc_GeneConsensus
 * Wise2_length_gsc_GeneConsensus
 * Wise2_flush_GeneConsensus
 * Wise2_add_GeneConsensus
 * Wise2_free_GeneConsensus [destructor]
 *
 */



/* Functions that create, manipulate or act on GeneSingleCons
 *
 * Wise2_hard_link_GeneSingleCons
 * Wise2_GeneSingleCons_alloc
 * Wise2_replace_string_GeneSingleCons
 * Wise2_access_string_GeneSingleCons
 * Wise2_replace_number_GeneSingleCons
 * Wise2_access_number_GeneSingleCons
 * Wise2_free_GeneSingleCons [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_read_GeneFrequency21_file
 * Wise2_read_GeneFrequency21
 *

/* API for object GeneFrequency21 */
/* Function:  Wise2_hard_link_GeneFrequency21(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_GeneFrequency21 *]
 *
 * Returns Undocumented return value [Wise2_GeneFrequency21 *]
 *
 */
Wise2_GeneFrequency21 * Wise2_hard_link_GeneFrequency21( Wise2_GeneFrequency21 * obj);

/* Function:  Wise2_GeneFrequency21_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_GeneFrequency21 *]
 *
 */
Wise2_GeneFrequency21 * Wise2_GeneFrequency21_alloc();

/* Function:  Wise2_replace_ss5_GeneFrequency21(obj,ss5)
 *
 * Descrip:    Replace member variable ss5
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GeneFrequency21 *]
 * Arg:        ss5          New value of the variable [Wise2_GeneConsensus *]
 *
 * Returns member variable ss5 [boolean]
 *
 */
boolean Wise2_replace_ss5_GeneFrequency21( Wise2_GeneFrequency21 * obj,Wise2_GeneConsensus * ss5);

/* Function:  Wise2_access_ss5_GeneFrequency21(obj)
 *
 * Descrip:    Access member variable ss5
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GeneFrequency21 *]
 *
 * Returns member variable ss5 [Wise2_GeneConsensus *]
 *
 */
Wise2_GeneConsensus * Wise2_access_ss5_GeneFrequency21( Wise2_GeneFrequency21 * obj);

/* Function:  Wise2_replace_ss3_GeneFrequency21(obj,ss3)
 *
 * Descrip:    Replace member variable ss3
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GeneFrequency21 *]
 * Arg:        ss3          New value of the variable [Wise2_GeneConsensus *]
 *
 * Returns member variable ss3 [boolean]
 *
 */
boolean Wise2_replace_ss3_GeneFrequency21( Wise2_GeneFrequency21 * obj,Wise2_GeneConsensus * ss3);

/* Function:  Wise2_access_ss3_GeneFrequency21(obj)
 *
 * Descrip:    Access member variable ss3
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GeneFrequency21 *]
 *
 * Returns member variable ss3 [Wise2_GeneConsensus *]
 *
 */
Wise2_GeneConsensus * Wise2_access_ss3_GeneFrequency21( Wise2_GeneFrequency21 * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_GeneFrequency21(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_GeneFrequency21 *]
 *
 * Returns Undocumented return value [Wise2_GeneFrequency21 *]
 *
 */
Wise2_GeneFrequency21 * Wise2_free_GeneFrequency21( Wise2_GeneFrequency21 * obj);

/* API for object GeneConsensus */
/* Function:  Wise2_hard_link_GeneConsensus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_GeneConsensus *]
 *
 * Returns Undocumented return value [Wise2_GeneConsensus *]
 *
 */
Wise2_GeneConsensus * Wise2_hard_link_GeneConsensus( Wise2_GeneConsensus * obj);

/* Function:  Wise2_GeneConsensus_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneConsensus_alloc_len(GeneConsensusLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_GeneConsensus *]
 *
 */
Wise2_GeneConsensus * Wise2_GeneConsensus_alloc_std();

/* Function:  Wise2_replace_center_GeneConsensus(obj,center)
 *
 * Descrip:    Replace member variable center
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GeneConsensus *]
 * Arg:        center       New value of the variable [int]
 *
 * Returns member variable center [boolean]
 *
 */
boolean Wise2_replace_center_GeneConsensus( Wise2_GeneConsensus * obj,int center);

/* Function:  Wise2_access_center_GeneConsensus(obj)
 *
 * Descrip:    Access member variable center
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GeneConsensus *]
 *
 * Returns member variable center [int]
 *
 */
int Wise2_access_center_GeneConsensus( Wise2_GeneConsensus * obj);

/* Function:  Wise2_access_gsc_GeneConsensus(obj,i)
 *
 * Descrip:    Access members stored in the gsc list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_GeneConsensus *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_GeneSingleCons *]
 *
 */
Wise2_GeneSingleCons * Wise2_access_gsc_GeneConsensus( Wise2_GeneConsensus * obj,int i);

/* Function:  Wise2_length_gsc_GeneConsensus(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_GeneConsensus *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_gsc_GeneConsensus( Wise2_GeneConsensus * obj);

/* Function:  Wise2_flush_GeneConsensus(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_GeneConsensus *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_GeneConsensus( Wise2_GeneConsensus * obj);

/* Function:  Wise2_add_GeneConsensus(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_GeneConsensus *]
 * Arg:        add          Object to add to the list [Wise2_GeneSingleCons *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GeneConsensus( Wise2_GeneConsensus * obj,Wise2_GeneSingleCons * add);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_GeneConsensus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_GeneConsensus *]
 *
 * Returns Undocumented return value [Wise2_GeneConsensus *]
 *
 */
Wise2_GeneConsensus * Wise2_free_GeneConsensus( Wise2_GeneConsensus * obj);

/* API for object GeneSingleCons */
/* Function:  Wise2_hard_link_GeneSingleCons(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_GeneSingleCons *]
 *
 * Returns Undocumented return value [Wise2_GeneSingleCons *]
 *
 */
Wise2_GeneSingleCons * Wise2_hard_link_GeneSingleCons( Wise2_GeneSingleCons * obj);

/* Function:  Wise2_GeneSingleCons_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_GeneSingleCons *]
 *
 */
Wise2_GeneSingleCons * Wise2_GeneSingleCons_alloc();

/* Function:  Wise2_replace_string_GeneSingleCons(obj,string)
 *
 * Descrip:    Replace member variable string
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GeneSingleCons *]
 * Arg:        string       New value of the variable [char *]
 *
 * Returns member variable string [boolean]
 *
 */
boolean Wise2_replace_string_GeneSingleCons( Wise2_GeneSingleCons * obj,char * string);

/* Function:  Wise2_access_string_GeneSingleCons(obj)
 *
 * Descrip:    Access member variable string
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GeneSingleCons *]
 *
 * Returns member variable string [char *]
 *
 */
char * Wise2_access_string_GeneSingleCons( Wise2_GeneSingleCons * obj);

/* Function:  Wise2_replace_number_GeneSingleCons(obj,number)
 *
 * Descrip:    Replace member variable number
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GeneSingleCons *]
 * Arg:        number       New value of the variable [double]
 *
 * Returns member variable number [boolean]
 *
 */
boolean Wise2_replace_number_GeneSingleCons( Wise2_GeneSingleCons * obj,double number);

/* Function:  Wise2_access_number_GeneSingleCons(obj)
 *
 * Descrip:    Access member variable number
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_GeneSingleCons *]
 *
 * Returns member variable number [double]
 *
 */
double Wise2_access_number_GeneSingleCons( Wise2_GeneSingleCons * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_GeneSingleCons(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_GeneSingleCons *]
 *
 * Returns Undocumented return value [Wise2_GeneSingleCons *]
 *
 */
Wise2_GeneSingleCons * Wise2_free_GeneSingleCons( Wise2_GeneSingleCons * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_read_GeneFrequency21_file(filename)
 *
 * Descrip:    Opens the file with /openfile
 *
 *             Reads in a GeneFrequency (Mor-Ewan style)
 *
 *
 *
 * Arg:        filename     will open from WISECONFIGDIR etc via openfile [char *]
 *
 * Returns a newly allocated structure [Wise2_GeneFrequency21 *]
 *
 */
Wise2_GeneFrequency21 * Wise2_read_GeneFrequency21_file( char * filename);

/* Function:  Wise2_read_GeneFrequency21(ifp)
 *
 * Descrip:    Reads in a GeneFrequency (Mor-Ewan style)
 *             file from ifp
 *
 *
 * Arg:        ifp          file pointer [FILE *]
 *
 * Returns a newly allocated structure [Wise2_GeneFrequency21 *]
 *
 */
Wise2_GeneFrequency21 * Wise2_read_GeneFrequency21( FILE * ifp);


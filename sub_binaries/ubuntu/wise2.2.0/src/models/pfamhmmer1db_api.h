

/* Functions that create, manipulate or act on PfamHmmer1DB
 *
 * Wise2_hard_link_PfamHmmer1DB
 * Wise2_PfamHmmer1DB_alloc_std
 * Wise2_access_en_PfamHmmer1DB
 * Wise2_length_en_PfamHmmer1DB
 * Wise2_flush_PfamHmmer1DB
 * Wise2_add_PfamHmmer1DB
 * Wise2_replace_dirname_PfamHmmer1DB
 * Wise2_access_dirname_PfamHmmer1DB
 * Wise2_replace_cur_PfamHmmer1DB
 * Wise2_access_cur_PfamHmmer1DB
 * Wise2_replace_def_PfamHmmer1DB
 * Wise2_access_def_PfamHmmer1DB
 * Wise2_free_PfamHmmer1DB [destructor]
 *
 */



/* Functions that create, manipulate or act on PfamHmmer1Entry
 *
 * Wise2_hard_link_PfamHmmer1Entry
 * Wise2_PfamHmmer1Entry_alloc
 * Wise2_replace_entryname_PfamHmmer1Entry
 * Wise2_access_entryname_PfamHmmer1Entry
 * Wise2_replace_is_random_PfamHmmer1Entry
 * Wise2_access_is_random_PfamHmmer1Entry
 * Wise2_replace_is_hmmls_PfamHmmer1Entry
 * Wise2_access_is_hmmls_PfamHmmer1Entry
 * Wise2_replace_bits_cutoff_PfamHmmer1Entry
 * Wise2_access_bits_cutoff_PfamHmmer1Entry
 * Wise2_free_PfamHmmer1Entry [destructor]
 *
 */

/* API for object PfamHmmer1DB */
/* Function:  Wise2_hard_link_PfamHmmer1DB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_PfamHmmer1DB *]
 *
 * Returns Undocumented return value [Wise2_PfamHmmer1DB *]
 *
 */
Wise2_PfamHmmer1DB * Wise2_hard_link_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj);

/* Function:  Wise2_PfamHmmer1DB_alloc_std(void)
 *
 * Descrip:    Equivalent to PfamHmmer1DB_alloc_len(PfamHmmer1DBLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_PfamHmmer1DB *]
 *
 */
Wise2_PfamHmmer1DB * Wise2_PfamHmmer1DB_alloc_std();

/* Function:  Wise2_access_en_PfamHmmer1DB(obj,i)
 *
 * Descrip:    Access members stored in the en list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_PfamHmmer1DB *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_PfamHmmer1Entry *]
 *
 */
Wise2_PfamHmmer1Entry * Wise2_access_en_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj,int i);

/* Function:  Wise2_length_en_PfamHmmer1DB(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_PfamHmmer1DB *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_en_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj);

/* Function:  Wise2_flush_PfamHmmer1DB(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_PfamHmmer1DB *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj);

/* Function:  Wise2_add_PfamHmmer1DB(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_PfamHmmer1DB *]
 * Arg:        add          Object to add to the list [Wise2_PfamHmmer1Entry *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj,Wise2_PfamHmmer1Entry * add);

/* Function:  Wise2_replace_dirname_PfamHmmer1DB(obj,dirname)
 *
 * Descrip:    Replace member variable dirname
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1DB *]
 * Arg:        dirname      New value of the variable [char *]
 *
 * Returns member variable dirname [boolean]
 *
 */
boolean Wise2_replace_dirname_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj,char * dirname);

/* Function:  Wise2_access_dirname_PfamHmmer1DB(obj)
 *
 * Descrip:    Access member variable dirname
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1DB *]
 *
 * Returns member variable dirname [char *]
 *
 */
char * Wise2_access_dirname_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj);

/* Function:  Wise2_replace_cur_PfamHmmer1DB(obj,cur)
 *
 * Descrip:    Replace member variable cur
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1DB *]
 * Arg:        cur          New value of the variable [int]
 *
 * Returns member variable cur [boolean]
 *
 */
boolean Wise2_replace_cur_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj,int cur);

/* Function:  Wise2_access_cur_PfamHmmer1DB(obj)
 *
 * Descrip:    Access member variable cur
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1DB *]
 *
 * Returns member variable cur [int]
 *
 */
int Wise2_access_cur_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj);

/* Function:  Wise2_replace_def_PfamHmmer1DB(obj,def)
 *
 * Descrip:    Replace member variable def
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1DB *]
 * Arg:        def          New value of the variable [Wise2_RandomModel *]
 *
 * Returns member variable def [boolean]
 *
 */
boolean Wise2_replace_def_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj,Wise2_RandomModel * def);

/* Function:  Wise2_access_def_PfamHmmer1DB(obj)
 *
 * Descrip:    Access member variable def
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1DB *]
 *
 * Returns member variable def [Wise2_RandomModel *]
 *
 */
Wise2_RandomModel * Wise2_access_def_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_PfamHmmer1DB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_PfamHmmer1DB *]
 *
 * Returns Undocumented return value [Wise2_PfamHmmer1DB *]
 *
 */
Wise2_PfamHmmer1DB * Wise2_free_PfamHmmer1DB( Wise2_PfamHmmer1DB * obj);

/* API for object PfamHmmer1Entry */
/* Function:  Wise2_hard_link_PfamHmmer1Entry(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_PfamHmmer1Entry *]
 *
 * Returns Undocumented return value [Wise2_PfamHmmer1Entry *]
 *
 */
Wise2_PfamHmmer1Entry * Wise2_hard_link_PfamHmmer1Entry( Wise2_PfamHmmer1Entry * obj);

/* Function:  Wise2_PfamHmmer1Entry_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_PfamHmmer1Entry *]
 *
 */
Wise2_PfamHmmer1Entry * Wise2_PfamHmmer1Entry_alloc();

/* Function:  Wise2_replace_entryname_PfamHmmer1Entry(obj,entryname)
 *
 * Descrip:    Replace member variable entryname
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1Entry *]
 * Arg:        entryname    New value of the variable [char *]
 *
 * Returns member variable entryname [boolean]
 *
 */
boolean Wise2_replace_entryname_PfamHmmer1Entry( Wise2_PfamHmmer1Entry * obj,char * entryname);

/* Function:  Wise2_access_entryname_PfamHmmer1Entry(obj)
 *
 * Descrip:    Access member variable entryname
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1Entry *]
 *
 * Returns member variable entryname [char *]
 *
 */
char * Wise2_access_entryname_PfamHmmer1Entry( Wise2_PfamHmmer1Entry * obj);

/* Function:  Wise2_replace_is_random_PfamHmmer1Entry(obj,is_random)
 *
 * Descrip:    Replace member variable is_random
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1Entry *]
 * Arg:        is_random    New value of the variable [boolean]
 *
 * Returns member variable is_random [boolean]
 *
 */
boolean Wise2_replace_is_random_PfamHmmer1Entry( Wise2_PfamHmmer1Entry * obj,boolean is_random);

/* Function:  Wise2_access_is_random_PfamHmmer1Entry(obj)
 *
 * Descrip:    Access member variable is_random
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1Entry *]
 *
 * Returns member variable is_random [boolean]
 *
 */
boolean Wise2_access_is_random_PfamHmmer1Entry( Wise2_PfamHmmer1Entry * obj);

/* Function:  Wise2_replace_is_hmmls_PfamHmmer1Entry(obj,is_hmmls)
 *
 * Descrip:    Replace member variable is_hmmls
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1Entry *]
 * Arg:        is_hmmls     New value of the variable [boolean]
 *
 * Returns member variable is_hmmls [boolean]
 *
 */
boolean Wise2_replace_is_hmmls_PfamHmmer1Entry( Wise2_PfamHmmer1Entry * obj,boolean is_hmmls);

/* Function:  Wise2_access_is_hmmls_PfamHmmer1Entry(obj)
 *
 * Descrip:    Access member variable is_hmmls
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1Entry *]
 *
 * Returns member variable is_hmmls [boolean]
 *
 */
boolean Wise2_access_is_hmmls_PfamHmmer1Entry( Wise2_PfamHmmer1Entry * obj);

/* Function:  Wise2_replace_bits_cutoff_PfamHmmer1Entry(obj,bits_cutoff)
 *
 * Descrip:    Replace member variable bits_cutoff
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1Entry *]
 * Arg:        bits_cutoff  New value of the variable [double]
 *
 * Returns member variable bits_cutoff [boolean]
 *
 */
boolean Wise2_replace_bits_cutoff_PfamHmmer1Entry( Wise2_PfamHmmer1Entry * obj,double bits_cutoff);

/* Function:  Wise2_access_bits_cutoff_PfamHmmer1Entry(obj)
 *
 * Descrip:    Access member variable bits_cutoff
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PfamHmmer1Entry *]
 *
 * Returns member variable bits_cutoff [double]
 *
 */
double Wise2_access_bits_cutoff_PfamHmmer1Entry( Wise2_PfamHmmer1Entry * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_PfamHmmer1Entry(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_PfamHmmer1Entry *]
 *
 * Returns Undocumented return value [Wise2_PfamHmmer1Entry *]
 *
 */
Wise2_PfamHmmer1Entry * Wise2_free_PfamHmmer1Entry( Wise2_PfamHmmer1Entry * obj);


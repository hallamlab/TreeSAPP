

/* Functions that create, manipulate or act on RandomProteinDB
 *
 * Wise2_hard_link_RandomProteinDB
 * Wise2_RandomProteinDB_alloc
 * Wise2_replace_use_flat_length_RandomProteinDB
 * Wise2_access_use_flat_length_RandomProteinDB
 * Wise2_replace_length_RandomProteinDB
 * Wise2_access_length_RandomProteinDB
 * Wise2_replace_length_dist_RandomProteinDB
 * Wise2_access_length_dist_RandomProteinDB
 * Wise2_replace_emission_RandomProteinDB
 * Wise2_access_emission_RandomProteinDB
 * Wise2_replace_num_RandomProteinDB
 * Wise2_access_num_RandomProteinDB
 * Wise2_free_RandomProteinDB [destructor]
 *
 */



/* Functions that create, manipulate or act on RandomDNADB
 *
 * Wise2_hard_link_RandomDNADB
 * Wise2_RandomDNADB_alloc
 * Wise2_replace_use_flat_length_RandomDNADB
 * Wise2_access_use_flat_length_RandomDNADB
 * Wise2_replace_length_RandomDNADB
 * Wise2_access_length_RandomDNADB
 * Wise2_replace_length_dist_RandomDNADB
 * Wise2_access_length_dist_RandomDNADB
 * Wise2_replace_emission_RandomDNADB
 * Wise2_access_emission_RandomDNADB
 * Wise2_replace_num_RandomDNADB
 * Wise2_access_num_RandomDNADB
 * Wise2_free_RandomDNADB [destructor]
 *
 */

/* API for object RandomProteinDB */
/* Function:  Wise2_hard_link_RandomProteinDB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_RandomProteinDB *]
 *
 * Returns Undocumented return value [Wise2_RandomProteinDB *]
 *
 */
Wise2_RandomProteinDB * Wise2_hard_link_RandomProteinDB( Wise2_RandomProteinDB * obj);

/* Function:  Wise2_RandomProteinDB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_RandomProteinDB *]
 *
 */
Wise2_RandomProteinDB * Wise2_RandomProteinDB_alloc();

/* Function:  Wise2_replace_use_flat_length_RandomProteinDB(obj,use_flat_length)
 *
 * Descrip:    Replace member variable use_flat_length
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomProteinDB *]
 * Arg:        use_flat_length New value of the variable [boolean]
 *
 * Returns member variable use_flat_length [boolean]
 *
 */
boolean Wise2_replace_use_flat_length_RandomProteinDB( Wise2_RandomProteinDB * obj,boolean use_flat_length);

/* Function:  Wise2_access_use_flat_length_RandomProteinDB(obj)
 *
 * Descrip:    Access member variable use_flat_length
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomProteinDB *]
 *
 * Returns member variable use_flat_length [boolean]
 *
 */
boolean Wise2_access_use_flat_length_RandomProteinDB( Wise2_RandomProteinDB * obj);

/* Function:  Wise2_replace_length_RandomProteinDB(obj,length)
 *
 * Descrip:    Replace member variable length
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomProteinDB *]
 * Arg:        length       New value of the variable [int]
 *
 * Returns member variable length [boolean]
 *
 */
boolean Wise2_replace_length_RandomProteinDB( Wise2_RandomProteinDB * obj,int length);

/* Function:  Wise2_access_length_RandomProteinDB(obj)
 *
 * Descrip:    Access member variable length
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomProteinDB *]
 *
 * Returns member variable length [int]
 *
 */
int Wise2_access_length_RandomProteinDB( Wise2_RandomProteinDB * obj);

/* Function:  Wise2_replace_length_dist_RandomProteinDB(obj,length_dist)
 *
 * Descrip:    Replace member variable length_dist
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomProteinDB *]
 * Arg:        length_dist  New value of the variable [Wise2_Histogram *]
 *
 * Returns member variable length_dist [boolean]
 *
 */
boolean Wise2_replace_length_dist_RandomProteinDB( Wise2_RandomProteinDB * obj,Wise2_Histogram * length_dist);

/* Function:  Wise2_access_length_dist_RandomProteinDB(obj)
 *
 * Descrip:    Access member variable length_dist
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomProteinDB *]
 *
 * Returns member variable length_dist [Wise2_Histogram *]
 *
 */
Wise2_Histogram * Wise2_access_length_dist_RandomProteinDB( Wise2_RandomProteinDB * obj);

/* Function:  Wise2_replace_emission_RandomProteinDB(obj,emission)
 *
 * Descrip:    Replace member variable emission
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomProteinDB *]
 * Arg:        emission     New value of the variable [Wise2_RandomModel *]
 *
 * Returns member variable emission [boolean]
 *
 */
boolean Wise2_replace_emission_RandomProteinDB( Wise2_RandomProteinDB * obj,Wise2_RandomModel * emission);

/* Function:  Wise2_access_emission_RandomProteinDB(obj)
 *
 * Descrip:    Access member variable emission
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomProteinDB *]
 *
 * Returns member variable emission [Wise2_RandomModel *]
 *
 */
Wise2_RandomModel * Wise2_access_emission_RandomProteinDB( Wise2_RandomProteinDB * obj);

/* Function:  Wise2_replace_num_RandomProteinDB(obj,num)
 *
 * Descrip:    Replace member variable num
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomProteinDB *]
 * Arg:        num          New value of the variable [int]
 *
 * Returns member variable num [boolean]
 *
 */
boolean Wise2_replace_num_RandomProteinDB( Wise2_RandomProteinDB * obj,int num);

/* Function:  Wise2_access_num_RandomProteinDB(obj)
 *
 * Descrip:    Access member variable num
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomProteinDB *]
 *
 * Returns member variable num [int]
 *
 */
int Wise2_access_num_RandomProteinDB( Wise2_RandomProteinDB * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_RandomProteinDB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_RandomProteinDB *]
 *
 * Returns Undocumented return value [Wise2_RandomProteinDB *]
 *
 */
Wise2_RandomProteinDB * Wise2_free_RandomProteinDB( Wise2_RandomProteinDB * obj);

/* API for object RandomDNADB */
/* Function:  Wise2_hard_link_RandomDNADB(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_RandomDNADB *]
 *
 * Returns Undocumented return value [Wise2_RandomDNADB *]
 *
 */
Wise2_RandomDNADB * Wise2_hard_link_RandomDNADB( Wise2_RandomDNADB * obj);

/* Function:  Wise2_RandomDNADB_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_RandomDNADB *]
 *
 */
Wise2_RandomDNADB * Wise2_RandomDNADB_alloc();

/* Function:  Wise2_replace_use_flat_length_RandomDNADB(obj,use_flat_length)
 *
 * Descrip:    Replace member variable use_flat_length
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomDNADB *]
 * Arg:        use_flat_length New value of the variable [boolean]
 *
 * Returns member variable use_flat_length [boolean]
 *
 */
boolean Wise2_replace_use_flat_length_RandomDNADB( Wise2_RandomDNADB * obj,boolean use_flat_length);

/* Function:  Wise2_access_use_flat_length_RandomDNADB(obj)
 *
 * Descrip:    Access member variable use_flat_length
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomDNADB *]
 *
 * Returns member variable use_flat_length [boolean]
 *
 */
boolean Wise2_access_use_flat_length_RandomDNADB( Wise2_RandomDNADB * obj);

/* Function:  Wise2_replace_length_RandomDNADB(obj,length)
 *
 * Descrip:    Replace member variable length
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomDNADB *]
 * Arg:        length       New value of the variable [int]
 *
 * Returns member variable length [boolean]
 *
 */
boolean Wise2_replace_length_RandomDNADB( Wise2_RandomDNADB * obj,int length);

/* Function:  Wise2_access_length_RandomDNADB(obj)
 *
 * Descrip:    Access member variable length
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomDNADB *]
 *
 * Returns member variable length [int]
 *
 */
int Wise2_access_length_RandomDNADB( Wise2_RandomDNADB * obj);

/* Function:  Wise2_replace_length_dist_RandomDNADB(obj,length_dist)
 *
 * Descrip:    Replace member variable length_dist
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomDNADB *]
 * Arg:        length_dist  New value of the variable [Wise2_Histogram *]
 *
 * Returns member variable length_dist [boolean]
 *
 */
boolean Wise2_replace_length_dist_RandomDNADB( Wise2_RandomDNADB * obj,Wise2_Histogram * length_dist);

/* Function:  Wise2_access_length_dist_RandomDNADB(obj)
 *
 * Descrip:    Access member variable length_dist
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomDNADB *]
 *
 * Returns member variable length_dist [Wise2_Histogram *]
 *
 */
Wise2_Histogram * Wise2_access_length_dist_RandomDNADB( Wise2_RandomDNADB * obj);

/* Function:  Wise2_replace_emission_RandomDNADB(obj,emission)
 *
 * Descrip:    Replace member variable emission
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomDNADB *]
 * Arg:        emission     New value of the variable [Wise2_RandomModelDNA *]
 *
 * Returns member variable emission [boolean]
 *
 */
boolean Wise2_replace_emission_RandomDNADB( Wise2_RandomDNADB * obj,Wise2_RandomModelDNA * emission);

/* Function:  Wise2_access_emission_RandomDNADB(obj)
 *
 * Descrip:    Access member variable emission
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomDNADB *]
 *
 * Returns member variable emission [Wise2_RandomModelDNA *]
 *
 */
Wise2_RandomModelDNA * Wise2_access_emission_RandomDNADB( Wise2_RandomDNADB * obj);

/* Function:  Wise2_replace_num_RandomDNADB(obj,num)
 *
 * Descrip:    Replace member variable num
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomDNADB *]
 * Arg:        num          New value of the variable [int]
 *
 * Returns member variable num [boolean]
 *
 */
boolean Wise2_replace_num_RandomDNADB( Wise2_RandomDNADB * obj,int num);

/* Function:  Wise2_access_num_RandomDNADB(obj)
 *
 * Descrip:    Access member variable num
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomDNADB *]
 *
 * Returns member variable num [int]
 *
 */
int Wise2_access_num_RandomDNADB( Wise2_RandomDNADB * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_RandomDNADB(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_RandomDNADB *]
 *
 * Returns Undocumented return value [Wise2_RandomDNADB *]
 *
 */
Wise2_RandomDNADB * Wise2_free_RandomDNADB( Wise2_RandomDNADB * obj);


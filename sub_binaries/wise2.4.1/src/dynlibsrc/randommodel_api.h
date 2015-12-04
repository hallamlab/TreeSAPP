

/* Functions that create, manipulate or act on RandomModelDNA
 *
 * Wise2_hard_link_RandomModelDNA
 * Wise2_RandomModelDNA_alloc
 * Wise2_replace_name_RandomModelDNA
 * Wise2_access_name_RandomModelDNA
 * Wise2_free_RandomModelDNA [destructor]
 *
 */



/* Functions that create, manipulate or act on RandomModel
 *
 * Wise2_hard_link_RandomModel
 * Wise2_RandomModel_alloc
 * Wise2_replace_name_RandomModel
 * Wise2_access_name_RandomModel
 * Wise2_free_RandomModel [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_RandomModelDNA_std
 * Wise2_default_RandomModel
 *

/* API for object RandomModelDNA */
/* Function:  Wise2_hard_link_RandomModelDNA(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_RandomModelDNA *]
 *
 * Returns Undocumented return value [Wise2_RandomModelDNA *]
 *
 */
Wise2_RandomModelDNA * Wise2_hard_link_RandomModelDNA( Wise2_RandomModelDNA * obj);

/* Function:  Wise2_RandomModelDNA_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_RandomModelDNA *]
 *
 */
Wise2_RandomModelDNA * Wise2_RandomModelDNA_alloc();

/* Function:  Wise2_replace_name_RandomModelDNA(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomModelDNA *]
 * Arg:        name         New value of the variable [char *]
 *
 * Returns member variable name [boolean]
 *
 */
boolean Wise2_replace_name_RandomModelDNA( Wise2_RandomModelDNA * obj,char * name);

/* Function:  Wise2_access_name_RandomModelDNA(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomModelDNA *]
 *
 * Returns member variable name [char *]
 *
 */
char * Wise2_access_name_RandomModelDNA( Wise2_RandomModelDNA * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_RandomModelDNA(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_RandomModelDNA *]
 *
 * Returns Undocumented return value [Wise2_RandomModelDNA *]
 *
 */
Wise2_RandomModelDNA * Wise2_free_RandomModelDNA( Wise2_RandomModelDNA * obj);

/* API for object RandomModel */
/* Function:  Wise2_hard_link_RandomModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_RandomModel *]
 *
 * Returns Undocumented return value [Wise2_RandomModel *]
 *
 */
Wise2_RandomModel * Wise2_hard_link_RandomModel( Wise2_RandomModel * obj);

/* Function:  Wise2_RandomModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_RandomModel *]
 *
 */
Wise2_RandomModel * Wise2_RandomModel_alloc();

/* Function:  Wise2_replace_name_RandomModel(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomModel *]
 * Arg:        name         New value of the variable [char *]
 *
 * Returns member variable name [boolean]
 *
 */
boolean Wise2_replace_name_RandomModel( Wise2_RandomModel * obj,char * name);

/* Function:  Wise2_access_name_RandomModel(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_RandomModel *]
 *
 * Returns member variable name [char *]
 *
 */
char * Wise2_access_name_RandomModel( Wise2_RandomModel * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_RandomModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_RandomModel *]
 *
 * Returns Undocumented return value [Wise2_RandomModel *]
 *
 */
Wise2_RandomModel * Wise2_free_RandomModel( Wise2_RandomModel * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_RandomModelDNA_std(void)
 *
 * Descrip:    Returns a structure with 0.25 set in each place
 *
 *
 *
 * Returns Undocumented return value [Wise2_RandomModelDNA *]
 *
 */
Wise2_RandomModelDNA * Wise2_RandomModelDNA_std();

/* Function:  Wise2_default_RandomModel(void)
 *
 * Descrip:    Gives a default random model numbers from
 *             swissprot34- via the HMMEr1 package
 *
 *
 *
 * Returns Undocumented return value [Wise2_RandomModel *]
 *
 */
Wise2_RandomModel * Wise2_default_RandomModel();




/* Functions that create, manipulate or act on cDNAParser
 *
 * Wise2_hard_link_cDNAParser
 * Wise2_cDNAParser_alloc
 * Wise2_free_cDNAParser [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_flat_cDNAParser
 *

/* API for object cDNAParser */
/* Function:  Wise2_hard_link_cDNAParser(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_cDNAParser *]
 *
 * Returns Undocumented return value [Wise2_cDNAParser *]
 *
 */
Wise2_cDNAParser * Wise2_hard_link_cDNAParser( Wise2_cDNAParser * obj);

/* Function:  Wise2_cDNAParser_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_cDNAParser *]
 *
 */
Wise2_cDNAParser * Wise2_cDNAParser_alloc();

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_cDNAParser(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_cDNAParser *]
 *
 * Returns Undocumented return value [Wise2_cDNAParser *]
 *
 */
Wise2_cDNAParser * Wise2_free_cDNAParser( Wise2_cDNAParser * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_flat_cDNAParser(p)
 *
 * Descrip:    Makes a flat (ie, indels of 1 or 2 == p)
 *             cDNA parser. This means that insertions
 *             and deletions of both 1 or 2 bases are
 *             all parameterised at the same probability
 *
 *
 *
 * Arg:        p            probability of an indel [Probability]
 *
 * Returns Undocumented return value [Wise2_cDNAParser *]
 *
 */
Wise2_cDNAParser * Wise2_flat_cDNAParser( Probability p);


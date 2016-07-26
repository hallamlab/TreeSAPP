

/* Functions that create, manipulate or act on GeneParameter21
 *
 * Wise2_hard_link_GeneParameter21
 * Wise2_GeneParameter21_alloc_std
 * Wise2_free_GeneParameter21 [destructor]
 *
 */

/* API for object GeneParameter21 */
/* Function:  Wise2_hard_link_GeneParameter21(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_GeneParameter21 *]
 *
 * Returns Undocumented return value [Wise2_GeneParameter21 *]
 *
 */
Wise2_GeneParameter21 * Wise2_hard_link_GeneParameter21( Wise2_GeneParameter21 * obj);

/* Function:  Wise2_GeneParameter21_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneParameter21_alloc_len(GeneParameter21LISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_GeneParameter21 *]
 *
 */
Wise2_GeneParameter21 * Wise2_GeneParameter21_alloc_std();

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_GeneParameter21(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_GeneParameter21 *]
 *
 * Returns Undocumented return value [Wise2_GeneParameter21 *]
 *
 */
Wise2_GeneParameter21 * Wise2_free_GeneParameter21( Wise2_GeneParameter21 * obj);


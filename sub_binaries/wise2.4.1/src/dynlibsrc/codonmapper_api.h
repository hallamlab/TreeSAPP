

/* Functions that create, manipulate or act on CodonMapper
 *
 * Wise2_sprinkle_errors_over_CodonMapper
 * Wise2_hard_link_CodonMapper
 * Wise2_CodonMapper_alloc
 * Wise2_replace_ct_CodonMapper
 * Wise2_access_ct_CodonMapper
 * Wise2_free_CodonMapper [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_flat_CodonMapper
 *

/* API for object CodonMapper */
/* Function:  Wise2_sprinkle_errors_over_CodonMapper(cm,error)
 *
 * Descrip:    Takes a codon mapper and assummes that the majority of errors
 *             are due to a single base change in the codon at probability error.
 *             Therefore, for each codon it adds error * prob(codon) * 0.25 to each 
 *             other codon one base away, taking away therefore the result.
 *
 *
 *
 * Arg:        cm           CodonMapper to be sprinkled [Wise2_CodonMapper *]
 * Arg:        error        substitution error rate [double]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_sprinkle_errors_over_CodonMapper( Wise2_CodonMapper * cm,double error);

/* Function:  Wise2_hard_link_CodonMapper(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_CodonMapper *]
 *
 * Returns Undocumented return value [Wise2_CodonMapper *]
 *
 */
Wise2_CodonMapper * Wise2_hard_link_CodonMapper( Wise2_CodonMapper * obj);

/* Function:  Wise2_CodonMapper_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_CodonMapper *]
 *
 */
Wise2_CodonMapper * Wise2_CodonMapper_alloc();

/* Function:  Wise2_replace_ct_CodonMapper(obj,ct)
 *
 * Descrip:    Replace member variable ct
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_CodonMapper *]
 * Arg:        ct           New value of the variable [Wise2_CodonTable *]
 *
 * Returns member variable ct [boolean]
 *
 */
boolean Wise2_replace_ct_CodonMapper( Wise2_CodonMapper * obj,Wise2_CodonTable * ct);

/* Function:  Wise2_access_ct_CodonMapper(obj)
 *
 * Descrip:    Access member variable ct
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_CodonMapper *]
 *
 * Returns member variable ct [Wise2_CodonTable *]
 *
 */
Wise2_CodonTable * Wise2_access_ct_CodonMapper( Wise2_CodonMapper * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_CodonMapper(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_CodonMapper *]
 *
 * Returns Undocumented return value [Wise2_CodonMapper *]
 *
 */
Wise2_CodonMapper * Wise2_free_CodonMapper( Wise2_CodonMapper * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_flat_CodonMapper(ct)
 *
 * Descrip:    Makes a CodonMapper with no codon bias
 *             or error possiblities from codon table
 *
 *
 *
 * Arg:        ct           Codon Table giving codon->aa info [Wise2_CodonTable *]
 *
 * Returns Undocumented return value [Wise2_CodonMapper *]
 *
 */
Wise2_CodonMapper * Wise2_flat_CodonMapper( Wise2_CodonTable * ct);


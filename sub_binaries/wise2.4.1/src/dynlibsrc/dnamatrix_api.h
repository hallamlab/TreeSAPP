

/* Functions that create, manipulate or act on DnaMatrix
 *
 * Wise2_hard_link_DnaMatrix
 * Wise2_DnaMatrix_alloc
 * Wise2_free_DnaMatrix [destructor]
 *
 */



/* Functions that create, manipulate or act on DnaProbMatrix
 *
 * Wise2_flat_null_DnaProbMatrix
 * Wise2_hard_link_DnaProbMatrix
 * Wise2_DnaProbMatrix_alloc
 * Wise2_free_DnaProbMatrix [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_identity_DnaMatrix
 * Wise2_DnaProbMatrix_from_match
 * Wise2_DnaMatrix_from_DnaProbMatrix
 *

/* API for object DnaMatrix */
/* Function:  Wise2_hard_link_DnaMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_DnaMatrix *]
 *
 * Returns Undocumented return value [Wise2_DnaMatrix *]
 *
 */
Wise2_DnaMatrix * Wise2_hard_link_DnaMatrix( Wise2_DnaMatrix * obj);

/* Function:  Wise2_DnaMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_DnaMatrix *]
 *
 */
Wise2_DnaMatrix * Wise2_DnaMatrix_alloc();

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_DnaMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_DnaMatrix *]
 *
 * Returns Undocumented return value [Wise2_DnaMatrix *]
 *
 */
Wise2_DnaMatrix * Wise2_free_DnaMatrix( Wise2_DnaMatrix * obj);

/* API for object DnaProbMatrix */
/* Function:  Wise2_flat_null_DnaProbMatrix(dpm)
 *
 * Descrip:    makes a odds of dpm via a 0.25 factor 
 *             into each base.
 *
 *
 * Arg:        dpm          Undocumented argument [Wise2_DnaProbMatrix *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_flat_null_DnaProbMatrix( Wise2_DnaProbMatrix * dpm);

/* Function:  Wise2_hard_link_DnaProbMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_DnaProbMatrix *]
 *
 * Returns Undocumented return value [Wise2_DnaProbMatrix *]
 *
 */
Wise2_DnaProbMatrix * Wise2_hard_link_DnaProbMatrix( Wise2_DnaProbMatrix * obj);

/* Function:  Wise2_DnaProbMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_DnaProbMatrix *]
 *
 */
Wise2_DnaProbMatrix * Wise2_DnaProbMatrix_alloc();

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_DnaProbMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_DnaProbMatrix *]
 *
 * Returns Undocumented return value [Wise2_DnaProbMatrix *]
 *
 */
Wise2_DnaProbMatrix * Wise2_free_DnaProbMatrix( Wise2_DnaProbMatrix * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_identity_DnaMatrix(id_score,mismatch)
 *
 * Descrip:    makes an idenity matrix wth id_score on the leading
 *             diagonal and mismatch elsewhere.
 *
 *
 *
 * Arg:        id_score     score of idenities [Score]
 * Arg:        mismatch     score of mistmatches [Score]
 *
 * Returns Undocumented return value [Wise2_DnaMatrix *]
 *
 */
Wise2_DnaMatrix * Wise2_identity_DnaMatrix( Score id_score,Score mismatch);

/* Function:  Wise2_DnaProbMatrix_from_match(match,nmask_type)
 *
 * Descrip:    Makes a probability matrix from simple match/mismatch 
 *             probabilities.
 *
 *
 *
 * Arg:        match        Undocumented argument [Probability]
 * Arg:        nmask_type   Undocumented argument [int]
 *
 * Returns Undocumented return value [Wise2_DnaProbMatrix *]
 *
 */
Wise2_DnaProbMatrix * Wise2_DnaProbMatrix_from_match( Probability match,int nmask_type);

/* Function:  Wise2_DnaMatrix_from_DnaProbMatrix(dpm)
 *
 * Descrip:    Maps probabilities to scores
 *
 *
 * Arg:        dpm          Undocumented argument [Wise2_DnaProbMatrix *]
 *
 * Returns Undocumented return value [Wise2_DnaMatrix *]
 *
 */
Wise2_DnaMatrix * Wise2_DnaMatrix_from_DnaProbMatrix( Wise2_DnaProbMatrix * dpm);


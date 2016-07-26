

/* Functions that create, manipulate or act on DnaStartEnd
 *
 * Wise2_hard_link_DnaStartEnd
 * Wise2_DnaStartEnd_alloc
 * Wise2_free_DnaStartEnd [destructor]
 *
 */



/* Helper functions in the module
 *
 * Wise2_make_align_dnaalign
 * Wise2_DnaStartEnd_from_policy
 *

/* API for object DnaStartEnd */
/* Function:  Wise2_hard_link_DnaStartEnd(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_DnaStartEnd *]
 *
 * Returns Undocumented return value [Wise2_DnaStartEnd *]
 *
 */
Wise2_DnaStartEnd * Wise2_hard_link_DnaStartEnd( Wise2_DnaStartEnd * obj);

/* Function:  Wise2_DnaStartEnd_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_DnaStartEnd *]
 *
 */
Wise2_DnaStartEnd * Wise2_DnaStartEnd_alloc();

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_DnaStartEnd(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_DnaStartEnd *]
 *
 * Returns Undocumented return value [Wise2_DnaStartEnd *]
 *
 */
Wise2_DnaStartEnd * Wise2_free_DnaStartEnd( Wise2_DnaStartEnd * obj);



/* These functions are not associated with an object */
/* Function:  Wise2_make_align_dnaalign(one,two,mat,se,qgap,qext,tgap,text,dpri)
 *
 * Descrip:    Makes an alignment out of two DNA sequences
 *
 *
 * Arg:        one          first sequence to align [Wise2_Sequence *]
 * Arg:        two          second sequence to align [Wise2_Sequence *]
 * Arg:        mat          DnaMatrix for the matching [Wise2_DnaMatrix *]
 * Arg:        se           DnaStartEnd policy [Wise2_DnaStartEnd *]
 * Arg:        qgap         gap open penalty in query (one) coordinate [int]
 * Arg:        qext         gap extension penalty in query (one) coordinate [int]
 * Arg:        tgap         gap open penalty in target (two) coordinate [int]
 * Arg:        text         gap extension penalty in target (two) coordinate [int]
 * Arg:        dpri         DPRunImpl structure [Wise2_DPRunImpl *]
 *
 * Returns an alb structure of the alignment [Wise2_AlnBlock *]
 *
 */
Wise2_AlnBlock * Wise2_make_align_dnaalign( Wise2_Sequence * one,Wise2_Sequence * two,Wise2_DnaMatrix * mat,Wise2_DnaStartEnd * se,int qgap,int qext,int tgap,int text,Wise2_DPRunImpl * dpri);

/* Function:  Wise2_DnaStartEnd_from_policy(policy)
 *
 * Descrip:    Makes a DnaStartEnd from a particular string.
 *             Possible strings are:
 *
 *             local - fully local
 *             global - fully global
 *             edge - aligns only to edges
 *
 *
 * Arg:        policy       Undocumented argument [char *]
 *
 * Returns Undocumented return value [Wise2_DnaStartEnd *]
 *
 */
Wise2_DnaStartEnd * Wise2_DnaStartEnd_from_policy( char * policy);


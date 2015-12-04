

/* Functions that create, manipulate or act on PackAln
 *
 * Wise2_show_simple_PackAln
 * Wise2_show_bits_and_cumlative_PackAln
 * Wise2_hard_link_PackAln
 * Wise2_PackAln_alloc_std
 * Wise2_access_pau_PackAln
 * Wise2_length_pau_PackAln
 * Wise2_flush_PackAln
 * Wise2_add_PackAln
 * Wise2_replace_score_PackAln
 * Wise2_access_score_PackAln
 * Wise2_free_PackAln [destructor]
 *
 */



/* Functions that create, manipulate or act on PackAlnUnit
 *
 * Wise2_hard_link_PackAlnUnit
 * Wise2_PackAlnUnit_alloc
 * Wise2_replace_i_PackAlnUnit
 * Wise2_access_i_PackAlnUnit
 * Wise2_replace_j_PackAlnUnit
 * Wise2_access_j_PackAlnUnit
 * Wise2_replace_state_PackAlnUnit
 * Wise2_access_state_PackAlnUnit
 * Wise2_replace_score_PackAlnUnit
 * Wise2_access_score_PackAlnUnit
 * Wise2_free_PackAlnUnit [destructor]
 *
 */

/* API for object PackAln */
/* Function:  Wise2_show_simple_PackAln(pal,ofp)
 *
 * Descrip:    shows packaln with a pretty verbose debugging 
 *             format
 *
 *
 * Arg:        pal          Undocumented argument [Wise2_PackAln *]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_show_simple_PackAln( Wise2_PackAln * pal,FILE * ofp);

/* Function:  Wise2_show_bits_and_cumlative_PackAln(pal,ofp)
 *
 * Descrip:    Shows packaln as: 
 *
 *             i,j,state,score,bits,cumlative-score,cumlative-bits
 *
 *             cumlative score and cumlative bits are useful sometimes
 *
 *
 * Arg:        pal          Undocumented argument [Wise2_PackAln *]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_show_bits_and_cumlative_PackAln( Wise2_PackAln * pal,FILE * ofp);

/* Function:  Wise2_hard_link_PackAln(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_PackAln *]
 *
 * Returns Undocumented return value [Wise2_PackAln *]
 *
 */
Wise2_PackAln * Wise2_hard_link_PackAln( Wise2_PackAln * obj);

/* Function:  Wise2_PackAln_alloc_std(void)
 *
 * Descrip:    Equivalent to PackAln_alloc_len(PackAlnLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_PackAln *]
 *
 */
Wise2_PackAln * Wise2_PackAln_alloc_std();

/* Function:  Wise2_access_pau_PackAln(obj,i)
 *
 * Descrip:    Access members stored in the pau list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_PackAln *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_PackAlnUnit *]
 *
 */
Wise2_PackAlnUnit * Wise2_access_pau_PackAln( Wise2_PackAln * obj,int i);

/* Function:  Wise2_length_pau_PackAln(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_PackAln *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_pau_PackAln( Wise2_PackAln * obj);

/* Function:  Wise2_flush_PackAln(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_PackAln *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_PackAln( Wise2_PackAln * obj);

/* Function:  Wise2_add_PackAln(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_PackAln *]
 * Arg:        add          Object to add to the list [Wise2_PackAlnUnit *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_PackAln( Wise2_PackAln * obj,Wise2_PackAlnUnit * add);

/* Function:  Wise2_replace_score_PackAln(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PackAln *]
 * Arg:        score        New value of the variable [int]
 *
 * Returns member variable score [boolean]
 *
 */
boolean Wise2_replace_score_PackAln( Wise2_PackAln * obj,int score);

/* Function:  Wise2_access_score_PackAln(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PackAln *]
 *
 * Returns member variable score [int]
 *
 */
int Wise2_access_score_PackAln( Wise2_PackAln * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_PackAln(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_PackAln *]
 *
 * Returns Undocumented return value [Wise2_PackAln *]
 *
 */
Wise2_PackAln * Wise2_free_PackAln( Wise2_PackAln * obj);

/* API for object PackAlnUnit */
/* Function:  Wise2_hard_link_PackAlnUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_PackAlnUnit *]
 *
 * Returns Undocumented return value [Wise2_PackAlnUnit *]
 *
 */
Wise2_PackAlnUnit * Wise2_hard_link_PackAlnUnit( Wise2_PackAlnUnit * obj);

/* Function:  Wise2_PackAlnUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_PackAlnUnit *]
 *
 */
Wise2_PackAlnUnit * Wise2_PackAlnUnit_alloc();

/* Function:  Wise2_replace_i_PackAlnUnit(obj,i)
 *
 * Descrip:    Replace member variable i
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PackAlnUnit *]
 * Arg:        i            New value of the variable [int]
 *
 * Returns member variable i [boolean]
 *
 */
boolean Wise2_replace_i_PackAlnUnit( Wise2_PackAlnUnit * obj,int i);

/* Function:  Wise2_access_i_PackAlnUnit(obj)
 *
 * Descrip:    Access member variable i
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PackAlnUnit *]
 *
 * Returns member variable i [int]
 *
 */
int Wise2_access_i_PackAlnUnit( Wise2_PackAlnUnit * obj);

/* Function:  Wise2_replace_j_PackAlnUnit(obj,j)
 *
 * Descrip:    Replace member variable j
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PackAlnUnit *]
 * Arg:        j            New value of the variable [int]
 *
 * Returns member variable j [boolean]
 *
 */
boolean Wise2_replace_j_PackAlnUnit( Wise2_PackAlnUnit * obj,int j);

/* Function:  Wise2_access_j_PackAlnUnit(obj)
 *
 * Descrip:    Access member variable j
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PackAlnUnit *]
 *
 * Returns member variable j [int]
 *
 */
int Wise2_access_j_PackAlnUnit( Wise2_PackAlnUnit * obj);

/* Function:  Wise2_replace_state_PackAlnUnit(obj,state)
 *
 * Descrip:    Replace member variable state
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PackAlnUnit *]
 * Arg:        state        New value of the variable [int]
 *
 * Returns member variable state [boolean]
 *
 */
boolean Wise2_replace_state_PackAlnUnit( Wise2_PackAlnUnit * obj,int state);

/* Function:  Wise2_access_state_PackAlnUnit(obj)
 *
 * Descrip:    Access member variable state
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PackAlnUnit *]
 *
 * Returns member variable state [int]
 *
 */
int Wise2_access_state_PackAlnUnit( Wise2_PackAlnUnit * obj);

/* Function:  Wise2_replace_score_PackAlnUnit(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PackAlnUnit *]
 * Arg:        score        New value of the variable [int]
 *
 * Returns member variable score [boolean]
 *
 */
boolean Wise2_replace_score_PackAlnUnit( Wise2_PackAlnUnit * obj,int score);

/* Function:  Wise2_access_score_PackAlnUnit(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_PackAlnUnit *]
 *
 * Returns member variable score [int]
 *
 */
int Wise2_access_score_PackAlnUnit( Wise2_PackAlnUnit * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_PackAlnUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_PackAlnUnit *]
 *
 * Returns Undocumented return value [Wise2_PackAlnUnit *]
 *
 */
Wise2_PackAlnUnit * Wise2_free_PackAlnUnit( Wise2_PackAlnUnit * obj);




/* Functions that create, manipulate or act on AlnRange
 *
 * Wise2_hard_link_AlnRange
 * Wise2_AlnRange_alloc
 * Wise2_replace_starti_AlnRange
 * Wise2_access_starti_AlnRange
 * Wise2_replace_startj_AlnRange
 * Wise2_access_startj_AlnRange
 * Wise2_replace_startstate_AlnRange
 * Wise2_access_startstate_AlnRange
 * Wise2_replace_stopi_AlnRange
 * Wise2_access_stopi_AlnRange
 * Wise2_replace_stopj_AlnRange
 * Wise2_access_stopj_AlnRange
 * Wise2_replace_stopstate_AlnRange
 * Wise2_access_stopstate_AlnRange
 * Wise2_replace_startscore_AlnRange
 * Wise2_access_startscore_AlnRange
 * Wise2_replace_stopscore_AlnRange
 * Wise2_access_stopscore_AlnRange
 * Wise2_free_AlnRange [destructor]
 *
 */



/* Functions that create, manipulate or act on AlnRangeSet
 *
 * Wise2_hard_link_AlnRangeSet
 * Wise2_AlnRangeSet_alloc_std
 * Wise2_replace_score_AlnRangeSet
 * Wise2_access_score_AlnRangeSet
 * Wise2_access_alr_AlnRangeSet
 * Wise2_length_alr_AlnRangeSet
 * Wise2_flush_AlnRangeSet
 * Wise2_add_AlnRangeSet
 * Wise2_free_AlnRangeSet [destructor]
 *
 */

/* API for object AlnRange */
/* Function:  Wise2_hard_link_AlnRange(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_AlnRange *]
 *
 * Returns Undocumented return value [Wise2_AlnRange *]
 *
 */
Wise2_AlnRange * Wise2_hard_link_AlnRange( Wise2_AlnRange * obj);

/* Function:  Wise2_AlnRange_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_AlnRange *]
 *
 */
Wise2_AlnRange * Wise2_AlnRange_alloc();

/* Function:  Wise2_replace_starti_AlnRange(obj,starti)
 *
 * Descrip:    Replace member variable starti
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 * Arg:        starti       New value of the variable [int]
 *
 * Returns member variable starti [boolean]
 *
 */
boolean Wise2_replace_starti_AlnRange( Wise2_AlnRange * obj,int starti);

/* Function:  Wise2_access_starti_AlnRange(obj)
 *
 * Descrip:    Access member variable starti
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 *
 * Returns member variable starti [int]
 *
 */
int Wise2_access_starti_AlnRange( Wise2_AlnRange * obj);

/* Function:  Wise2_replace_startj_AlnRange(obj,startj)
 *
 * Descrip:    Replace member variable startj
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 * Arg:        startj       New value of the variable [int]
 *
 * Returns member variable startj [boolean]
 *
 */
boolean Wise2_replace_startj_AlnRange( Wise2_AlnRange * obj,int startj);

/* Function:  Wise2_access_startj_AlnRange(obj)
 *
 * Descrip:    Access member variable startj
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 *
 * Returns member variable startj [int]
 *
 */
int Wise2_access_startj_AlnRange( Wise2_AlnRange * obj);

/* Function:  Wise2_replace_startstate_AlnRange(obj,startstate)
 *
 * Descrip:    Replace member variable startstate
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 * Arg:        startstate   New value of the variable [int]
 *
 * Returns member variable startstate [boolean]
 *
 */
boolean Wise2_replace_startstate_AlnRange( Wise2_AlnRange * obj,int startstate);

/* Function:  Wise2_access_startstate_AlnRange(obj)
 *
 * Descrip:    Access member variable startstate
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 *
 * Returns member variable startstate [int]
 *
 */
int Wise2_access_startstate_AlnRange( Wise2_AlnRange * obj);

/* Function:  Wise2_replace_stopi_AlnRange(obj,stopi)
 *
 * Descrip:    Replace member variable stopi
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 * Arg:        stopi        New value of the variable [int]
 *
 * Returns member variable stopi [boolean]
 *
 */
boolean Wise2_replace_stopi_AlnRange( Wise2_AlnRange * obj,int stopi);

/* Function:  Wise2_access_stopi_AlnRange(obj)
 *
 * Descrip:    Access member variable stopi
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 *
 * Returns member variable stopi [int]
 *
 */
int Wise2_access_stopi_AlnRange( Wise2_AlnRange * obj);

/* Function:  Wise2_replace_stopj_AlnRange(obj,stopj)
 *
 * Descrip:    Replace member variable stopj
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 * Arg:        stopj        New value of the variable [int]
 *
 * Returns member variable stopj [boolean]
 *
 */
boolean Wise2_replace_stopj_AlnRange( Wise2_AlnRange * obj,int stopj);

/* Function:  Wise2_access_stopj_AlnRange(obj)
 *
 * Descrip:    Access member variable stopj
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 *
 * Returns member variable stopj [int]
 *
 */
int Wise2_access_stopj_AlnRange( Wise2_AlnRange * obj);

/* Function:  Wise2_replace_stopstate_AlnRange(obj,stopstate)
 *
 * Descrip:    Replace member variable stopstate
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 * Arg:        stopstate    New value of the variable [int]
 *
 * Returns member variable stopstate [boolean]
 *
 */
boolean Wise2_replace_stopstate_AlnRange( Wise2_AlnRange * obj,int stopstate);

/* Function:  Wise2_access_stopstate_AlnRange(obj)
 *
 * Descrip:    Access member variable stopstate
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 *
 * Returns member variable stopstate [int]
 *
 */
int Wise2_access_stopstate_AlnRange( Wise2_AlnRange * obj);

/* Function:  Wise2_replace_startscore_AlnRange(obj,startscore)
 *
 * Descrip:    Replace member variable startscore
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 * Arg:        startscore   New value of the variable [int]
 *
 * Returns member variable startscore [boolean]
 *
 */
boolean Wise2_replace_startscore_AlnRange( Wise2_AlnRange * obj,int startscore);

/* Function:  Wise2_access_startscore_AlnRange(obj)
 *
 * Descrip:    Access member variable startscore
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 *
 * Returns member variable startscore [int]
 *
 */
int Wise2_access_startscore_AlnRange( Wise2_AlnRange * obj);

/* Function:  Wise2_replace_stopscore_AlnRange(obj,stopscore)
 *
 * Descrip:    Replace member variable stopscore
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 * Arg:        stopscore    New value of the variable [int]
 *
 * Returns member variable stopscore [boolean]
 *
 */
boolean Wise2_replace_stopscore_AlnRange( Wise2_AlnRange * obj,int stopscore);

/* Function:  Wise2_access_stopscore_AlnRange(obj)
 *
 * Descrip:    Access member variable stopscore
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRange *]
 *
 * Returns member variable stopscore [int]
 *
 */
int Wise2_access_stopscore_AlnRange( Wise2_AlnRange * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_AlnRange(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_AlnRange *]
 *
 * Returns Undocumented return value [Wise2_AlnRange *]
 *
 */
Wise2_AlnRange * Wise2_free_AlnRange( Wise2_AlnRange * obj);

/* API for object AlnRangeSet */
/* Function:  Wise2_hard_link_AlnRangeSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_AlnRangeSet *]
 *
 * Returns Undocumented return value [Wise2_AlnRangeSet *]
 *
 */
Wise2_AlnRangeSet * Wise2_hard_link_AlnRangeSet( Wise2_AlnRangeSet * obj);

/* Function:  Wise2_AlnRangeSet_alloc_std(void)
 *
 * Descrip:    Equivalent to AlnRangeSet_alloc_len(AlnRangeSetLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_AlnRangeSet *]
 *
 */
Wise2_AlnRangeSet * Wise2_AlnRangeSet_alloc_std();

/* Function:  Wise2_replace_score_AlnRangeSet(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRangeSet *]
 * Arg:        score        New value of the variable [int]
 *
 * Returns member variable score [boolean]
 *
 */
boolean Wise2_replace_score_AlnRangeSet( Wise2_AlnRangeSet * obj,int score);

/* Function:  Wise2_access_score_AlnRangeSet(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_AlnRangeSet *]
 *
 * Returns member variable score [int]
 *
 */
int Wise2_access_score_AlnRangeSet( Wise2_AlnRangeSet * obj);

/* Function:  Wise2_access_alr_AlnRangeSet(obj,i)
 *
 * Descrip:    Access members stored in the alr list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_AlnRangeSet *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_AlnRange *]
 *
 */
Wise2_AlnRange * Wise2_access_alr_AlnRangeSet( Wise2_AlnRangeSet * obj,int i);

/* Function:  Wise2_length_alr_AlnRangeSet(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_AlnRangeSet *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_alr_AlnRangeSet( Wise2_AlnRangeSet * obj);

/* Function:  Wise2_flush_AlnRangeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_AlnRangeSet *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_AlnRangeSet( Wise2_AlnRangeSet * obj);

/* Function:  Wise2_add_AlnRangeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_AlnRangeSet *]
 * Arg:        add          Object to add to the list [Wise2_AlnRange *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AlnRangeSet( Wise2_AlnRangeSet * obj,Wise2_AlnRange * add);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_AlnRangeSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_AlnRangeSet *]
 *
 * Returns Undocumented return value [Wise2_AlnRangeSet *]
 *
 */
Wise2_AlnRangeSet * Wise2_free_AlnRangeSet( Wise2_AlnRangeSet * obj);


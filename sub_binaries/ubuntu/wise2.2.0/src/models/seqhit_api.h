

/* Functions that create, manipulate or act on DnaSequenceHitList
 *
 * Wise2_show_DnaSequenceHitList
 * Wise2_read_MSPcrunch_DnaSequenceHitList
 * Wise2_hard_link_DnaSequenceHitList
 * Wise2_DnaSequenceHitList_alloc
 * Wise2_replace_forward_DnaSequenceHitList
 * Wise2_access_forward_DnaSequenceHitList
 * Wise2_replace_backward_DnaSequenceHitList
 * Wise2_access_backward_DnaSequenceHitList
 * Wise2_free_DnaSequenceHitList [destructor]
 *
 */



/* Functions that create, manipulate or act on SegmentHitList
 *
 * Wise2_hard_link_SegmentHitList
 * Wise2_SegmentHitList_alloc_std
 * Wise2_access_seghit_SegmentHitList
 * Wise2_length_seghit_SegmentHitList
 * Wise2_flush_SegmentHitList
 * Wise2_add_SegmentHitList
 * Wise2_free_SegmentHitList [destructor]
 *
 */



/* Functions that create, manipulate or act on SegmentHit
 *
 * Wise2_hard_link_SegmentHit
 * Wise2_SegmentHit_alloc
 * Wise2_replace_name_SegmentHit
 * Wise2_access_name_SegmentHit
 * Wise2_replace_qstart_SegmentHit
 * Wise2_access_qstart_SegmentHit
 * Wise2_replace_qend_SegmentHit
 * Wise2_access_qend_SegmentHit
 * Wise2_replace_tstart_SegmentHit
 * Wise2_access_tstart_SegmentHit
 * Wise2_replace_tend_SegmentHit
 * Wise2_access_tend_SegmentHit
 * Wise2_replace_score_SegmentHit
 * Wise2_access_score_SegmentHit
 * Wise2_replace_next_hit_SegmentHit
 * Wise2_access_next_hit_SegmentHit
 * Wise2_free_SegmentHit [destructor]
 *
 */

/* API for object DnaSequenceHitList */
/* Function:  Wise2_show_DnaSequenceHitList(dsl,ofp)
 *
 * Descrip:    shows a DnaSequenceHitsList -
 *
 *             only really useful for debugging
 *
 *
 * Arg:        dsl          Undocumented argument [Wise2_DnaSequenceHitList *]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_show_DnaSequenceHitList( Wise2_DnaSequenceHitList * dsl,FILE * ofp);

/* Function:  Wise2_read_MSPcrunch_DnaSequenceHitList(ifp)
 *
 * Descrip:    Reads a MSPcrunch -x output file 
 *
 *
 * Arg:        ifp          input file to read [FILE *]
 *
 * Returns newly allocated structure [Wise2_DnaSequenceHitList *]
 *
 */
Wise2_DnaSequenceHitList * Wise2_read_MSPcrunch_DnaSequenceHitList( FILE * ifp);

/* Function:  Wise2_hard_link_DnaSequenceHitList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_DnaSequenceHitList *]
 *
 * Returns Undocumented return value [Wise2_DnaSequenceHitList *]
 *
 */
Wise2_DnaSequenceHitList * Wise2_hard_link_DnaSequenceHitList( Wise2_DnaSequenceHitList * obj);

/* Function:  Wise2_DnaSequenceHitList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_DnaSequenceHitList *]
 *
 */
Wise2_DnaSequenceHitList * Wise2_DnaSequenceHitList_alloc();

/* Function:  Wise2_replace_forward_DnaSequenceHitList(obj,forward)
 *
 * Descrip:    Replace member variable forward
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DnaSequenceHitList *]
 * Arg:        forward      New value of the variable [Wise2_SegmentHitList *]
 *
 * Returns member variable forward [boolean]
 *
 */
boolean Wise2_replace_forward_DnaSequenceHitList( Wise2_DnaSequenceHitList * obj,Wise2_SegmentHitList * forward);

/* Function:  Wise2_access_forward_DnaSequenceHitList(obj)
 *
 * Descrip:    Access member variable forward
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DnaSequenceHitList *]
 *
 * Returns member variable forward [Wise2_SegmentHitList *]
 *
 */
Wise2_SegmentHitList * Wise2_access_forward_DnaSequenceHitList( Wise2_DnaSequenceHitList * obj);

/* Function:  Wise2_replace_backward_DnaSequenceHitList(obj,backward)
 *
 * Descrip:    Replace member variable backward
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DnaSequenceHitList *]
 * Arg:        backward     New value of the variable [Wise2_SegmentHitList *]
 *
 * Returns member variable backward [boolean]
 *
 */
boolean Wise2_replace_backward_DnaSequenceHitList( Wise2_DnaSequenceHitList * obj,Wise2_SegmentHitList * backward);

/* Function:  Wise2_access_backward_DnaSequenceHitList(obj)
 *
 * Descrip:    Access member variable backward
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_DnaSequenceHitList *]
 *
 * Returns member variable backward [Wise2_SegmentHitList *]
 *
 */
Wise2_SegmentHitList * Wise2_access_backward_DnaSequenceHitList( Wise2_DnaSequenceHitList * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_DnaSequenceHitList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_DnaSequenceHitList *]
 *
 * Returns Undocumented return value [Wise2_DnaSequenceHitList *]
 *
 */
Wise2_DnaSequenceHitList * Wise2_free_DnaSequenceHitList( Wise2_DnaSequenceHitList * obj);

/* API for object SegmentHitList */
/* Function:  Wise2_hard_link_SegmentHitList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_SegmentHitList *]
 *
 * Returns Undocumented return value [Wise2_SegmentHitList *]
 *
 */
Wise2_SegmentHitList * Wise2_hard_link_SegmentHitList( Wise2_SegmentHitList * obj);

/* Function:  Wise2_SegmentHitList_alloc_std(void)
 *
 * Descrip:    Equivalent to SegmentHitList_alloc_len(SegmentHitListLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_SegmentHitList *]
 *
 */
Wise2_SegmentHitList * Wise2_SegmentHitList_alloc_std();

/* Function:  Wise2_access_seghit_SegmentHitList(obj,i)
 *
 * Descrip:    Access members stored in the seghit list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_SegmentHitList *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_SegmentHit *]
 *
 */
Wise2_SegmentHit * Wise2_access_seghit_SegmentHitList( Wise2_SegmentHitList * obj,int i);

/* Function:  Wise2_length_seghit_SegmentHitList(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_SegmentHitList *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_seghit_SegmentHitList( Wise2_SegmentHitList * obj);

/* Function:  Wise2_flush_SegmentHitList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_SegmentHitList *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_SegmentHitList( Wise2_SegmentHitList * obj);

/* Function:  Wise2_add_SegmentHitList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_SegmentHitList *]
 * Arg:        add          Object to add to the list [Wise2_SegmentHit *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SegmentHitList( Wise2_SegmentHitList * obj,Wise2_SegmentHit * add);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_SegmentHitList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_SegmentHitList *]
 *
 * Returns Undocumented return value [Wise2_SegmentHitList *]
 *
 */
Wise2_SegmentHitList * Wise2_free_SegmentHitList( Wise2_SegmentHitList * obj);

/* API for object SegmentHit */
/* Function:  Wise2_hard_link_SegmentHit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_SegmentHit *]
 *
 * Returns Undocumented return value [Wise2_SegmentHit *]
 *
 */
Wise2_SegmentHit * Wise2_hard_link_SegmentHit( Wise2_SegmentHit * obj);

/* Function:  Wise2_SegmentHit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_SegmentHit *]
 *
 */
Wise2_SegmentHit * Wise2_SegmentHit_alloc();

/* Function:  Wise2_replace_name_SegmentHit(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 * Arg:        name         New value of the variable [char *]
 *
 * Returns member variable name [boolean]
 *
 */
boolean Wise2_replace_name_SegmentHit( Wise2_SegmentHit * obj,char * name);

/* Function:  Wise2_access_name_SegmentHit(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 *
 * Returns member variable name [char *]
 *
 */
char * Wise2_access_name_SegmentHit( Wise2_SegmentHit * obj);

/* Function:  Wise2_replace_qstart_SegmentHit(obj,qstart)
 *
 * Descrip:    Replace member variable qstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 * Arg:        qstart       New value of the variable [int]
 *
 * Returns member variable qstart [boolean]
 *
 */
boolean Wise2_replace_qstart_SegmentHit( Wise2_SegmentHit * obj,int qstart);

/* Function:  Wise2_access_qstart_SegmentHit(obj)
 *
 * Descrip:    Access member variable qstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 *
 * Returns member variable qstart [int]
 *
 */
int Wise2_access_qstart_SegmentHit( Wise2_SegmentHit * obj);

/* Function:  Wise2_replace_qend_SegmentHit(obj,qend)
 *
 * Descrip:    Replace member variable qend
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 * Arg:        qend         New value of the variable [int]
 *
 * Returns member variable qend [boolean]
 *
 */
boolean Wise2_replace_qend_SegmentHit( Wise2_SegmentHit * obj,int qend);

/* Function:  Wise2_access_qend_SegmentHit(obj)
 *
 * Descrip:    Access member variable qend
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 *
 * Returns member variable qend [int]
 *
 */
int Wise2_access_qend_SegmentHit( Wise2_SegmentHit * obj);

/* Function:  Wise2_replace_tstart_SegmentHit(obj,tstart)
 *
 * Descrip:    Replace member variable tstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 * Arg:        tstart       New value of the variable [int]
 *
 * Returns member variable tstart [boolean]
 *
 */
boolean Wise2_replace_tstart_SegmentHit( Wise2_SegmentHit * obj,int tstart);

/* Function:  Wise2_access_tstart_SegmentHit(obj)
 *
 * Descrip:    Access member variable tstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 *
 * Returns member variable tstart [int]
 *
 */
int Wise2_access_tstart_SegmentHit( Wise2_SegmentHit * obj);

/* Function:  Wise2_replace_tend_SegmentHit(obj,tend)
 *
 * Descrip:    Replace member variable tend
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 * Arg:        tend         New value of the variable [int]
 *
 * Returns member variable tend [boolean]
 *
 */
boolean Wise2_replace_tend_SegmentHit( Wise2_SegmentHit * obj,int tend);

/* Function:  Wise2_access_tend_SegmentHit(obj)
 *
 * Descrip:    Access member variable tend
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 *
 * Returns member variable tend [int]
 *
 */
int Wise2_access_tend_SegmentHit( Wise2_SegmentHit * obj);

/* Function:  Wise2_replace_score_SegmentHit(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 * Arg:        score        New value of the variable [double]
 *
 * Returns member variable score [boolean]
 *
 */
boolean Wise2_replace_score_SegmentHit( Wise2_SegmentHit * obj,double score);

/* Function:  Wise2_access_score_SegmentHit(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 *
 * Returns member variable score [double]
 *
 */
double Wise2_access_score_SegmentHit( Wise2_SegmentHit * obj);

/* Function:  Wise2_replace_next_hit_SegmentHit(obj,next_hit)
 *
 * Descrip:    Replace member variable next_hit
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 * Arg:        next_hit     New value of the variable [Wise2_SegmentHit *]
 *
 * Returns member variable next_hit [boolean]
 *
 */
boolean Wise2_replace_next_hit_SegmentHit( Wise2_SegmentHit * obj,Wise2_SegmentHit * next_hit);

/* Function:  Wise2_access_next_hit_SegmentHit(obj)
 *
 * Descrip:    Access member variable next_hit
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_SegmentHit *]
 *
 * Returns member variable next_hit [Wise2_SegmentHit *]
 *
 */
Wise2_SegmentHit * Wise2_access_next_hit_SegmentHit( Wise2_SegmentHit * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_SegmentHit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_SegmentHit *]
 *
 * Returns Undocumented return value [Wise2_SegmentHit *]
 *
 */
Wise2_SegmentHit * Wise2_free_SegmentHit( Wise2_SegmentHit * obj);


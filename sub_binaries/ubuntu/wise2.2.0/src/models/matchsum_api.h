

/* Functions that create, manipulate or act on MatchSummarySet
 *
 * Wise2_MatchSummarySet_from_AlnBlock_estwise
 * Wise2_MatchSummarySet_from_AlnBlock_genewise
 * Wise2_hard_link_MatchSummarySet
 * Wise2_MatchSummarySet_alloc_std
 * Wise2_access_ms_MatchSummarySet
 * Wise2_length_ms_MatchSummarySet
 * Wise2_flush_MatchSummarySet
 * Wise2_add_MatchSummarySet
 * Wise2_free_MatchSummarySet [destructor]
 *
 */



/* Functions that create, manipulate or act on MatchSummary
 *
 * Wise2_hard_link_MatchSummary
 * Wise2_MatchSummary_alloc
 * Wise2_replace_bits_MatchSummary
 * Wise2_access_bits_MatchSummary
 * Wise2_replace_qname_MatchSummary
 * Wise2_access_qname_MatchSummary
 * Wise2_replace_tname_MatchSummary
 * Wise2_access_tname_MatchSummary
 * Wise2_replace_qstart_MatchSummary
 * Wise2_access_qstart_MatchSummary
 * Wise2_replace_qend_MatchSummary
 * Wise2_access_qend_MatchSummary
 * Wise2_replace_tstart_MatchSummary
 * Wise2_access_tstart_MatchSummary
 * Wise2_replace_tend_MatchSummary
 * Wise2_access_tend_MatchSummary
 * Wise2_replace_qintron_MatchSummary
 * Wise2_access_qintron_MatchSummary
 * Wise2_replace_qframeshift_MatchSummary
 * Wise2_access_qframeshift_MatchSummary
 * Wise2_replace_tintron_MatchSummary
 * Wise2_access_tintron_MatchSummary
 * Wise2_replace_tframeshift_MatchSummary
 * Wise2_access_tframeshift_MatchSummary
 * Wise2_free_MatchSummary [destructor]
 *
 */

/* API for object MatchSummarySet */
/* Function:  Wise2_MatchSummarySet_from_AlnBlock_estwise(alb,qname,offset,target)
 *
 * Descrip:    Builds a MatchSummarySet from a
 *             EstWise alignment. this makes
 *             alot of assumptions about the labels
 *             setc in alb, so make sure it was a 
 *             estwise alignment  - however as you
 *             can notice this is exactly the same 
 *             labels as found in genewise set
 *
 *
 * Arg:        alb          Undocumented argument [Wise2_AlnBlock *]
 * Arg:        qname        Undocumented argument [char *]
 * Arg:        offset       Undocumented argument [int]
 * Arg:        target       Undocumented argument [Wise2_Sequence *]
 *
 * Returns Undocumented return value [Wise2_MatchSummarySet *]
 *
 */
Wise2_MatchSummarySet * Wise2_MatchSummarySet_from_AlnBlock_estwise( Wise2_AlnBlock * alb,char * qname,int offset,Wise2_Sequence * target);

/* Function:  Wise2_MatchSummarySet_from_AlnBlock_genewise(alb,qname,protoff,target)
 *
 * Descrip:    Builds a MatchSummarySet from a
 *             GeneWise alignment. this makes
 *             alot of assumptions about the labels
 *             setc in alb, so make sure it was a 
 *             genewise alignment 
 *
 *
 * Arg:        alb          Undocumented argument [Wise2_AlnBlock *]
 * Arg:        qname        Undocumented argument [char *]
 * Arg:        protoff      Undocumented argument [int]
 * Arg:        target       Undocumented argument [Wise2_Sequence *]
 *
 * Returns Undocumented return value [Wise2_MatchSummarySet *]
 *
 */
Wise2_MatchSummarySet * Wise2_MatchSummarySet_from_AlnBlock_genewise( Wise2_AlnBlock * alb,char * qname,int protoff,Wise2_Sequence * target);

/* Function:  Wise2_hard_link_MatchSummarySet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_MatchSummarySet *]
 *
 * Returns Undocumented return value [Wise2_MatchSummarySet *]
 *
 */
Wise2_MatchSummarySet * Wise2_hard_link_MatchSummarySet( Wise2_MatchSummarySet * obj);

/* Function:  Wise2_MatchSummarySet_alloc_std(void)
 *
 * Descrip:    Equivalent to MatchSummarySet_alloc_len(MatchSummarySetLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_MatchSummarySet *]
 *
 */
Wise2_MatchSummarySet * Wise2_MatchSummarySet_alloc_std();

/* Function:  Wise2_access_ms_MatchSummarySet(obj,i)
 *
 * Descrip:    Access members stored in the ms list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_MatchSummarySet *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_MatchSummary *]
 *
 */
Wise2_MatchSummary * Wise2_access_ms_MatchSummarySet( Wise2_MatchSummarySet * obj,int i);

/* Function:  Wise2_length_ms_MatchSummarySet(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_MatchSummarySet *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_ms_MatchSummarySet( Wise2_MatchSummarySet * obj);

/* Function:  Wise2_flush_MatchSummarySet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_MatchSummarySet *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_MatchSummarySet( Wise2_MatchSummarySet * obj);

/* Function:  Wise2_add_MatchSummarySet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_MatchSummarySet *]
 * Arg:        add          Object to add to the list [Wise2_MatchSummary *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_MatchSummarySet( Wise2_MatchSummarySet * obj,Wise2_MatchSummary * add);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_MatchSummarySet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_MatchSummarySet *]
 *
 * Returns Undocumented return value [Wise2_MatchSummarySet *]
 *
 */
Wise2_MatchSummarySet * Wise2_free_MatchSummarySet( Wise2_MatchSummarySet * obj);

/* API for object MatchSummary */
/* Function:  Wise2_hard_link_MatchSummary(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_MatchSummary *]
 *
 * Returns Undocumented return value [Wise2_MatchSummary *]
 *
 */
Wise2_MatchSummary * Wise2_hard_link_MatchSummary( Wise2_MatchSummary * obj);

/* Function:  Wise2_MatchSummary_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_MatchSummary *]
 *
 */
Wise2_MatchSummary * Wise2_MatchSummary_alloc();

/* Function:  Wise2_replace_bits_MatchSummary(obj,bits)
 *
 * Descrip:    Replace member variable bits
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 * Arg:        bits         New value of the variable [double]
 *
 * Returns member variable bits [boolean]
 *
 */
boolean Wise2_replace_bits_MatchSummary( Wise2_MatchSummary * obj,double bits);

/* Function:  Wise2_access_bits_MatchSummary(obj)
 *
 * Descrip:    Access member variable bits
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 *
 * Returns member variable bits [double]
 *
 */
double Wise2_access_bits_MatchSummary( Wise2_MatchSummary * obj);

/* Function:  Wise2_replace_qname_MatchSummary(obj,qname)
 *
 * Descrip:    Replace member variable qname
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 * Arg:        qname        New value of the variable [char *]
 *
 * Returns member variable qname [boolean]
 *
 */
boolean Wise2_replace_qname_MatchSummary( Wise2_MatchSummary * obj,char * qname);

/* Function:  Wise2_access_qname_MatchSummary(obj)
 *
 * Descrip:    Access member variable qname
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 *
 * Returns member variable qname [char *]
 *
 */
char * Wise2_access_qname_MatchSummary( Wise2_MatchSummary * obj);

/* Function:  Wise2_replace_tname_MatchSummary(obj,tname)
 *
 * Descrip:    Replace member variable tname
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 * Arg:        tname        New value of the variable [char *]
 *
 * Returns member variable tname [boolean]
 *
 */
boolean Wise2_replace_tname_MatchSummary( Wise2_MatchSummary * obj,char * tname);

/* Function:  Wise2_access_tname_MatchSummary(obj)
 *
 * Descrip:    Access member variable tname
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 *
 * Returns member variable tname [char *]
 *
 */
char * Wise2_access_tname_MatchSummary( Wise2_MatchSummary * obj);

/* Function:  Wise2_replace_qstart_MatchSummary(obj,qstart)
 *
 * Descrip:    Replace member variable qstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 * Arg:        qstart       New value of the variable [int]
 *
 * Returns member variable qstart [boolean]
 *
 */
boolean Wise2_replace_qstart_MatchSummary( Wise2_MatchSummary * obj,int qstart);

/* Function:  Wise2_access_qstart_MatchSummary(obj)
 *
 * Descrip:    Access member variable qstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 *
 * Returns member variable qstart [int]
 *
 */
int Wise2_access_qstart_MatchSummary( Wise2_MatchSummary * obj);

/* Function:  Wise2_replace_qend_MatchSummary(obj,qend)
 *
 * Descrip:    Replace member variable qend
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 * Arg:        qend         New value of the variable [int]
 *
 * Returns member variable qend [boolean]
 *
 */
boolean Wise2_replace_qend_MatchSummary( Wise2_MatchSummary * obj,int qend);

/* Function:  Wise2_access_qend_MatchSummary(obj)
 *
 * Descrip:    Access member variable qend
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 *
 * Returns member variable qend [int]
 *
 */
int Wise2_access_qend_MatchSummary( Wise2_MatchSummary * obj);

/* Function:  Wise2_replace_tstart_MatchSummary(obj,tstart)
 *
 * Descrip:    Replace member variable tstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 * Arg:        tstart       New value of the variable [int]
 *
 * Returns member variable tstart [boolean]
 *
 */
boolean Wise2_replace_tstart_MatchSummary( Wise2_MatchSummary * obj,int tstart);

/* Function:  Wise2_access_tstart_MatchSummary(obj)
 *
 * Descrip:    Access member variable tstart
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 *
 * Returns member variable tstart [int]
 *
 */
int Wise2_access_tstart_MatchSummary( Wise2_MatchSummary * obj);

/* Function:  Wise2_replace_tend_MatchSummary(obj,tend)
 *
 * Descrip:    Replace member variable tend
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 * Arg:        tend         New value of the variable [int]
 *
 * Returns member variable tend [boolean]
 *
 */
boolean Wise2_replace_tend_MatchSummary( Wise2_MatchSummary * obj,int tend);

/* Function:  Wise2_access_tend_MatchSummary(obj)
 *
 * Descrip:    Access member variable tend
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 *
 * Returns member variable tend [int]
 *
 */
int Wise2_access_tend_MatchSummary( Wise2_MatchSummary * obj);

/* Function:  Wise2_replace_qintron_MatchSummary(obj,qintron)
 *
 * Descrip:    Replace member variable qintron
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 * Arg:        qintron      New value of the variable [int]
 *
 * Returns member variable qintron [boolean]
 *
 */
boolean Wise2_replace_qintron_MatchSummary( Wise2_MatchSummary * obj,int qintron);

/* Function:  Wise2_access_qintron_MatchSummary(obj)
 *
 * Descrip:    Access member variable qintron
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 *
 * Returns member variable qintron [int]
 *
 */
int Wise2_access_qintron_MatchSummary( Wise2_MatchSummary * obj);

/* Function:  Wise2_replace_qframeshift_MatchSummary(obj,qframeshift)
 *
 * Descrip:    Replace member variable qframeshift
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 * Arg:        qframeshift  New value of the variable [int]
 *
 * Returns member variable qframeshift [boolean]
 *
 */
boolean Wise2_replace_qframeshift_MatchSummary( Wise2_MatchSummary * obj,int qframeshift);

/* Function:  Wise2_access_qframeshift_MatchSummary(obj)
 *
 * Descrip:    Access member variable qframeshift
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 *
 * Returns member variable qframeshift [int]
 *
 */
int Wise2_access_qframeshift_MatchSummary( Wise2_MatchSummary * obj);

/* Function:  Wise2_replace_tintron_MatchSummary(obj,tintron)
 *
 * Descrip:    Replace member variable tintron
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 * Arg:        tintron      New value of the variable [int]
 *
 * Returns member variable tintron [boolean]
 *
 */
boolean Wise2_replace_tintron_MatchSummary( Wise2_MatchSummary * obj,int tintron);

/* Function:  Wise2_access_tintron_MatchSummary(obj)
 *
 * Descrip:    Access member variable tintron
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 *
 * Returns member variable tintron [int]
 *
 */
int Wise2_access_tintron_MatchSummary( Wise2_MatchSummary * obj);

/* Function:  Wise2_replace_tframeshift_MatchSummary(obj,tframeshift)
 *
 * Descrip:    Replace member variable tframeshift
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 * Arg:        tframeshift  New value of the variable [int]
 *
 * Returns member variable tframeshift [boolean]
 *
 */
boolean Wise2_replace_tframeshift_MatchSummary( Wise2_MatchSummary * obj,int tframeshift);

/* Function:  Wise2_access_tframeshift_MatchSummary(obj)
 *
 * Descrip:    Access member variable tframeshift
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_MatchSummary *]
 *
 * Returns member variable tframeshift [int]
 *
 */
int Wise2_access_tframeshift_MatchSummary( Wise2_MatchSummary * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_MatchSummary(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_MatchSummary *]
 *
 * Returns Undocumented return value [Wise2_MatchSummary *]
 *
 */
Wise2_MatchSummary * Wise2_free_MatchSummary( Wise2_MatchSummary * obj);


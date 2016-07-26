

/* Functions that create, manipulate or act on Exon
 *
 * Wise2_hard_link_Exon
 * Wise2_Exon_alloc_std
 * Wise2_replace_start_Exon
 * Wise2_access_start_Exon
 * Wise2_replace_end_Exon
 * Wise2_access_end_Exon
 * Wise2_replace_used_Exon
 * Wise2_access_used_Exon
 * Wise2_replace_score_Exon
 * Wise2_access_score_Exon
 * Wise2_access_sf_Exon
 * Wise2_length_sf_Exon
 * Wise2_flush_Exon
 * Wise2_add_Exon
 * Wise2_replace_phase_Exon
 * Wise2_access_phase_Exon
 * Wise2_free_Exon [destructor]
 *
 */



/* Functions that create, manipulate or act on Transcript
 *
 * Wise2_get_cDNA_from_Transcript
 * Wise2_hard_link_Transcript
 * Wise2_Transcript_alloc_std
 * Wise2_access_exon_Transcript
 * Wise2_length_exon_Transcript
 * Wise2_flush_ex_Transcript
 * Wise2_add_ex_Transcript
 * Wise2_replace_parent_Transcript
 * Wise2_access_parent_Transcript
 * Wise2_access_translation_Transcript
 * Wise2_length_translation_Transcript
 * Wise2_flush_Transcript
 * Wise2_add_Transcript
 * Wise2_replace_cDNA_Transcript
 * Wise2_access_cDNA_Transcript
 * Wise2_free_Transcript [destructor]
 *
 */

/* API for object Exon */
/* Function:  Wise2_hard_link_Exon(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_Exon *]
 *
 * Returns Undocumented return value [Wise2_Exon *]
 *
 */
Wise2_Exon * Wise2_hard_link_Exon( Wise2_Exon * obj);

/* Function:  Wise2_Exon_alloc_std(void)
 *
 * Descrip:    Equivalent to Exon_alloc_len(ExonLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_Exon *]
 *
 */
Wise2_Exon * Wise2_Exon_alloc_std();

/* Function:  Wise2_replace_start_Exon(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Exon *]
 * Arg:        start        New value of the variable [int]
 *
 * Returns member variable start [boolean]
 *
 */
boolean Wise2_replace_start_Exon( Wise2_Exon * obj,int start);

/* Function:  Wise2_access_start_Exon(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Exon *]
 *
 * Returns member variable start [int]
 *
 */
int Wise2_access_start_Exon( Wise2_Exon * obj);

/* Function:  Wise2_replace_end_Exon(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Exon *]
 * Arg:        end          New value of the variable [int]
 *
 * Returns member variable end [boolean]
 *
 */
boolean Wise2_replace_end_Exon( Wise2_Exon * obj,int end);

/* Function:  Wise2_access_end_Exon(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Exon *]
 *
 * Returns member variable end [int]
 *
 */
int Wise2_access_end_Exon( Wise2_Exon * obj);

/* Function:  Wise2_replace_used_Exon(obj,used)
 *
 * Descrip:    Replace member variable used
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Exon *]
 * Arg:        used         New value of the variable [boolean]
 *
 * Returns member variable used [boolean]
 *
 */
boolean Wise2_replace_used_Exon( Wise2_Exon * obj,boolean used);

/* Function:  Wise2_access_used_Exon(obj)
 *
 * Descrip:    Access member variable used
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Exon *]
 *
 * Returns member variable used [boolean]
 *
 */
boolean Wise2_access_used_Exon( Wise2_Exon * obj);

/* Function:  Wise2_replace_score_Exon(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Exon *]
 * Arg:        score        New value of the variable [double]
 *
 * Returns member variable score [boolean]
 *
 */
boolean Wise2_replace_score_Exon( Wise2_Exon * obj,double score);

/* Function:  Wise2_access_score_Exon(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Exon *]
 *
 * Returns member variable score [double]
 *
 */
double Wise2_access_score_Exon( Wise2_Exon * obj);

/* Function:  Wise2_access_sf_Exon(obj,i)
 *
 * Descrip:    Access members stored in the sf list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_Exon *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_SupportingFeature *]
 *
 */
Wise2_SupportingFeature * Wise2_access_sf_Exon( Wise2_Exon * obj,int i);

/* Function:  Wise2_length_sf_Exon(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_Exon *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_sf_Exon( Wise2_Exon * obj);

/* Function:  Wise2_flush_Exon(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_Exon *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_Exon( Wise2_Exon * obj);

/* Function:  Wise2_add_Exon(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_Exon *]
 * Arg:        add          Object to add to the list [Wise2_SupportingFeature *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Exon( Wise2_Exon * obj,Wise2_SupportingFeature * add);

/* Function:  Wise2_replace_phase_Exon(obj,phase)
 *
 * Descrip:    Replace member variable phase
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Exon *]
 * Arg:        phase        New value of the variable [int]
 *
 * Returns member variable phase [boolean]
 *
 */
boolean Wise2_replace_phase_Exon( Wise2_Exon * obj,int phase);

/* Function:  Wise2_access_phase_Exon(obj)
 *
 * Descrip:    Access member variable phase
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Exon *]
 *
 * Returns member variable phase [int]
 *
 */
int Wise2_access_phase_Exon( Wise2_Exon * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_Exon(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_Exon *]
 *
 * Returns Undocumented return value [Wise2_Exon *]
 *
 */
Wise2_Exon * Wise2_free_Exon( Wise2_Exon * obj);

/* API for object Transcript */
/* Function:  Wise2_get_cDNA_from_Transcript(trs)
 *
 * Descrip:    gets the cDNA associated with this transcript,
 *             if necessary, building it from the exon information
 *             provided.
 *
 *             returns a soft-linked object. If you want to ensure
 *             that this cDNA object remains in memory use
 *             /hard_link_cDNA on the object.
 *
 *
 * Arg:        trs          transcript to get cDNA from [Wise2_Transcript *]
 *
 * Returns cDNA of the transcript [Wise2_cDNA *]
 *
 */
Wise2_cDNA * Wise2_get_cDNA_from_Transcript( Wise2_Transcript * trs);

/* Function:  Wise2_hard_link_Transcript(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_Transcript *]
 *
 * Returns Undocumented return value [Wise2_Transcript *]
 *
 */
Wise2_Transcript * Wise2_hard_link_Transcript( Wise2_Transcript * obj);

/* Function:  Wise2_Transcript_alloc_std(void)
 *
 * Descrip:    Equivalent to Transcript_alloc_len(TranscriptLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_Transcript *]
 *
 */
Wise2_Transcript * Wise2_Transcript_alloc_std();

/* Function:  Wise2_access_exon_Transcript(obj,i)
 *
 * Descrip:    Access members stored in the exon list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_Transcript *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_Exon *]
 *
 */
Wise2_Exon * Wise2_access_exon_Transcript( Wise2_Transcript * obj,int i);

/* Function:  Wise2_length_exon_Transcript(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_Transcript *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_exon_Transcript( Wise2_Transcript * obj);

/* Function:  Wise2_flush_ex_Transcript(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_Transcript *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_ex_Transcript( Wise2_Transcript * obj);

/* Function:  Wise2_add_ex_Transcript(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_Transcript *]
 * Arg:        add          Object to add to the list [Wise2_Exon *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ex_Transcript( Wise2_Transcript * obj,Wise2_Exon * add);

/* Function:  Wise2_replace_parent_Transcript(obj,parent)
 *
 * Descrip:    Replace member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Transcript *]
 * Arg:        parent       New value of the variable [Wise2_Gene *]
 *
 * Returns member variable parent [boolean]
 *
 */
boolean Wise2_replace_parent_Transcript( Wise2_Transcript * obj,Wise2_Gene * parent);

/* Function:  Wise2_access_parent_Transcript(obj)
 *
 * Descrip:    Access member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Transcript *]
 *
 * Returns member variable parent [Wise2_Gene *]
 *
 */
Wise2_Gene * Wise2_access_parent_Transcript( Wise2_Transcript * obj);

/* Function:  Wise2_access_translation_Transcript(obj,i)
 *
 * Descrip:    Access members stored in the translation list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_Transcript *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_Translation *]
 *
 */
Wise2_Translation * Wise2_access_translation_Transcript( Wise2_Transcript * obj,int i);

/* Function:  Wise2_length_translation_Transcript(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_Transcript *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_translation_Transcript( Wise2_Transcript * obj);

/* Function:  Wise2_flush_Transcript(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_Transcript *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_Transcript( Wise2_Transcript * obj);

/* Function:  Wise2_add_Transcript(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_Transcript *]
 * Arg:        add          Object to add to the list [Wise2_Translation *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Transcript( Wise2_Transcript * obj,Wise2_Translation * add);

/* Function:  Wise2_replace_cDNA_Transcript(obj,cDNA)
 *
 * Descrip:    Replace member variable cDNA
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Transcript *]
 * Arg:        cDNA         New value of the variable [Wise2_cDNA *]
 *
 * Returns member variable cDNA [boolean]
 *
 */
boolean Wise2_replace_cDNA_Transcript( Wise2_Transcript * obj,Wise2_cDNA * cDNA);

/* Function:  Wise2_access_cDNA_Transcript(obj)
 *
 * Descrip:    Access member variable cDNA
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Transcript *]
 *
 * Returns member variable cDNA [Wise2_cDNA *]
 *
 */
Wise2_cDNA * Wise2_access_cDNA_Transcript( Wise2_Transcript * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_Transcript(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_Transcript *]
 *
 * Returns Undocumented return value [Wise2_Transcript *]
 *
 */
Wise2_Transcript * Wise2_free_Transcript( Wise2_Transcript * obj);


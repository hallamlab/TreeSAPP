

/* Functions that create, manipulate or act on Gene
 *
 * Wise2_get_Genomic_from_Gene
 * Wise2_show_pretty_Gene
 * Wise2_hard_link_Gene
 * Wise2_Gene_alloc_std
 * Wise2_replace_start_Gene
 * Wise2_access_start_Gene
 * Wise2_replace_end_Gene
 * Wise2_access_end_Gene
 * Wise2_replace_parent_Gene
 * Wise2_access_parent_Gene
 * Wise2_replace_genomic_Gene
 * Wise2_access_genomic_Gene
 * Wise2_access_transcript_Gene
 * Wise2_length_transcript_Gene
 * Wise2_flush_Gene
 * Wise2_add_Gene
 * Wise2_replace_name_Gene
 * Wise2_access_name_Gene
 * Wise2_replace_bits_Gene
 * Wise2_access_bits_Gene
 * Wise2_replace_seqname_Gene
 * Wise2_access_seqname_Gene
 * Wise2_replace_ispseudo_Gene
 * Wise2_access_ispseudo_Gene
 * Wise2_free_Gene [destructor]
 *
 */

/* API for object Gene */
/* Function:  Wise2_get_Genomic_from_Gene(gene)
 *
 * Descrip:    Gives back a Genomic sequence type
 *             from a gene.
 *
 *
 * Arg:        gene         gene to get Genomic from [Wise2_Gene *]
 *
 * Returns Genomic DNA data structure [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_get_Genomic_from_Gene( Wise2_Gene * gene);

/* Function:  Wise2_show_pretty_Gene(ge,show_supporting,ofp)
 *
 * Descrip:    Shows a gene in the biologically accepted form
 *
 *
 * Arg:        ge           Undocumented argument [Wise2_Gene *]
 * Arg:        show_supporting Undocumented argument [boolean]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Returns Undocumented return value [void]
 *
 */
void Wise2_show_pretty_Gene( Wise2_Gene * ge,boolean show_supporting,FILE * ofp);

/* Function:  Wise2_hard_link_Gene(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_Gene *]
 *
 * Returns Undocumented return value [Wise2_Gene *]
 *
 */
Wise2_Gene * Wise2_hard_link_Gene( Wise2_Gene * obj);

/* Function:  Wise2_Gene_alloc_std(void)
 *
 * Descrip:    Equivalent to Gene_alloc_len(GeneLISTLENGTH)
 *
 *
 *
 * Returns Undocumented return value [Wise2_Gene *]
 *
 */
Wise2_Gene * Wise2_Gene_alloc_std();

/* Function:  Wise2_replace_start_Gene(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 * Arg:        start        New value of the variable [int]
 *
 * Returns member variable start [boolean]
 *
 */
boolean Wise2_replace_start_Gene( Wise2_Gene * obj,int start);

/* Function:  Wise2_access_start_Gene(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 *
 * Returns member variable start [int]
 *
 */
int Wise2_access_start_Gene( Wise2_Gene * obj);

/* Function:  Wise2_replace_end_Gene(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 * Arg:        end          New value of the variable [int]
 *
 * Returns member variable end [boolean]
 *
 */
boolean Wise2_replace_end_Gene( Wise2_Gene * obj,int end);

/* Function:  Wise2_access_end_Gene(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 *
 * Returns member variable end [int]
 *
 */
int Wise2_access_end_Gene( Wise2_Gene * obj);

/* Function:  Wise2_replace_parent_Gene(obj,parent)
 *
 * Descrip:    Replace member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 * Arg:        parent       New value of the variable [Wise2_GenomicRegion *]
 *
 * Returns member variable parent [boolean]
 *
 */
boolean Wise2_replace_parent_Gene( Wise2_Gene * obj,Wise2_GenomicRegion * parent);

/* Function:  Wise2_access_parent_Gene(obj)
 *
 * Descrip:    Access member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 *
 * Returns member variable parent [Wise2_GenomicRegion *]
 *
 */
Wise2_GenomicRegion * Wise2_access_parent_Gene( Wise2_Gene * obj);

/* Function:  Wise2_replace_genomic_Gene(obj,genomic)
 *
 * Descrip:    Replace member variable genomic
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 * Arg:        genomic      New value of the variable [Wise2_Genomic *]
 *
 * Returns member variable genomic [boolean]
 *
 */
boolean Wise2_replace_genomic_Gene( Wise2_Gene * obj,Wise2_Genomic * genomic);

/* Function:  Wise2_access_genomic_Gene(obj)
 *
 * Descrip:    Access member variable genomic
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 *
 * Returns member variable genomic [Wise2_Genomic *]
 *
 */
Wise2_Genomic * Wise2_access_genomic_Gene( Wise2_Gene * obj);

/* Function:  Wise2_access_transcript_Gene(obj,i)
 *
 * Descrip:    Access members stored in the transcript list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_Gene *]
 * Arg:        i            Position in the list [int]
 *
 * Returns Element of the list [Wise2_Transcript *]
 *
 */
Wise2_Transcript * Wise2_access_transcript_Gene( Wise2_Gene * obj,int i);

/* Function:  Wise2_length_transcript_Gene(obj)
 *
 * Descrip:    discover the length of the list
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the list [Wise2_Gene *]
 *
 * Returns length of the list [int]
 *
 */
int Wise2_length_transcript_Gene( Wise2_Gene * obj);

/* Function:  Wise2_flush_Gene(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj          Object which contains the list  [Wise2_Gene *]
 *
 * Returns Undocumented return value [int]
 *
 */
int Wise2_flush_Gene( Wise2_Gene * obj);

/* Function:  Wise2_add_Gene(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj          Object which contains the list [Wise2_Gene *]
 * Arg:        add          Object to add to the list [Wise2_Transcript *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Gene( Wise2_Gene * obj,Wise2_Transcript * add);

/* Function:  Wise2_replace_name_Gene(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 * Arg:        name         New value of the variable [char *]
 *
 * Returns member variable name [boolean]
 *
 */
boolean Wise2_replace_name_Gene( Wise2_Gene * obj,char * name);

/* Function:  Wise2_access_name_Gene(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 *
 * Returns member variable name [char *]
 *
 */
char * Wise2_access_name_Gene( Wise2_Gene * obj);

/* Function:  Wise2_replace_bits_Gene(obj,bits)
 *
 * Descrip:    Replace member variable bits
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 * Arg:        bits         New value of the variable [double]
 *
 * Returns member variable bits [boolean]
 *
 */
boolean Wise2_replace_bits_Gene( Wise2_Gene * obj,double bits);

/* Function:  Wise2_access_bits_Gene(obj)
 *
 * Descrip:    Access member variable bits
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 *
 * Returns member variable bits [double]
 *
 */
double Wise2_access_bits_Gene( Wise2_Gene * obj);

/* Function:  Wise2_replace_seqname_Gene(obj,seqname)
 *
 * Descrip:    Replace member variable seqname
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 * Arg:        seqname      New value of the variable [char *]
 *
 * Returns member variable seqname [boolean]
 *
 */
boolean Wise2_replace_seqname_Gene( Wise2_Gene * obj,char * seqname);

/* Function:  Wise2_access_seqname_Gene(obj)
 *
 * Descrip:    Access member variable seqname
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 *
 * Returns member variable seqname [char *]
 *
 */
char * Wise2_access_seqname_Gene( Wise2_Gene * obj);

/* Function:  Wise2_replace_ispseudo_Gene(obj,ispseudo)
 *
 * Descrip:    Replace member variable ispseudo
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 * Arg:        ispseudo     New value of the variable [boolean]
 *
 * Returns member variable ispseudo [boolean]
 *
 */
boolean Wise2_replace_ispseudo_Gene( Wise2_Gene * obj,boolean ispseudo);

/* Function:  Wise2_access_ispseudo_Gene(obj)
 *
 * Descrip:    Access member variable ispseudo
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Gene *]
 *
 * Returns member variable ispseudo [boolean]
 *
 */
boolean Wise2_access_ispseudo_Gene( Wise2_Gene * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_Gene(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_Gene *]
 *
 * Returns Undocumented return value [Wise2_Gene *]
 *
 */
Wise2_Gene * Wise2_free_Gene( Wise2_Gene * obj);


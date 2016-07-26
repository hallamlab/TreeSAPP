

/* Helper functions in the module
 *
 * Wise2_write_pretty_str_align
 * Wise2_write_pretty_seq_align
 * Wise2_write_pretty_Protein_align
 *



/* These functions are not associated with an object */
/* Function:  Wise2_write_pretty_str_align(alb,qname,query,tname,target,name,main,ofp)
 *
 * Descrip:    This gives an interface into the alignment
 *             display using strings and files.
 *
 *
 * Arg:        alb          alignment structure [Wise2_AlnBlock *]
 * Arg:        qname        name of first sequence [char *]
 * Arg:        query        first sequence [char *]
 * Arg:        tname        name of second sequence [char *]
 * Arg:        target       second sequence [char *]
 * Arg:        name         length of the name block [int]
 * Arg:        main         length of the main block [int]
 * Arg:        ofp          output file [FILE *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_write_pretty_str_align( Wise2_AlnBlock * alb,char * qname,char * query,char * tname,char * target,int name,int main,FILE * ofp);

/* Function:  Wise2_write_pretty_seq_align(alb,q,t,name,main,ofp)
 *
 * Descrip:    This gives an interface into the alignment
 *             display using sequences and files. A more
 *             generic function is write_pretty_str_align
 *
 *
 * Arg:        alb          alignment structure [Wise2_AlnBlock *]
 * Arg:        q            first sequence [Wise2_Sequence *]
 * Arg:        t            second sequence  [Wise2_Sequence *]
 * Arg:        name         length of the name block [int]
 * Arg:        main         length of the main block [int]
 * Arg:        ofp          output file [FILE *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_write_pretty_seq_align( Wise2_AlnBlock * alb,Wise2_Sequence * q,Wise2_Sequence * t,int name,int main,FILE * ofp);

/* Function:  Wise2_write_pretty_Protein_align(alb,q,t,name,main,ofp)
 *
 * Descrip:    This gives an interface into the
 *             alignment display using Protein
 *             objects
 *
 *
 * Arg:        alb          alignment structure [Wise2_AlnBlock *]
 * Arg:        q            first sequence [Wise2_Protein *]
 * Arg:        t            second sequence  [Wise2_Protein *]
 * Arg:        name         length of the name block [int]
 * Arg:        main         length of the main block [int]
 * Arg:        ofp          output file [FILE *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_write_pretty_Protein_align( Wise2_AlnBlock * alb,Wise2_Protein * q,Wise2_Protein * t,int name,int main,FILE * ofp);


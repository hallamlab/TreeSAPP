

/* Helper functions in the module
 *
 * Wise2_protein2genomic_ascii_display
 * Wise2_protcdna_ascii_display
 *



/* These functions are not associated with an object */
/* Function:  Wise2_protein2genomic_ascii_display(alb,p,gen,ct,name,main,ofp)
 *
 * Descrip:    shows the alignment in alb between protsequence and protname
 *             with genomic into ofp with pretty formatting
 *
 *
 * Arg:        alb          logical alignment [Wise2_AlnBlock *]
 * Arg:        p            protein sequence [Wise2_Protein *]
 * Arg:        gen          genomic dna to do the comparison [Wise2_Genomic *]
 * Arg:        ct           codon table for translation [Wise2_CodonTable *]
 * Arg:        name         length of name block [int]
 * Arg:        main         length of main block [int]
 * Arg:        ofp          output file [FILE *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_protein2genomic_ascii_display( Wise2_AlnBlock * alb,Wise2_Protein * p,Wise2_Genomic * gen,Wise2_CodonTable * ct,int name,int main,FILE * ofp);

/* Function:  Wise2_protcdna_ascii_display(alb,protsequence,protname,protoff,cdna,ct,name,main,mult,ofp)
 *
 * Descrip:    shows the alignment in alb between protsequence and protname
 *             with cdna into ofp with pretty formatting
 *
 *
 * Arg:        alb          logical alignment [Wise2_AlnBlock *]
 * Arg:        protsequence protein sequence - either real or an artifical consensus [char *]
 * Arg:        protname     name of the protein [char *]
 * Arg:        protoff      offset of the alb from the protein [int]
 * Arg:        cdna         cdna of the match [Wise2_cDNA *]
 * Arg:        ct           codon table for translation [Wise2_CodonTable *]
 * Arg:        name         length of name block [int]
 * Arg:        main         length of main block [int]
 * Arg:        mult         is multi-match [boolean]
 * Arg:        ofp          output file [FILE *]
 *
 * Returns Undocumented return value [boolean]
 *
 */
boolean Wise2_protcdna_ascii_display( Wise2_AlnBlock * alb,char * protsequence,char * protname,int protoff,Wise2_cDNA * cdna,Wise2_CodonTable * ct,int name,int main,boolean mult,FILE * ofp);


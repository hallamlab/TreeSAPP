#ifndef DYNAMITEgenedisplayHEADERFILE
#define DYNAMITEgenedisplayHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "geneutil.h"



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  protein2genomic_ascii_display(alb,p,gen,ct,name,main,ofp)
 *
 * Descrip:    shows the alignment in alb between protsequence and protname
 *             with genomic into ofp with pretty formatting
 *
 *
 * Arg:         alb [UNKN ] logical alignment [AlnBlock *]
 * Arg:           p [UNKN ] protein sequence [Protein *]
 * Arg:         gen [UNKN ] genomic dna to do the comparison [Genomic *]
 * Arg:          ct [UNKN ] codon table for translation [CodonTable *]
 * Arg:        name [UNKN ] length of name block [int]
 * Arg:        main [UNKN ] length of main block [int]
 * Arg:         ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_protein2genomic_ascii_display(AlnBlock * alb,Protein * p,Genomic * gen,CodonTable * ct,int name,int main,FILE * ofp);
#define protein2genomic_ascii_display Wise2_protein2genomic_ascii_display


/* Function:  protgene_ascii_display(alb,protsequence,protname,protoff,gen,ct,name,main,mult,ofp)
 *
 * Descrip:    shows the alignment in alb between protsequence and protname
 *             with genomic into ofp with pretty formatting
 *
 *
 * Arg:                 alb [UNKN ] logical alignment [AlnBlock *]
 * Arg:        protsequence [UNKN ] protein sequence - either real or an artifical consensus [char *]
 * Arg:            protname [UNKN ] name of the protein [char *]
 * Arg:             protoff [UNKN ] offset of the alb from the protein [int]
 * Arg:                 gen [UNKN ] genomic dna to do the comparison [Genomic *]
 * Arg:                  ct [UNKN ] codon table for translation [CodonTable *]
 * Arg:                name [UNKN ] length of name block [int]
 * Arg:                main [UNKN ] length of main block [int]
 * Arg:                mult [UNKN ] is multi-match [boolean]
 * Arg:                 ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_protgene_ascii_display(AlnBlock * alb,char * protsequence,char * protname,int protoff,Genomic * gen,CodonTable * ct,int name,int main,boolean mult,FILE * ofp);
#define protgene_ascii_display Wise2_protgene_ascii_display


/* Function:  protcdna_ascii_display(alb,protsequence,protname,protoff,cdna,ct,name,main,mult,ofp)
 *
 * Descrip:    shows the alignment in alb between protsequence and protname
 *             with cdna into ofp with pretty formatting
 *
 *
 * Arg:                 alb [UNKN ] logical alignment [AlnBlock *]
 * Arg:        protsequence [UNKN ] protein sequence - either real or an artifical consensus [char *]
 * Arg:            protname [UNKN ] name of the protein [char *]
 * Arg:             protoff [UNKN ] offset of the alb from the protein [int]
 * Arg:                cdna [UNKN ] cdna of the match [cDNA *]
 * Arg:                  ct [UNKN ] codon table for translation [CodonTable *]
 * Arg:                name [UNKN ] length of name block [int]
 * Arg:                main [UNKN ] length of main block [int]
 * Arg:                mult [UNKN ] is multi-match [boolean]
 * Arg:                 ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_protcdna_ascii_display(AlnBlock * alb,char * protsequence,char * protname,int protoff,cDNA * cdna,CodonTable * ct,int name,int main,boolean mult,FILE * ofp);
#define protcdna_ascii_display Wise2_protcdna_ascii_display


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
char Wise2_match_central_line_std(char hmm,int score,char seq);
#define match_central_line_std Wise2_match_central_line_std
boolean Wise2_protdna_btc_display(AlnBlock * alb,char * protsequence,char * protname_in,int protoff,Sequence * dna,CodonTable * ct,int name,int main,btCanvas * btc,char (*match_central_line)(char,int,char),boolean multalign);
#define protdna_btc_display Wise2_protdna_btc_display
boolean Wise2_write_alignment_separator(btCanvas * btc,int aln,int score);
#define write_alignment_separator Wise2_write_alignment_separator
boolean Wise2_write_name_start_stuff(btCanvas * btc,char * protname,int protoff,char * dnaname,Sequence * dna,int name_len,AlnColumn * alc);
#define write_name_start_stuff Wise2_write_name_start_stuff
boolean Wise2_write_intron_desc(btPasteArea * btp,int start,int stop,int in_number,boolean is_split,char prot,char trans,char * dna);
#define write_intron_desc Wise2_write_intron_desc
boolean Wise2_write_3intron_match(btPasteArea * btp,int phase,int length,char * seq);
#define write_3intron_match Wise2_write_3intron_match
boolean Wise2_write_5intron_match(btPasteArea * btp,int phase,int length,char * seq);
#define write_5intron_match Wise2_write_5intron_match
boolean Wise2_write_codon_match(btPasteArea * btp,char match_letter,char midline,int c_start,char aa,char * seq) ;
#define write_codon_match Wise2_write_codon_match


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif

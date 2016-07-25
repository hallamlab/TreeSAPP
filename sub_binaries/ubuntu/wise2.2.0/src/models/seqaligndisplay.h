#ifndef DYNAMITEseqaligndisplayHEADERFILE
#define DYNAMITEseqaligndisplayHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  write_pretty_Protein_align(alb,q,t,name,main,ofp)
 *
 * Descrip:    This gives an interface into the
 *             alignment display using Protein
 *             objects
 *
 *
 * Arg:         alb [UNKN ] alignment structure [AlnBlock *]
 * Arg:           q [UNKN ] first sequence [Protein *]
 * Arg:           t [UNKN ] second sequence  [Protein *]
 * Arg:        name [UNKN ] length of the name block [int]
 * Arg:        main [UNKN ] length of the main block [int]
 * Arg:         ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_write_pretty_Protein_align(AlnBlock * alb,Protein * q,Protein * t,int name,int main,FILE * ofp);
#define write_pretty_Protein_align Wise2_write_pretty_Protein_align


/* Function:  write_pretty_seq_align(alb,q,t,name,main,ofp)
 *
 * Descrip:    This gives an interface into the alignment
 *             display using sequences and files. A more
 *             generic function is write_pretty_str_align
 *
 *
 * Arg:         alb [UNKN ] alignment structure [AlnBlock *]
 * Arg:           q [UNKN ] first sequence [Sequence *]
 * Arg:           t [UNKN ] second sequence  [Sequence *]
 * Arg:        name [UNKN ] length of the name block [int]
 * Arg:        main [UNKN ] length of the main block [int]
 * Arg:         ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_write_pretty_seq_align(AlnBlock * alb,Sequence * q,Sequence * t,int name,int main,FILE * ofp);
#define write_pretty_seq_align Wise2_write_pretty_seq_align


/* Function:  write_pretty_str_align(alb,qname,query,tname,target,name,main,ofp)
 *
 * Descrip:    This gives an interface into the alignment
 *             display using strings and files.
 *
 *
 * Arg:           alb [UNKN ] alignment structure [AlnBlock *]
 * Arg:         qname [UNKN ] name of first sequence [char *]
 * Arg:         query [UNKN ] first sequence [char *]
 * Arg:         tname [UNKN ] name of second sequence [char *]
 * Arg:        target [UNKN ] second sequence [char *]
 * Arg:          name [UNKN ] length of the name block [int]
 * Arg:          main [UNKN ] length of the main block [int]
 * Arg:           ofp [UNKN ] output file [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_write_pretty_str_align(AlnBlock * alb,char * qname,char * query,char * tname,char * target,int name,int main,FILE * ofp);
#define write_pretty_str_align Wise2_write_pretty_str_align


/* Function:  write_pretty_str_align_btc(alb,qname,query,tname,target,btc)
 *
 * Descrip:    This function writes precisely
 *             what you expect for a a simple alignment.
 *
 *             We can reuse this routine all over the place because 
 *             we dont use any hard coded structure for the
 *             query or the target sequence letters. ... but crap
 *             type checking it has to be said!
 *
 *             Also we use a generic btCanvas that could have
 *             any implementation underneath (eg, ASCII, postscript etc).
 *
 *
 * Arg:           alb [UNKN ] Undocumented argument [AlnBlock *]
 * Arg:         qname [UNKN ] Undocumented argument [char *]
 * Arg:         query [UNKN ] Undocumented argument [char *]
 * Arg:         tname [UNKN ] Undocumented argument [char *]
 * Arg:        target [UNKN ] Undocumented argument [char *]
 * Arg:           btc [UNKN ] Undocumented argument [btCanvas *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_write_pretty_str_align_btc(AlnBlock * alb,char * qname,char * query,char * tname,char * target,btCanvas * btc);
#define write_pretty_str_align_btc Wise2_write_pretty_str_align_btc


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif

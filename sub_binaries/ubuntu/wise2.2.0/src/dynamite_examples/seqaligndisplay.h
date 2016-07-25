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



/* Function:  write_pretty_seq_align(alb,q,t,name,main,ofp)
 *
 * Descrip:    This gives an interface into the alignment
 *             display using sequences and files. A more
 *             generic function is write_pretty_str_align
 *
 *
 * Arg:        alb          Undocumented argument [AlnBlock *]
 * Arg:        q            Undocumented argument [Sequence *]
 * Arg:        t            Undocumented argument [Sequence *]
 * Arg:        name         Undocumented argument [int]
 * Arg:        main         Undocumented argument [int]
 * Arg:        ofp          Undocumented argument [FILE *]
 *
 * Return    Undocumented return value [boolean]
 *
 */
boolean write_pretty_seq_align(AlnBlock * alb,Sequence * q,Sequence * t,int name,int main,FILE * ofp);


/* Function:  write_pretty_str_align(alb,qname,query,tname,target,btc)
 *
 * Descrip:    This function writes precisely
 *             what you expect for a a simple alignment.
 *
 *             We can reuse this routine all over the place because 
 *             we dont use any hard coded structure for the
 *             query or the target sequence letters.
 *
 *             Also we use a generic btCanvas that could have
 *             any implementation underneath (eg, ASCII, postscript etc).
 *
 *
 * Arg:        alb          Undocumented argument [AlnBlock *]
 * Arg:        qname        Undocumented argument [char *]
 * Arg:        query        Undocumented argument [char *]
 * Arg:        tname        Undocumented argument [char *]
 * Arg:        target       Undocumented argument [char *]
 * Arg:        btc          Undocumented argument [btCanvas *]
 *
 * Return    Undocumented return value [boolean]
 *
 */
boolean write_pretty_str_align(AlnBlock * alb,char * qname,char * query,char * tname,char * target,btCanvas * btc);


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

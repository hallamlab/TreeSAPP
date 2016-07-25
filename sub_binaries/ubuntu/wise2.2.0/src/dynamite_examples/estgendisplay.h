#ifndef DYNAMITEestgendisplayHEADERFILE
#define DYNAMITEestgendisplayHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  write_pretty_estgen_seq_align(alb,q,t,name,main,ofp)
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
boolean write_pretty_estgen_seq_align(AlnBlock * alb,Sequence * q,Sequence * t,int name,int main,FILE * ofp);


/* Function:  write_pretty_estgen_align(alb,qname,query,q_start,tname,target,t_start,name,btc)
 *
 * Descrip:    This is tied to the labels in cdna2genomic.dy 
 *
 *             Notice the use of start/end points.
 *
 *
 * Arg:        alb          Undocumented argument [AlnBlock *]
 * Arg:        qname        Undocumented argument [char *]
 * Arg:        query        Undocumented argument [char *]
 * Arg:        q_start      Undocumented argument [int]
 * Arg:        tname        Undocumented argument [char *]
 * Arg:        target       Undocumented argument [char *]
 * Arg:        t_start      Undocumented argument [int]
 * Arg:        name         Undocumented argument [int]
 * Arg:        btc          Undocumented argument [btCanvas *]
 *
 * Return    Undocumented return value [boolean]
 *
 */
boolean write_pretty_estgen_align(AlnBlock * alb,char * qname,char * query,int q_start,char * tname,char * target,int t_start,int name,btCanvas * btc);


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

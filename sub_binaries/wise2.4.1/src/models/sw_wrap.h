#ifndef DYNAMITEsw_wrapHEADERFILE
#define DYNAMITEsw_wrapHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "proteinsw.h"
#include "abc.h"
#include "pba.h"



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  Align_strings_ProteinSmithWaterman(one,two,comp,gap,ext,dpenv,dpri)
 *
 * Descrip:    This is the most *stupidly* abstracted view of two sequences
 *             getting aligned, being two strings.
 *
 *             It would be much better if you used Sequence objects or Protein
 *             objects to carry the proteins.
 *
 *
 * Arg:          one [UNKN ] string of the first sequence [char *]
 * Arg:          two [UNKN ] string of the second sequence [char *]
 * Arg:         comp [UNKN ] Comparison Matrix [CompMat *]
 * Arg:          gap [UNKN ] gap penalty [int]
 * Arg:          ext [UNKN ] extension penalty [int]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:         dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_Align_strings_ProteinSmithWaterman(char * one,char * two,CompMat * comp,int gap,int ext,DPEnvelope * dpenv,DPRunImpl * dpri);
#define Align_strings_ProteinSmithWaterman Wise2_Align_strings_ProteinSmithWaterman


/* Function:  Align_Sequences_ProteinSmithWaterman(one,two,comp,gap,ext,dpenv,dpri)
 *
 * Descrip:    This function is a mid-level abstraction of
 *             comparing two sequences, which could be
 *             generic types (eg DNA!). This is tested
 *             for and warnings are given but the alignment
 *             is still calculated. To prevent this test
 *             warning either make sure the Sequence types
 *             are set to PROTEIN or, better still, use the
 *             high level abstraction Align_Proteins_SmithWaterman
 *
 *             Otherwise this performs a standard smith waterman
 *             protein alignment...
 *
 *             To display the alignment use  write_pretty_seq_align
 *
 *
 * Arg:          one [READ ] First sequence to compare [Sequence *]
 * Arg:          two [READ ] Second sequecne to compare [Sequence *]
 * Arg:         comp [READ ] Comparison matrix to use [CompMat *]
 * Arg:          gap [UNKN ] gap penalty. Must be negative or 0 [int]
 * Arg:          ext [UNKN ] ext penalty. Must be negative or 0 [int]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:         dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [OWNER]  new AlnBlock structure representing the alignment [AlnBlock *]
 *
 */
AlnBlock * Wise2_Align_Sequences_ProteinSmithWaterman(Sequence * one,Sequence * two,CompMat * comp,int gap,int ext,DPEnvelope * dpenv,DPRunImpl * dpri);
#define Align_Sequences_ProteinSmithWaterman Wise2_Align_Sequences_ProteinSmithWaterman


/* Function:  Align_Sequences_ProteinBlockAligner(one,two,comp,bentry,bexit,b1exit,b_self,b_on,dpenv,dpri)
 *
 * Descrip:    This is a way of aligning two sequences using
 *             the ProteinBlockAligner algorithm
 *
 *
 *
 * Arg:           one [UNKN ] Sequence to align [Sequence *]
 * Arg:           two [UNKN ] Sequence to align [Sequence *]
 * Arg:          comp [UNKN ] Comparison Matrix [CompMat *]
 * Arg:        bentry [UNKN ] entry into a block [int]
 * Arg:         bexit [UNKN ] exit for a block [int]
 * Arg:        b1exit [UNKN ] exit for a block of length one [int]
 * Arg:        b_self [UNKN ] self transition [int]
 * Arg:          b_on [UNKN ] onwards transition [int]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_Align_Sequences_ProteinBlockAligner(Sequence * one,Sequence * two,CompMat * comp,int bentry,int bexit,int b1exit,int b_self,int b_on,DPEnvelope * dpenv,DPRunImpl * dpri);
#define Align_Sequences_ProteinBlockAligner Wise2_Align_Sequences_ProteinBlockAligner


/* Function:  Align_Sequences_ProteinABC(one,two,comp,a,b,c,dpenv,dpri)
 *
 * Descrip:    Align_Sequences_ProteinABC
 *             this function is analogous to Align_Sequences_ProteinSmithWaterman
 *             but using the abc model
 *
 *
 * Arg:          one [UNKN ] Sequence to align [Sequence *]
 * Arg:          two [UNKN ] Sequence to align [Sequence *]
 * Arg:         comp [UNKN ] Comparison Matrix [CompMat *]
 * Arg:            a [UNKN ] genearlized affine gap cost  [int]
 * Arg:            b [UNKN ] genearlized affine gap cost  [int]
 * Arg:            c [UNKN ] genearlized affine gap cost  [int]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:         dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_Align_Sequences_ProteinABC(Sequence * one,Sequence * two,CompMat * comp,int a, int b, int c,DPEnvelope * dpenv,DPRunImpl * dpri);
#define Align_Sequences_ProteinABC Wise2_Align_Sequences_ProteinABC


/* Function:  Align_Proteins_ABC(one,two,comp,a,b,c,dpenv,dpri)
 *
 * Descrip:    Analogous to Align_Proteins_SmithWaterman for ABC model
 *
 *
 * Arg:          one [UNKN ] protein to align [Protein *]
 * Arg:          two [UNKN ] protein to align [Protein *]
 * Arg:         comp [UNKN ] comparison matrix [CompMat *]
 * Arg:            a [UNKN ] generalized affine gap cost a [int]
 * Arg:            b [UNKN ] generalized affine gap cost b [int]
 * Arg:            c [UNKN ] generalized affine gap cost c [int]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:         dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_Align_Proteins_ABC(Protein * one,Protein * two,CompMat * comp,int a,int b,int c,DPEnvelope * dpenv,DPRunImpl * dpri) ;
#define Align_Proteins_ABC Wise2_Align_Proteins_ABC


/* Function:  Align_Proteins_SmithWaterman(one,two,comp,gap,ext,dpenv,dpri)
 *
 * Descrip:    This is the most correct way of aligning two Proteins,
 *             using Protein objects, which can be assummed to be
 *             proteins with no objections
 *
 *             To display the alignment use write_pretty_Protein_align
 *
 *
 *
 * Arg:          one [UNKN ] Protein to align [Protein *]
 * Arg:          two [UNKN ] Protein to align [Protein *]
 * Arg:         comp [UNKN ] Comparison Matrix [CompMat *]
 * Arg:          gap [UNKN ] gap penalty [int]
 * Arg:          ext [UNKN ] extension penalty [int]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:         dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_Align_Proteins_SmithWaterman(Protein * one,Protein * two,CompMat * comp,int gap,int ext,DPEnvelope * dpenv,DPRunImpl * dpri);
#define Align_Proteins_SmithWaterman Wise2_Align_Proteins_SmithWaterman


/* Function:  Hscore_from_ProteinSW(querydb,targetdb,comp,gap,ext,bits_cutoff,report_level,die_on_error,dbsi)
 *
 * Descrip:    Runs a database psw search 
 *
 *
 * Arg:             querydb [UNKN ] query database  [ProteinDB*]
 * Arg:            targetdb [UNKN ] target database [ProteinDB*]
 * Arg:                comp [UNKN ] comparison matrix [CompMat*]
 * Arg:                 gap [UNKN ] gap penalty [int]
 * Arg:                 ext [UNKN ] extension penalty [int]
 * Arg:         bits_cutoff [UNKN ]  [double]
 * Arg:        report_level [UNKN ]  [int]
 * Arg:        die_on_error [UNKN ]  [boolean]
 * Arg:                dbsi [UNKN ]  [DBSearchImpl*]
 *
 * Return [UNKN ]  Undocumented return value [Hscore *]
 *
 */
Hscore * Wise2_Hscore_from_ProteinSW(ProteinDB* querydb,ProteinDB* targetdb,CompMat* comp,int gap,int ext,double bits_cutoff,int report_level,boolean die_on_error,DBSearchImpl* dbsi);
#define Hscore_from_ProteinSW Wise2_Hscore_from_ProteinSW


/* Function:  Hscore_from_ProteinABC(querydb,targetdb,comp,a,b,c,bits_cutoff,report_level,die_on_error,dbsi)
 *
 * Descrip:    Runs a database abc search 
 *
 *
 * Arg:             querydb [UNKN ] query database  [ProteinDB*]
 * Arg:            targetdb [UNKN ] target database [ProteinDB*]
 * Arg:                comp [UNKN ] comparison matrix [CompMat*]
 * Arg:                   a [UNKN ] generalized affine gap cost a [int]
 * Arg:                   b [UNKN ] generalized affine gap cost b [int]
 * Arg:                   c [UNKN ] generalized affine gap cost c [int]
 * Arg:         bits_cutoff [UNKN ]  [double]
 * Arg:        report_level [UNKN ]  [int]
 * Arg:        die_on_error [UNKN ]  [boolean]
 * Arg:                dbsi [UNKN ]  [DBSearchImpl*]
 *
 * Return [UNKN ]  Undocumented return value [Hscore *]
 *
 */
Hscore * Wise2_Hscore_from_ProteinABC(ProteinDB* querydb,ProteinDB* targetdb,CompMat* comp,int a,int b,int c,double bits_cutoff,int report_level,boolean die_on_error,DBSearchImpl* dbsi);
#define Hscore_from_ProteinABC Wise2_Hscore_from_ProteinABC


/* Function:  Hscore_from_ProteinBA(querydb,targetdb,comp,bentry,bexit,bfor_trans,b_self_trans,b3exit,bits_cutoff,report_level,dbsi)
 *
 * Descrip:    Runs a database pba search
 *
 *
 * Arg:             querydb [UNKN ] query database [ProteinDB*]
 * Arg:            targetdb [UNKN ] target database [ProteinDB*]
 * Arg:                comp [UNKN ] comparison matrix [CompMat*]
 * Arg:              bentry [UNKN ]  [Score]
 * Arg:               bexit [UNKN ]  [Score]
 * Arg:          bfor_trans [UNKN ]  [Score]
 * Arg:        b_self_trans [UNKN ]  [Score]
 * Arg:              b3exit [UNKN ]  [Score]
 * Arg:         bits_cutoff [UNKN ]  [double]
 * Arg:        report_level [UNKN ]  [int]
 * Arg:                dbsi [UNKN ]  [DBSearchImpl*]
 *
 * Return [UNKN ]  Undocumented return value [Hscore *]
 *
 */
Hscore * Wise2_Hscore_from_ProteinBA(ProteinDB* querydb,ProteinDB* targetdb,CompMat* comp,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit,double bits_cutoff,int report_level,DBSearchImpl* dbsi);
#define Hscore_from_ProteinBA Wise2_Hscore_from_ProteinBA


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

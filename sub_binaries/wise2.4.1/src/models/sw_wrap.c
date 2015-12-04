#ifdef _cplusplus
extern "C" {
#endif
#include "sw_wrap.h"

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
# line 37 "sw_wrap.dy"
AlnBlock * Align_strings_ProteinSmithWaterman(char * one,char * two,CompMat * comp,int gap,int ext,DPEnvelope * dpenv,DPRunImpl * dpri)
{
  Sequence * one_s;
  Sequence * two_s;
  AlnBlock * out;

  /* error check the strings? */

  one_s = new_Sequence_from_strings(NULL,one);
  if( one_s == NULL ) {
    warn("Cannot make new sequence...\n");
    return NULL;
  }

  two_s = new_Sequence_from_strings(NULL,two);
  if( two_s == NULL ) {
    warn("Cannot make new sequence...\n");
    return NULL;
  }

  out = Align_Sequences_ProteinSmithWaterman(one_s,two_s,comp,gap,ext,dpenv,dpri);

  free_Sequence(one_s);
  free_Sequence(two_s);

  return out;
}
  

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
# line 88 "sw_wrap.dy"
AlnBlock * Align_Sequences_ProteinSmithWaterman(Sequence * one,Sequence * two,CompMat * comp,int gap,int ext,DPEnvelope * dpenv,DPRunImpl * dpri)
{
  AlnBlock * out = NULL;
  ComplexSequenceEvalSet * evalfunc = NULL;
  ComplexSequence * query_cs = NULL;
  ComplexSequence * target_cs = NULL;
  PackAln * pal = NULL;

  if( one == NULL || two == NULL || comp == NULL ) {
    warn("Passed in NULL objects into Align_Sequences_ProteinSmithWaterman!");
    return NULL;
  }

  if( one->type != SEQUENCE_PROTEIN ) {
    warn("Sequence %s is not typed as protein... ignoring!\n",one->name);
  }

  if( two->type != SEQUENCE_PROTEIN ) {
    warn("Sequence %s is not typed as protein... ignoring!\n",two->name);
  }


  if( gap > 0 || ext > 0 ) {
    warn("Gap penalties %d,%d only make sense if they are negative",gap,ext);
    return NULL;
  }

  evalfunc = default_aminoacid_ComplexSequenceEvalSet();
  
  query_cs = new_ComplexSequence(one,evalfunc);
  if( query_cs == NULL )
    goto cleanup;
  target_cs = new_ComplexSequence(two,evalfunc);
  if( target_cs == NULL )
    goto cleanup;

  pal = PackAln_bestmemory_ProteinSW(query_cs,target_cs,comp,gap,ext,dpenv,dpri);
  if( pal == NULL ) 
    goto cleanup;

  out = convert_PackAln_to_AlnBlock_ProteinSW(pal);
  
  goto cleanup;

  cleanup :

    if( query_cs != NULL )
      free_ComplexSequence(query_cs);

  if( target_cs != NULL )
    free_ComplexSequence(target_cs);

  if( pal != NULL )
    free_PackAln(pal);

  if( evalfunc != NULL )
    free_ComplexSequenceEvalSet(evalfunc);

  return out;

}

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
# line 164 "sw_wrap.dy"
AlnBlock * Align_Sequences_ProteinBlockAligner(Sequence * one,Sequence * two,CompMat * comp,int bentry,int bexit,int b1exit,int b_self,int b_on,DPEnvelope * dpenv,DPRunImpl * dpri)
{
  ComplexSequence * cone=NULL;
  ComplexSequence * ctwo=NULL;
  AlnBlock * alb= NULL;
  PackAln * pal=NULL;
  ComplexSequenceEvalSet * evalfunc = NULL;

  if( one == NULL || two == NULL || comp == NULL ) {
    warn("Passed in NULL objects into Align_Sequences_ProteinSmithWaterman!");
    return NULL;
  }

  if( one->type != SEQUENCE_PROTEIN ) {
    warn("Sequence %s is not typed as protein... ignoring!\n",one->name);
  }

  if( two->type != SEQUENCE_PROTEIN ) {
    warn("Sequence %s is not typed as protein... ignoring!\n",two->name);
  }


  evalfunc = default_aminoacid_ComplexSequenceEvalSet();
  
  cone = new_ComplexSequence(one,evalfunc);
  if( cone == NULL )
    goto cleanup;
  ctwo = new_ComplexSequence(two,evalfunc);
  if( ctwo == NULL )
    goto cleanup;

  pal = PackAln_bestmemory_ProteinBlockAligner(cone,ctwo,comp,bentry,b1exit,b_on,b_self,bexit,dpenv,dpri);

  if( pal == NULL ) 
    goto cleanup;

  alb = convert_PackAln_to_AlnBlock_ProteinSW(pal);
  
  goto cleanup;

  cleanup :

  if( cone != NULL )
      free_ComplexSequence(cone);

  if( ctwo != NULL )
    free_ComplexSequence(ctwo);

  if( pal != NULL )
    free_PackAln(pal);

  if( evalfunc != NULL )
    free_ComplexSequenceEvalSet(evalfunc);

  return alb;

}

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
# line 234 "sw_wrap.dy"
AlnBlock * Align_Sequences_ProteinABC(Sequence * one,Sequence * two,CompMat * comp,int a, int b, int c,DPEnvelope * dpenv,DPRunImpl * dpri)
{
   ComplexSequence * cone=NULL;
   ComplexSequence * ctwo=NULL;
   AlnBlock * alb= NULL;
   PackAln * pal=NULL;
   ComplexSequenceEvalSet * evalfunc = NULL;

   if( one == NULL || two == NULL || comp == NULL ) {
     warn("Passed in NULL objects into Align_Sequences_ProteinSmithWaterman!");
     return NULL;
   }

   if( one->type != SEQUENCE_PROTEIN ) {
     warn("Sequence %s is not typed as protein... ignoring!\n",one->name);
   }

   if( two->type != SEQUENCE_PROTEIN ) {
     warn("Sequence %s is not typed as protein... ignoring!\n",two->name);
   }

   evalfunc = default_aminoacid_ComplexSequenceEvalSet();

   cone = new_ComplexSequence(one,evalfunc);
   if( cone == NULL )
      goto cleanup;
    ctwo = new_ComplexSequence(two,evalfunc);
   if( ctwo == NULL )
     goto cleanup;

   pal = PackAln_bestmemory_abc(cone,ctwo,comp,a,b,c,NULL,dpri);

   if( pal == NULL )
     goto cleanup;

   alb = convert_PackAln_to_AlnBlock_abc(pal);

   goto cleanup;

   cleanup :

   if( cone != NULL )
       free_ComplexSequence(cone);

   if( ctwo != NULL )
     free_ComplexSequence(ctwo);

   if( pal != NULL )
     free_PackAln(pal);

   if( evalfunc != NULL )
     free_ComplexSequenceEvalSet(evalfunc);

   return alb;

}
  
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
# line 301 "sw_wrap.dy"
AlnBlock * Align_Proteins_ABC(Protein * one,Protein * two,CompMat * comp,int a,int b,int c,DPEnvelope * dpenv,DPRunImpl * dpri) 
{

  if( one == NULL || two == NULL || comp == NULL ) {
    warn("Passed in NULL objects into Align_Proteins_ABC!");
    return NULL;
  }

  return Align_Sequences_ProteinABC(one->baseseq,two->baseseq,comp,a,b,c,dpenv,dpri);
}

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
# line 326 "sw_wrap.dy"
AlnBlock * Align_Proteins_SmithWaterman(Protein * one,Protein * two,CompMat * comp,int gap,int ext,DPEnvelope * dpenv,DPRunImpl * dpri)
{
  if( one == NULL || two == NULL || comp == NULL ) {
    warn("Passed in NULL objects into Align_Proteins_SmithWaterman!");
    return NULL;
  }
    
  
  return Align_Sequences_ProteinSmithWaterman(one->baseseq,two->baseseq,comp,gap,ext,dpenv,dpri);
}

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
# line 350 "sw_wrap.dy"
Hscore * Hscore_from_ProteinSW(ProteinDB* querydb,ProteinDB* targetdb,CompMat* comp,int gap,int ext,double bits_cutoff,int report_level,boolean die_on_error,DBSearchImpl* dbsi)
{

   Hscore * out = NULL;
 
   Search_Return_Type ret;

   ret = SEARCH_ERROR;

   out = std_score_Hscore(bits_cutoff,report_level);

   if( dbsi == NULL ) {
      warn("Passed a NULL dbsi search implementaion object. Exiting without searching");
      goto exit;
   }

   if( querydb == NULL ) {
      warn("Passed a NULL querydb. Exiting without searching");
      goto exit;
   }

   if( targetdb == NULL ) {
      warn("Passed a NULL targetdb. Exiting without searching");
      goto exit;
   }

   if( comp == NULL ) {
      warn("Passed a NULL comparison matrix. Exiting without searching");
      goto exit;
   }

   ret = search_ProteinSW(dbsi,out,querydb,targetdb,comp,gap,ext);

   exit:

   return out;
}

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
# line 402 "sw_wrap.dy"
Hscore * Hscore_from_ProteinABC(ProteinDB* querydb,ProteinDB* targetdb,CompMat* comp,int a,int b,int c,double bits_cutoff,int report_level,boolean die_on_error,DBSearchImpl* dbsi)
{

   Hscore * out = NULL;
 
   Search_Return_Type ret;

   ret = SEARCH_ERROR;

   out = std_score_Hscore(bits_cutoff,report_level);

   if( dbsi == NULL ) {
      warn("Passed a NULL dbsi search implementaion object. Exiting without searching");
      goto exit;
   }

   if( querydb == NULL ) {
      warn("Passed a NULL querydb. Exiting without searching");
      goto exit;
   }

   if( targetdb == NULL ) {
      warn("Passed a NULL targetdb. Exiting without searching");
      goto exit;
   }

   if( comp == NULL ) {
      warn("Passed a NULL comparison matrix. Exiting without searching");
      goto exit;
   }

   ret = search_abc(dbsi,out,querydb,targetdb,comp,a,b,c);

   exit:

   return out;
}

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
# line 455 "sw_wrap.dy"
Hscore * Hscore_from_ProteinBA(ProteinDB* querydb,ProteinDB* targetdb,CompMat* comp,Score bentry,Score bexit,Score bfor_trans,Score b_self_trans,Score b3exit,double bits_cutoff,int report_level,DBSearchImpl* dbsi)
{

   Hscore * out = NULL;

   Search_Return_Type ret;

   ret = SEARCH_ERROR;

   out = std_score_Hscore(bits_cutoff,report_level);

   if( dbsi == NULL ) {
      warn("Passed a NULL dbsi search implementaion object. Exiting without searching");
      goto exit;
   }

   if( querydb == NULL ) {
      warn("Passed a NULL querydb. Exiting without searching");
      goto exit;
   }

   if( targetdb == NULL ) {
      warn("Passed a NULL targetdb. Exiting without searching");
      goto exit;
   }

   if( comp == NULL ) {
      warn("Passed a NULL comparison matrix. Exiting without searching");
      goto exit;
   }

   ret = search_ProteinBlockAligner(dbsi,out,querydb,targetdb,comp,bentry,bexit,bfor_trans,b_self_trans,b3exit);

   exit:

   return out;
}

# line 521 "sw_wrap.c"

#ifdef _cplusplus
}
#endif

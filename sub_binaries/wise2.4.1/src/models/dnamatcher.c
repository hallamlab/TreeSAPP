#ifdef _cplusplus
extern "C" {
#endif
#include "dnamatcher.h"


# line 20 "dnamatcher.dy"
void show_help_DnaMatchPara(FILE * ofp)
{
  fprintf(ofp,"Dna Matching Parameters\n");
  fprintf(ofp,"   -dm_match    [5]   match score\n");
  fprintf(ofp,"   -dm_mismatch [-4]  mismatch score\n");
  fprintf(ofp,"   -dm_gappen   [5]   gap open penalty\n");
  fprintf(ofp,"   -dm_extpen   [1]   gap extension penalty\n");
  
  show_help_DPRunImpl(ofp);
}

# line 31 "dnamatcher.dy"
DnaMatchPara * new_DnaMatchPara_from_argv(int * argc,char ** argv)
{
  DnaMatchPara * out;
  int match    = 5;
  int mismatch = -10;

  out = DnaMatchPara_alloc();

  strip_out_integer_argument(argc,argv,"dm_match",&match);
  strip_out_integer_argument(argc,argv,"dm_mismatch",&mismatch);

  assert(mismatch < 0 );

  out->dpri = new_DPRunImpl_from_argv(argc,argv);
  out->mat = identity_DnaMatrix(match,mismatch);
  out->dse = DnaStartEnd_from_policy("local");

  out->gap = 30;
  out->ext = 20;

  strip_out_integer_argument(argc,argv,"dm_gappen",&out->gap);
  strip_out_integer_argument(argc,argv,"dm_extpen",&out->ext);


  return out;
}


# line 59 "dnamatcher.dy"
HitList * HitList_from_Sequence_SequenceSet_DNA(Sequence * query,SequenceSet * set,DnaMatchPara * p)
{
  int i;
  HitList * out;
  HitPair * pair;
  HitAln * aln;

  AlnBlock * forward;
  AlnBlock * reverse;

  Sequence * rev;
  
  char buffer[512];

  out = HitList_alloc_std();

  for(i=0;i<set->len;i++) {
    
    rev = reverse_complement_Sequence(set->set[i]);

    ckfree(rev->name);
    sprintf(buffer,"%s.reverse",set->set[i]->name);
    rev->name  = stringalloc(buffer);

    pair = HitPair_alloc_std();

    aln = HitAln_alloc();

    forward = make_align_dnaalign(query,set->set[i],p->mat,p->dse,-p->gap,-p->ext,-p->gap,-p->ext,p->dpri);
    reverse = make_align_dnaalign(query,rev,p->mat,p->dse,-p->gap,-p->ext,-p->gap,-p->ext,p->dpri);

    if( forward->score > reverse->score ) {
      pair->query  = hard_link_Sequence(query);
      pair->target = hard_link_Sequence(set->set[i]);
      aln->alb = hard_link_AlnBlock(forward);
    } else {
      pair->query  = hard_link_Sequence(query);
      pair->target = hard_link_Sequence(rev);
      aln->alb = hard_link_AlnBlock(reverse);
    }



    add_HitPair(pair,aln);

    aln->raw_score = pair->raw_score = aln->alb->score;
    add_HitList(out,pair);


    free_AlnBlock(forward);
    free_AlnBlock(reverse);

    free_Sequence(rev);

  }

  return out;
}
# line 106 "dnamatcher.c"
/* Function:  hard_link_DnaMatchPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaMatchPara *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatchPara *]
 *
 */
DnaMatchPara * hard_link_DnaMatchPara(DnaMatchPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaMatchPara object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaMatchPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaMatchPara *]
 *
 */
DnaMatchPara * DnaMatchPara_alloc(void) 
{
    DnaMatchPara * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaMatchPara *) ckalloc (sizeof(DnaMatchPara))) == NULL)    {  
      warn("DnaMatchPara_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->dpri = NULL;    
    out->mat = NULL; 
    out->dse = NULL; 
    out->gap = 0;    
    out->ext = 0;    


    return out;  
}    


/* Function:  free_DnaMatchPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaMatchPara *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatchPara *]
 *
 */
DnaMatchPara * free_DnaMatchPara(DnaMatchPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaMatchPara obj. Should be trappable");  
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->dpri != NULL)   
      free_DPRunImpl(obj->dpri);     
    if( obj->mat != NULL)    
      free_DnaMatrix(obj->mat);  
    if( obj->dse != NULL)    
      free_DnaStartEnd(obj->dse);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

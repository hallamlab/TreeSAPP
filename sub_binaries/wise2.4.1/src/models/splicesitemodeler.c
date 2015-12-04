#ifdef _cplusplus
extern "C" {
#endif
#include "splicesitemodeler.h"




# line 33 "splicesitemodeler.dy"
int SpliceSiteModeler_ComplexSequence_eval_func(int type,void * data,char * seq)
{
  if( type != SPLICESITEMODELER_TYPE ) {
    warn("In evaluating a ComplexSequence Eval function using a SpliceSite modeler, the type was wrong. This is a bad bug, indicating a bad data problem which will probably repeat");
  }

  return SpliceSiteModel_score((SpliceSiteModel *) data,seq);
}


# line 43 "splicesitemodeler.dy"
ComplexSequenceEval * ComplexSequenceEval_from_SpliceSiteModel(SpliceSiteModel * ssm)
{
  ComplexSequenceEval * out;

  out = ComplexSequenceEval_alloc();

  /* shouldn't really add ones, but this is ok anyway. 
     Yukky hack due to not understanding a bug in the window 
     determination
     */

  /*  out->left_window  = ssm->offset + ssm->pre_splice_site +1; */
  out->left_window  = 11;
  /*  out->right_window = ssm->offset + ssm->post_splice_site +1; */
  out->right_window = 8;
  out->outside_score= NEGI;
  out->data_type    = SPLICESITEMODELER_TYPE;
  out->data         = ssm;
  out->type         = SEQUENCE_GENOMIC;
  out->eval_func    = SpliceSiteModeler_ComplexSequence_eval_func;

  return out;
}



/***
  This *must* have been passed a *seq with
  appropiate lengths on the left and right
  otherwise there is going to be tears!

  Really this function is only to be called from
  trusted function, for example, SpliceSite_ComplexSequence_eval


***/
# line 79 "splicesitemodeler.dy"
Score SpliceSiteModel_score(SpliceSiteModel * ssm,char * seq)
{
  int len;
  int i;
  int score;
  char * be = seq;
  base b;

  /* check I have enough sequence */

  /*  fprintf(stderr,"Being passed sequence %c%c%c\n",seq[0],seq[1],seq[2]); */

  /* first calculate the CC score */

  score = score_from_ComplexConsensi(seq- ssm->offset - ssm->pre_splice_site,ssm->cc);
  
  /* now move over the random score          */
  /* random score is subtracted - ie divided */
  /* out from the model                      */

  len = ssm->start_random - ssm->stop_random +1;

  for(i=0,seq = seq - ssm->start_random+1;i<len;i++,seq++) {

    if( *seq == '\0' ) {
      warn("You are attempting to score an impossible base (%d from SS) [%s] in a splice site",(int)(seq - be),be);
      return NEGI;
    }

    b = base_from_char(*seq);
    score -= ssm->rmds->base[b];
  }

  /* this is for the possibility of errors/non splice consensus etc */

  if( score < ssm->error_pos ) 
    score = ssm->error_pos;

  return score;
}

# line 120 "splicesitemodeler.dy"
SpliceSiteModel * std_5SS_SpliceSiteModel(int offset,ComplexConsensi * cc,RandomModelDNAScore * rmds)
{
  return new_SpliceSiteModel(offset,3,7,7,0,cc,rmds,0.0000001);
}

# line 125 "splicesitemodeler.dy"
SpliceSiteModel * std_3SS_SpliceSiteModel(int offset,ComplexConsensi * cc,RandomModelDNAScore * rmds)
{

  /** mystical -1 everywhere because naturally we call
    caGcag - G as being the 3' SS

    however, if you think about it, it is the

    cacCag which is really the 3'SS (the residue after the splice site).

    Rather than enforcing these rules in the calling function, we subtract one from the offset here

    **/

  return new_SpliceSiteModel(offset-1,3,3,offset+2,offset-1,cc,rmds,0.0000001);
}

# line 142 "splicesitemodeler.dy"
SpliceSiteModel * new_SpliceSiteModel(int offset,int pre_length,int post_length,int start,int stop,ComplexConsensi * cc,RandomModelDNAScore * rmds,Probability error)
{
  SpliceSiteModel * out;

  out = SpliceSiteModel_alloc();
  if( out == NULL )
    return NULL;

  out->offset = offset;
  out->pre_splice_site = pre_length;
  out->post_splice_site = post_length;
  out->start_random = start;
  out->stop_random  = stop;
  out->cc  = hard_link_ComplexConsensi(cc);
  out->rmds = hard_link_RandomModelDNAScore(rmds);
  out->error_pos = Probability2Score(error);

  return out;
}

# line 143 "splicesitemodeler.c"
/* Function:  hard_link_SpliceSiteModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SpliceSiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteModel *]
 *
 */
SpliceSiteModel * hard_link_SpliceSiteModel(SpliceSiteModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SpliceSiteModel object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SpliceSiteModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteModel *]
 *
 */
SpliceSiteModel * SpliceSiteModel_alloc(void) 
{
    SpliceSiteModel * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SpliceSiteModel *) ckalloc (sizeof(SpliceSiteModel))) == NULL)  {  
      warn("SpliceSiteModel_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->offset = 0; 
    out->pre_splice_site = 0;    
    out->post_splice_site = 0;   
    out->start_random = 0;   
    out->stop_random = 0;    
    out->cc = NULL;  
    out->rmds = NULL;    
    out->error_pos = 0;  


    return out;  
}    


/* Function:  free_SpliceSiteModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SpliceSiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteModel *]
 *
 */
SpliceSiteModel * free_SpliceSiteModel(SpliceSiteModel * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SpliceSiteModel obj. Should be trappable");   
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
    if( obj->cc != NULL) 
      free_ComplexConsensi(obj->cc);     
    if( obj->rmds != NULL)   
      free_RandomModelDNAScore(obj->rmds);   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

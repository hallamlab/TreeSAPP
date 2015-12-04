#ifdef _cplusplus
extern "C" {
#endif
#include "geneparameter.h"

# line 58 "geneparameter.dy"
GeneWiseCodonModel * GeneWiseCodonModel_from_GeneFrequencies(double * cds,GeneConsensus * donor,GeneConsensus * acceptor)
{
  int i;
  int j,k,l;
  int codon,perm;
  int one,two,three;

  GeneWiseCodonModel * out;

  out = GeneWiseCodonModel_alloc();

  for(i=0;i<64;i++) {
    out->in_donor[i] = out->in_acceptor[i] = out->in_cds[i] = 0.0000001;

    if( cds[i] < 0.00000001 ) {
      out->in_cds[i] = 0.0000001;
    } else {
      out->in_cds[i] = cds[i];
    }
  }

  /** done cds **/

  /** for splice site, need to figure out if any entry matches, and add it to the total **/

  
  for(l=0;l<donor->len;l++) {
    one   = base_from_char(donor->gsc[l]->string[0] == '-' ? 'N' : donor->gsc[l]->string[0]);
    two   = base_from_char(donor->gsc[l]->string[1] == '-' ? 'N' : donor->gsc[l]->string[1]);
    three = base_from_char(donor->gsc[l]->string[2] == '-' ? 'N' : donor->gsc[l]->string[2]);
    codon = one*25 + two *5 + three;

    for(i=0;i<4;i++) 
      for(j=0;j<4;j++) 
	for(k=0;k<4;k++) {
	  perm = permute_possible_random_bases(codon,i,j,k);

	  /** now add this number /64 to its list **/

	  out->in_donor[perm] += donor->gsc[l]->number / 64.0;
	}

  }
  
  for(l=0;l<acceptor->len;l++) {
    one   = base_from_char(acceptor->gsc[l]->string[3] == '-' ? 'N' : acceptor->gsc[l]->string[3]);
    two   = base_from_char(acceptor->gsc[l]->string[4] == '-' ? 'N' : acceptor->gsc[l]->string[4]);
    three = base_from_char(acceptor->gsc[l]->string[5] == '-' ? 'N' : acceptor->gsc[l]->string[5]);
    codon = one*25 + two *5 + three;

    for(i=0;i<4;i++) 
      for(j=0;j<4;j++) 
	for(k=0;k<4;k++) {
	  perm = permute_possible_random_bases(codon,i,j,k);

	  /** now add this number /64 to its list **/

	  out->in_acceptor[perm] = acceptor->gsc[l]->number / 64.0;
	}

  }

  return out;
}

/* Function:  GeneParameter21_from_GeneModel(gm,ct,rnd_loop,cds_loop,rnd_to_model,link_loop,link_to_model,subs_error,indel_error)
 *
 * Descrip:    This actually makes the GeneParameter21 stuff from the
 *             new statistics
 *
 *
 * Arg:                   gm [UNKN ] Undocumented argument [GeneModel *]
 * Arg:                   ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:             rnd_loop [UNKN ] Undocumented argument [Probability]
 * Arg:             cds_loop [UNKN ] Undocumented argument [Probability]
 * Arg:         rnd_to_model [UNKN ] Undocumented argument [Probability]
 * Arg:            link_loop [UNKN ] Undocumented argument [Probability]
 * Arg:        link_to_model [UNKN ] Undocumented argument [Probability]
 * Arg:           subs_error [UNKN ] Undocumented argument [Probability]
 * Arg:          indel_error [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
# line 127 "geneparameter.dy"
GeneParameter21 * GeneParameter21_from_GeneModel(GeneModel * gm,CodonTable * ct,Probability rnd_loop,Probability cds_loop,Probability rnd_to_model,Probability link_loop,Probability link_to_model,Probability subs_error,Probability indel_error)
{
  GeneParameter21 * out;
  CodonFrequency  * cf;
  RandomModelDNAScore * rmds;
  ComplexSequenceEval * cse;
  int i;

  out = GeneParameter21_alloc_len(4);

  cf = CodonFrequence_from_raw_counts(gm->codon,ct);
  out->cm = new_CodonMapper(ct,cf);
  free_CodonFrequency(cf);

  out->ct = hard_link_CodonTable(ct);

  out->gp = std_GeneParser21();
	
  for(i=0;i<5;i++) 
    out->gp->central[i] = gm->rnd->base[i];


  GeneParser21_fold_in_RandomModelDNA(out->gp,gm->rnd);

  out->gp->transition[GP21_CDS2CDS] = cds_loop;
  out->gp->transition[GP21_CDS2RND] = (1-cds_loop);
  out->gp->transition[GP21_RND2RND] = rnd_loop;
  /*  fprintf(stderr,"Score is %f\n",out->transition[GP21_RND2RND]); */
  out->gp->transition[GP21_RND2CDS] = (1-rnd_loop-rnd_to_model);
  out->gp->transition[GP21_RND2MODEL] = rnd_to_model;
  out->gp->transition[GP21_LINK2MODEL] = link_to_model;
  out->gp->transition[GP21_LINK2LINK] = link_loop;
  out->gp->transition[GP21_LINK2RND] = (1- link_loop - link_to_model) ;



  /** build random codon stuff, for soaking up "unused" cds **/

  out->rc = RandomCodon_from_raw_CodonFrequency(gm->codon,ct);

  out->cses = new_ComplexSequenceEvalSet_from_GeneModel(gm);


  /*** errors ***/

  sprinkle_errors_over_CodonMapper(out->cm,subs_error);

  add_flat_error_probabilities_GeneParser21(out->gp,indel_error);


  fold_in_RandomModelDNA_into_RandomCodon(out->rc,gm->rnd);


  return out;
}

# line 183 "geneparameter.dy"
GeneParameter21 * GeneParameter21_from_GeneFrequency21(GeneFrequency21 * gf,CodonTable * ct,RandomModelDNA * rmd,Probability rnd_loop,Probability cds_loop,Probability rnd_to_model,Probability link_loop,Probability link_to_model)
{
  GeneParameter21 * out;
  CodonFrequency  * cf;
  SpliceSiteModel * ssm;
  ComplexConsensi * cc;
  RandomModelDNAScore * rmds;
  ComplexSequenceEval * cse;

  out = GeneParameter21_alloc_len(4);

  out->gp = GeneParser21_from_GeneFrequency21_cds(gf,rnd_loop,cds_loop,rnd_to_model,link_loop,link_to_model);
  
  /** build a codon frequency, and then from that a codon mapper **/

  cf = CodonFrequency_from_GeneFrequency21(gf,ct);

  out->cm = new_CodonMapper(ct,cf);
  out->ct = hard_link_CodonTable(ct);

  free_CodonFrequency(cf);

  out->gwcm = GeneWiseCodonModel_from_GeneFrequencies(gf->codon,gf->ss5,gf->ss3);

  /** build random codon stuff, for soaking up "unused" cds **/

  out->rc = RandomCodon_from_raw_CodonFrequency(gf->codon,ct);


  /** make a new ComplexSequenceEvalSet **/

  out->cses = ComplexSequenceEvalSet_alloc_len(6);
  out->cses->type = SEQUENCE_GENOMIC;

  /** put in the base/codon eval functions **/

  add_ComplexSequenceEvalSet(out->cses,base_number_ComplexSequenceEval());
  add_ComplexSequenceEvalSet(out->cses,codon_number_ComplexSequenceEval());


  /** make a RandomModelDNAScore **/

  rmds = RandomModelDNAScore_from_RandomModelDNA(rmd);

  /** for each splice site, build a complex consensi  **/
  /** model, then a splice site model for each offset **/
  /** then both add it to the SpliceSite model list,  **/
  /** and attach a ComplexSequenceEval to the set     **/

  /** 5'SS **/

  cc = ComplexConsensi_5SS_from_GeneFrequency(gf);

  /** only one offset, at 7 **/

  ssm = std_5SS_SpliceSiteModel(0,cc,rmds);

  cse = ComplexSequenceEval_from_SpliceSiteModel(ssm);
  cse->left_lookback = 10;

  /** add to Set **/

  add_GeneParameter21(out,ssm);

  /** add complexeval to cses **/

  add_ComplexSequenceEvalSet(out->cses,cse);
  
  /** ok, free ComplexConsensi. Remember has been hard linked in std_5SS_Splice etc **/

  free_ComplexConsensi(cc);


  /** 3'SS **/

  cc = ComplexConsensi_3SS_from_GeneFrequency(gf);


  ssm = std_3SS_SpliceSiteModel(0,cc,rmds);
  cse = ComplexSequenceEval_from_SpliceSiteModel(ssm);
  cse->left_lookback = 6;

  add_GeneParameter21(out,ssm);
  add_ComplexSequenceEvalSet(out->cses,cse);


  /** ok, we can free the complex consensi for 3'SS **/

  free_ComplexConsensi(cc);

  /** and free the randommodel DNA score, as that gets hard-linked as well **/

  free_RandomModelDNAScore(rmds);


  /* 
   * ok, here we would add the necessary repeat and coding info, 
   * but for now... add flat zeros
   */


  add_ComplexSequenceEvalSet(out->cses,flat_zero());
  add_ComplexSequenceEvalSet(out->cses,flat_zero());

  /** c'est tout **/


  return out;
}


# line 257 "geneparameter.c"
/* Function:  hard_link_GeneWiseCodonModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseCodonModel *]
 *
 */
GeneWiseCodonModel * hard_link_GeneWiseCodonModel(GeneWiseCodonModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneWiseCodonModel object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneWiseCodonModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseCodonModel *]
 *
 */
GeneWiseCodonModel * GeneWiseCodonModel_alloc(void) 
{
    GeneWiseCodonModel * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWiseCodonModel *) ckalloc (sizeof(GeneWiseCodonModel))) == NULL)    {  
      warn("GeneWiseCodonModel_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* in_donor[64] is an array: no default possible */ 
    /* in_acceptor[64] is an array: no default possible */ 
    /* in_cds[64] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_GeneWiseCodonModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseCodonModel *]
 *
 */
GeneWiseCodonModel * free_GeneWiseCodonModel(GeneWiseCodonModel * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneWiseCodonModel obj. Should be trappable");    
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_GeneParameter21(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GeneParameter21
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SpliceSiteModel **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GeneParameter21(SpliceSiteModel ** list,int i,int j)  
{
    SpliceSiteModel * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GeneParameter21(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GeneParameter21 which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SpliceSiteModel **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GeneParameter21(SpliceSiteModel ** list,int left,int right,int (*comp)(SpliceSiteModel * ,SpliceSiteModel * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GeneParameter21(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GeneParameter21 (list,++last,i);    
      }  
    swap_GeneParameter21 (list,left,last);   
    qsort_GeneParameter21(list,left,last-1,comp);    
    qsort_GeneParameter21(list,last+1,right,comp);   
}    


/* Function:  sort_GeneParameter21(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GeneParameter21
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GeneParameter21 *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GeneParameter21(GeneParameter21 * obj,int (*comp)(SpliceSiteModel *, SpliceSiteModel *)) 
{
    qsort_GeneParameter21(obj->ss,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_GeneParameter21(obj,len)
 *
 * Descrip:    Really an internal function for add_GeneParameter21
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneParameter21 *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GeneParameter21(GeneParameter21 * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GeneParameter21 called with no need");    
      return TRUE;   
      }  


    if( (obj->ss = (SpliceSiteModel ** ) ckrealloc (obj->ss,sizeof(SpliceSiteModel *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GeneParameter21, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GeneParameter21(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneParameter21 *]
 * Arg:        add [OWNER] Object to add to the list [SpliceSiteModel *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GeneParameter21(GeneParameter21 * obj,SpliceSiteModel * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GeneParameter21(obj,obj->len + GeneParameter21LISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->ss[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_GeneParameter21(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GeneParameter21 *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GeneParameter21(GeneParameter21 * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->ss[i] != NULL)    {  
        free_SpliceSiteModel(obj->ss[i]);    
        obj->ss[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GeneParameter21_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneParameter21_alloc_len(GeneParameter21LISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
GeneParameter21 * GeneParameter21_alloc_std(void) 
{
    return GeneParameter21_alloc_len(GeneParameter21LISTLENGTH); 
}    


/* Function:  GeneParameter21_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
GeneParameter21 * GeneParameter21_alloc_len(int len) 
{
    GeneParameter21 * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GeneParameter21_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->ss = (SpliceSiteModel ** ) ckcalloc (len,sizeof(SpliceSiteModel *))) == NULL)   {  
      warn("Warning, ckcalloc failed in GeneParameter21_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GeneParameter21(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneParameter21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
GeneParameter21 * hard_link_GeneParameter21(GeneParameter21 * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneParameter21 object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneParameter21_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
GeneParameter21 * GeneParameter21_alloc(void) 
{
    GeneParameter21 * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneParameter21 *) ckalloc (sizeof(GeneParameter21))) == NULL)  {  
      warn("GeneParameter21_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->gp = NULL;  
    out->cm = NULL;  
    out->cses = NULL;    
    out->ss = NULL;  
    out->len = out->maxlen = 0;  
    out->rc = NULL;  
    out->gwcm = NULL;    
    out->ct = NULL;  
    out->modelled_splice = TRUE; 
    out->gms = NULL; 


    return out;  
}    


/* Function:  free_GeneParameter21(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneParameter21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParameter21 *]
 *
 */
GeneParameter21 * free_GeneParameter21(GeneParameter21 * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneParameter21 obj. Should be trappable");   
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
    if( obj->gp != NULL) 
      free_GeneParser21(obj->gp);    
    if( obj->cm != NULL) 
      free_CodonMapper(obj->cm);     
    if( obj->cses != NULL)   
      free_ComplexSequenceEvalSet(obj->cses);    
    if( obj->ss != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->ss[i] != NULL)  
          free_SpliceSiteModel(obj->ss[i]);  
        }  
      ckfree(obj->ss);   
      }  
    if( obj->rc != NULL) 
      free_RandomCodon(obj->rc);     
    if( obj->gwcm != NULL)   
      free_GeneWiseCodonModel(obj->gwcm);    
    if( obj->ct != NULL) 
      free_CodonTable(obj->ct);  
    if( obj->gms != NULL)    
      free_GeneModel(obj->gms);  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

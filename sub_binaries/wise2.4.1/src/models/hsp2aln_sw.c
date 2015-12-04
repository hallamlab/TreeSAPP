#ifdef _cplusplus
extern "C" {
#endif
#include "hsp2aln_sw.h"


#ifdef PTHREAD

# line 56 "hsp2aln_sw.dy"
HitList * HitList_from_LinearHSPmanager_heuristic_threaded(LinearHSPmanager * lm,DPRunImpl * dpri,int thr_no,HSPset2HitPairPara * para)
{
  struct hsp2aln_thread_manager tm;
  pthread_attr_t pat;
  int i;

  pthread_attr_init(&pat);     

#ifndef __sgi /* SGI can't set system scope ... */   
#ifdef  HAS_PTHREAD_SETSCOPE 
  pthread_attr_setscope(&pat, PTHREAD_SCOPE_SYSTEM);   
#endif /* set scope */   
#endif /* sgi */ 
  /* Give thread libraries a hint that there are num of threads to run */ 
#ifdef HAS_PTHREAD_SETCONCURRENCY    
  pthread_setconcurrency(thr_no+1);    
#endif /* set concurrency */ 

  tm.input = lm;
  tm.current_pos = 0;
  tm.output = HitList_alloc_std();
  tm.dpri = dpri;
  tm.para = para;
  tm.topscore = lm->set[0]->score;

  if( pthread_mutex_init(&(tm.input_lock),NULL) != 0 ) {
    fatal("Unable to make input mutex");
  }
  if( pthread_mutex_init(&(tm.output_lock),NULL) != 0 ) {
    fatal("Unable to make output mutex");
  }
  tm.pool = ckcalloc(thr_no,sizeof(pthread_t));

  for(i=0;i<thr_no;i++)    {  
    if( pthread_create(tm.pool+i,&pat,worker_thread_LM2HitList,(void *)&tm) ) 
      fatal("Unable to create a thread [thread %d]!",i); 
  }  
  
  /* Now - wait for all the threads to exit */ 
  for(i=0;i<thr_no;i++)    {  
    if( pthread_join(tm.pool[i],NULL) != 0 )  
      fatal("Unable to join a thread!");   
  }  

  ckfree(tm.pool);

  return tm.output;
}

# line 105 "hsp2aln_sw.dy"
void * worker_thread_LM2HitList(void * p)
{
  struct hsp2aln_thread_manager * tm = (struct hsp2aln_thread_manager *)p;
  HSPset * set;
  HitPair * pair;
  DPRunImpl * thread_dpri;

  thread_dpri = clone_DPRunImpl(tm->dpri);

  while(1) {
    /* get input lock, die if at end */
    if( pthread_mutex_lock(&(tm->input_lock)) != 0 )
      fatal("bad error getting input lock");
    
    if( tm->current_pos >= tm->input->len || tm->para->no_hitalns > 0 && tm->current_pos > tm->para->no_hitalns ||
	(tm->para->best_hit == TRUE && (((double)(tm->topscore - tm->input->set[tm->current_pos]->score))*100.0/tm->topscore) > tm->para->perc_hit_dropoff)
	) {
      /* end of this thread */

      if( pthread_mutex_unlock(&(tm->input_lock))!= 0 )    
	fatal("Error in releasing input lock for ProteinSW");  
      break;
    } else {
      /* get a HSPset */
      set= tm->input->set[tm->current_pos++];
      /* release the lock */
      if( pthread_mutex_unlock(&(tm->input_lock))!= 0 )    
	fatal("Error in releasing input lock for ProteinSW");  
    }

    /* now got a set in HSPset */

    pair = HitPair_from_HSPset_heuristic(set,thread_dpri,tm->input->mat,tm->para);
    
    if( pthread_mutex_lock(&(tm->output_lock))!= 0 )   
      fatal("Error on getting output lock");

    add_HitList(tm->output,pair);

    if( pthread_mutex_unlock(&(tm->output_lock))!= 0 )   
      fatal("Error on getting output lock");
    
  }

  free_DPRunImpl(thread_dpri);

  return NULL;
}

#else
# line 155 "hsp2aln_sw.dy"
HitList * HitList_from_LinearHSPmanager_heuristic_threaded(LinearHSPmanager * lm,DPRunImpl * dpri,int thr_no,HSPset2HitPairPara * para)
{
  fatal("Not compiled with pthreads");
  return NULL;
}
#endif

 
# line 163 "hsp2aln_sw.dy"
HitList * HitList_from_LinearHSPmanager_heuristic(LinearHSPmanager * lm,DPRunImpl * dpri,HSPset2HitPairPara * para)
{
  HitList * out;
  int i;
  int topscore;
  
  out = HitList_alloc_std();
  out->mat = hard_link_CompMat(lm->mat);

  if( lm->len <= 0 ) {
    return out;
  }
 
  topscore = lm->set[0]->score;

  for(i=0;i<lm->len;i++) {


    if( para->no_hitalns > 0 && i > para->no_hitalns ) {
      break;
    }
    if(  para->best_hit == TRUE && (((double)(topscore - lm->set[i]->score))*100.0/topscore) > para->perc_hit_dropoff ) {
      break;
    }
    add_HitList(out,HitPair_from_HSPset_heuristic(lm->set[i],dpri,lm->mat,para));
  }

  return out;
}


# line 194 "hsp2aln_sw.dy"
HitPair * HitPair_from_HSPset_heuristic(HSPset * set,DPRunImpl * dpri,CompMat * mat,HSPset2HitPairPara *p)
{
  HitPair * out;
  HitAln * aln;
  int i;
  int total_score = 0;
  Hsp2AlnHelper * helper;

  out = HitPair_alloc_std();
  out->query  = hard_link_Sequence(set->hsp[0]->query);
  out->target = hard_link_Sequence(set->hsp[0]->target);

  out->query->type = SEQUENCE_PROTEIN;
  out->target->type = SEQUENCE_PROTEIN;


  helper = build_HSP2AlnHelper(set,p->hsp_width,p->hsp_length,p->poor_score,p->poor_score_factor);
  
  
  for(i=0;i<helper->len;i++) {
    if( p->no_subalns != 0 && i > p->no_subalns ) {
      break;
    }

    aln = HitAln_alloc();
  
    if( p->debug == TRUE ) {
      fprintf(stdout,"For %s to %s, DPENV is\n",out->query->name,out->target->name);
      show_DPEnvelope(helper->dpenv[i],stdout);
      fprintf(stdout,"\n-----------\n");
    }

    aln->alb = Align_Sequences_ProteinSmithWaterman(out->query,out->target,mat,-12,-2,helper->dpenv[i],dpri);
    aln->raw_score = aln->alb->score;
    total_score += aln->raw_score;
    add_HitPair(out,aln);
    break;
  }

  free_Hsp2AlnHelper(helper);
  out->raw_score = total_score;
  return out;
}


# line 239 "hsp2aln_sw.dy"
Hsp2AlnHelper * build_HSP2AlnHelper(HSPset * set,int width,int tail,int min_score,int small_factor)
{
  Hsp2AlnHelper * out;
  DPEnvelope * dpenv;
  DPUnit * dpunit;
  int i;
  int j;
  int k;
  int eaten;
  int factor = 1;


  out = Hsp2AlnHelper_alloc_std();

  sort_HSPset_by_score(set);


  for(i=0;i<set->len;i++) {
    dpunit = DPUnit_alloc();

    

    dpunit->starti = set->hsp[i]->query_start - (tail*factor);
    dpunit->startj = set->hsp[i]->target_start - (tail*factor);
    dpunit->type = DPENV_DIAG;
    dpunit->height = (width*factor);
    dpunit->length = set->hsp[i]->length + 2*(tail*factor);
    

    eaten = 0;

    for(j=0;j<out->len;j++) {
      for(k=0;k<out->dpenv[j]->len;k++) {
	if( overlap_DPUnit(out->dpenv[j]->dpu[k],dpunit) == TRUE ) {
	  add_DPEnvelope(out->dpenv[j],dpunit);
	  eaten = 1;
	  break;
	}
      }
      if( eaten == 1 ) {
	break;
      }
    }

    if( eaten == 0 ) {
      dpenv = DPEnvelope_alloc_std();
      add_DPEnvelope(dpenv,dpunit);
      add_Hsp2AlnHelper(out,dpenv);
    }

  }

  return out;
}


# line 295 "hsp2aln_sw.dy"
void show_help_HSPset2HitPairPara(FILE * ofp)
{
  fprintf(ofp,"Converting HSP sets to HitPair heuristic parameters\n");
  fprintf(ofp,"  -hsp2hit_width  [no] width around each HSP to consider\n");
  fprintf(ofp,"  -hsp2hit_length [no] length around each HSP to consider\n");
  fprintf(ofp,"  -hsp2hit_subaln [no] number of HSP subalignments to consider (disabled)\n");
  fprintf(ofp,"  -hsp2hit_hitaln [no] number of hitpairs to assess\n");
  fprintf(ofp,"  -[no]hsp2hit_best    use best-in-search truncation (default no)\n");
  fprintf(ofp,"  -hsp2hit_best_perc [10] percentage off best score taken in best-in-search truncation\n");
  fprintf(ofp,"  -[no]hsp2hit_debug   print debugging features on stdout (default no)\n");
}

# line 307 "hsp2aln_sw.dy"
HSPset2HitPairPara * new_HSPset2HitPairPara_from_argv(int * argc,char ** argv)
{
  HSPset2HitPairPara * out;

  out = HSPset2HitPairPara_alloc();

  strip_out_integer_argument(argc,argv,"hsp2hit_width",&out->hsp_width);
  strip_out_integer_argument(argc,argv,"hsp2hit_length",&out->hsp_length);
  strip_out_integer_argument(argc,argv,"hsp2hit_subaln",&out->no_subalns);
  strip_out_integer_argument(argc,argv,"hsp2hit_hitaln",&out->no_hitalns);

  strip_out_boolean_def_argument(argc,argv,"hsp2hit_best",&out->best_hit);
  strip_out_boolean_def_argument(argc,argv,"hsp2hit_debug",&out->debug);
  
  strip_out_float_argument(argc,argv,"hsp2hit_best_perc",&out->perc_hit_dropoff);
  return out;
}



# line 287 "hsp2aln_sw.c"
/* Function:  hard_link_HSPset2HitPairPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [HSPset2HitPairPara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPset2HitPairPara *]
 *
 */
HSPset2HitPairPara * hard_link_HSPset2HitPairPara(HSPset2HitPairPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a HSPset2HitPairPara object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  HSPset2HitPairPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [HSPset2HitPairPara *]
 *
 */
HSPset2HitPairPara * HSPset2HitPairPara_alloc(void) 
{
    HSPset2HitPairPara * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(HSPset2HitPairPara *) ckalloc (sizeof(HSPset2HitPairPara))) == NULL)    {  
      warn("HSPset2HitPairPara_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->hsp_width = 40; 
    out->hsp_length = 90;    
    out->no_subalns = 0; 
    out->no_hitalns = 0; 
    out->best_hit = FALSE;   
    out->perc_hit_dropoff = 10.0;    
    out->debug = FALSE;  
    out->poor_score_factor = 2;  
    out->poor_score = 100;   


    return out;  
}    


/* Function:  free_HSPset2HitPairPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [HSPset2HitPairPara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPset2HitPairPara *]
 *
 */
HSPset2HitPairPara * free_HSPset2HitPairPara(HSPset2HitPairPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a HSPset2HitPairPara obj. Should be trappable");    
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


/* Function:  swap_Hsp2AlnHelper(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Hsp2AlnHelper
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DPEnvelope **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Hsp2AlnHelper(DPEnvelope ** list,int i,int j)  
{
    DPEnvelope * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Hsp2AlnHelper(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Hsp2AlnHelper which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DPEnvelope **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Hsp2AlnHelper(DPEnvelope ** list,int left,int right,int (*comp)(DPEnvelope * ,DPEnvelope * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Hsp2AlnHelper(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Hsp2AlnHelper (list,++last,i);  
      }  
    swap_Hsp2AlnHelper (list,left,last); 
    qsort_Hsp2AlnHelper(list,left,last-1,comp);  
    qsort_Hsp2AlnHelper(list,last+1,right,comp); 
}    


/* Function:  sort_Hsp2AlnHelper(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Hsp2AlnHelper
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Hsp2AlnHelper *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Hsp2AlnHelper(Hsp2AlnHelper * obj,int (*comp)(DPEnvelope *, DPEnvelope *)) 
{
    qsort_Hsp2AlnHelper(obj->dpenv,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_Hsp2AlnHelper(obj,len)
 *
 * Descrip:    Really an internal function for add_Hsp2AlnHelper
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Hsp2AlnHelper *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Hsp2AlnHelper(Hsp2AlnHelper * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Hsp2AlnHelper called with no need");  
      return TRUE;   
      }  


    if( (obj->dpenv = (DPEnvelope ** ) ckrealloc (obj->dpenv,sizeof(DPEnvelope *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Hsp2AlnHelper, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Hsp2AlnHelper(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Hsp2AlnHelper *]
 * Arg:        add [OWNER] Object to add to the list [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Hsp2AlnHelper(Hsp2AlnHelper * obj,DPEnvelope * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Hsp2AlnHelper(obj,obj->len + Hsp2AlnHelperLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->dpenv[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_Hsp2AlnHelper(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Hsp2AlnHelper *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Hsp2AlnHelper(Hsp2AlnHelper * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->dpenv[i] != NULL) {  
        free_DPEnvelope(obj->dpenv[i]);  
        obj->dpenv[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Hsp2AlnHelper_alloc_std(void)
 *
 * Descrip:    Equivalent to Hsp2AlnHelper_alloc_len(Hsp2AlnHelperLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnHelper *]
 *
 */
Hsp2AlnHelper * Hsp2AlnHelper_alloc_std(void) 
{
    return Hsp2AlnHelper_alloc_len(Hsp2AlnHelperLISTLENGTH); 
}    


/* Function:  Hsp2AlnHelper_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnHelper *]
 *
 */
Hsp2AlnHelper * Hsp2AlnHelper_alloc_len(int len) 
{
    Hsp2AlnHelper * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Hsp2AlnHelper_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->dpenv = (DPEnvelope ** ) ckcalloc (len,sizeof(DPEnvelope *))) == NULL)  {  
      warn("Warning, ckcalloc failed in Hsp2AlnHelper_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Hsp2AlnHelper(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Hsp2AlnHelper *]
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnHelper *]
 *
 */
Hsp2AlnHelper * hard_link_Hsp2AlnHelper(Hsp2AlnHelper * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Hsp2AlnHelper object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Hsp2AlnHelper_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnHelper *]
 *
 */
Hsp2AlnHelper * Hsp2AlnHelper_alloc(void) 
{
    Hsp2AlnHelper * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Hsp2AlnHelper *) ckalloc (sizeof(Hsp2AlnHelper))) == NULL)  {  
      warn("Hsp2AlnHelper_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->dpenv = NULL;   
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_Hsp2AlnHelper(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Hsp2AlnHelper *]
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnHelper *]
 *
 */
Hsp2AlnHelper * free_Hsp2AlnHelper(Hsp2AlnHelper * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Hsp2AlnHelper obj. Should be trappable"); 
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
    if( obj->dpenv != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->dpenv[i] != NULL)   
          free_DPEnvelope(obj->dpenv[i]);    
        }  
      ckfree(obj->dpenv);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_Hsp2AlnPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Hsp2AlnPara *]
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnPara *]
 *
 */
Hsp2AlnPara * hard_link_Hsp2AlnPara(Hsp2AlnPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Hsp2AlnPara object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Hsp2AlnPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnPara *]
 *
 */
Hsp2AlnPara * Hsp2AlnPara_alloc(void) 
{
    Hsp2AlnPara * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Hsp2AlnPara *) ckalloc (sizeof(Hsp2AlnPara))) == NULL)  {  
      warn("Hsp2AlnPara_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->align_size = 0; 


    return out;  
}    


/* Function:  free_Hsp2AlnPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Hsp2AlnPara *]
 *
 * Return [UNKN ]  Undocumented return value [Hsp2AlnPara *]
 *
 */
Hsp2AlnPara * free_Hsp2AlnPara(Hsp2AlnPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Hsp2AlnPara obj. Should be trappable");   
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



#ifdef _cplusplus
}
#endif

#ifdef _cplusplus
extern "C" {
#endif
#include "hsp2hitscan.h"

#include <sys/time.h>
#include <sys/resource.h>


/* Function:  new_SeedHit(posi,posj)
 *
 * Descrip:    Makes a new SeedHit
 *
 *
 * Arg:        posi [UNKN ] Undocumented argument [int]
 * Arg:        posj [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [SeedHit *]
 *
 */
# line 40 "hsp2hitscan.dy"
SeedHit * new_SeedHit(int posi,int posj)
{
  SeedHit * out;

  out = SeedHit_alloc();
  out->diagonal = posi - posj;
  out->posi = posi;
  out->posj = posj;
  out->used = 0;

  return out;
}

/* Function:  add_SeedHit_SeedHitManager(shm,t,posi,posj,score)
 *
 * Descrip:    Adds a new position to the manger, allocating any
 *             new datastructures
 *
 *
 * Arg:          shm [UNKN ] Undocumented argument [SeedHitManager *]
 * Arg:            t [UNKN ] Undocumented argument [Sequence *]
 * Arg:         posi [UNKN ] Undocumented argument [int]
 * Arg:         posj [UNKN ] Undocumented argument [int]
 * Arg:        score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 57 "hsp2hitscan.dy"
boolean add_SeedHit_SeedHitManager(SeedHitManager * shm,Sequence * t,int posi,int posj,int score)
{
  SeedHitSet * sh = NULL;

  if( (sh = g_hash_table_lookup(shm->hash,(gpointer)t)) == NULL ) {
    sh = SeedHitSet_alloc_std();
    sh->target = t;
    g_hash_table_insert(shm->hash,(gpointer)t,sh);
    add_SeedHitManager(shm,sh);
  }

  sh->score += score;

  add_SeedHitSet(sh,new_SeedHit(posi,posj));

  return TRUE;
}



/* Function:  populate_HSP_from_SeedHitManager(query,shm,para,p)
 *
 * Descrip:    Chooses which Seeds to use for HSPs. If there are more than factor x max_results
 *             cases, restricted it to factor x max_results scored by seed hits. Then askes
 *             for two diagonal seeds within a wobble factor to start the HSP process
 *
 *
 * Arg:        query [UNKN ] Undocumented argument [Sequence *]
 * Arg:          shm [UNKN ] Undocumented argument [SeedHitManager *]
 * Arg:         para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 * Arg:            p [UNKN ] Undocumented argument [HSPScanPara *]
 *
 * Return [UNKN ]  Undocumented return value [HSPmanager *]
 *
 */
# line 82 "hsp2hitscan.dy"
HSPmanager * populate_HSP_from_SeedHitManager(Sequence * query,SeedHitManager * shm,HSPScanInterfacePara * para,HSPScanPara * p)
{
  HSPmanager * out;
  int max_depth;
  int i,j,k;

  assert(shm);


  out = new_HSPmanager(query,p->mat,p->drop_off);
  
  if( p->seed_factor * para->max_results > shm->len ) {
    max_depth = shm->len;
  } else {
    sort_SeedHitManager_by_score(shm);
    max_depth = p->seed_factor * para->max_results;
  }

  fprintf(stderr,"Looking at %d seend hits from %d\n",max_depth,shm->len);

  for(i=0;i<max_depth;i++) {
    sort_SeedHitSet_by_diagonal(shm->set[i]);
    for(j=0;j<shm->set[i]->len;j++) {
      if( shm->set[i]->sh[j]->used == 1 ) {
	continue;
      }
      for(k=j+1;k<shm->set[i]->len;k++) {
	if( shm->set[i]->sh[j]->diagonal + p->twohit_wobble < shm->set[i]->sh[k]->diagonal ) {
	  break;
	}
	/* else - we have two potential seeds */

	if( shm->set[i]->sh[j]->used != 1 ) { 
	  add_pair_HSPmanager(out,shm->set[i]->target,shm->set[i]->sh[j]->posi,shm->set[i]->sh[j]->posj);
	}
	if( shm->set[i]->sh[k]->used != 1 ) { 
	  add_pair_HSPmanager(out,shm->set[i]->target,shm->set[i]->sh[k]->posi,shm->set[i]->sh[k]->posj);
	}

	shm->set[i]->sh[j]->used = 1;
	shm->set[i]->sh[k]->used = 1;

      }
    }
  }


  return out;
}



/* Function:  new_twohit_HSPScanInterface(sli,mat,drop_off,score_cutoff)
 *
 * Descrip:    Builds a 2 hit search model for protein searches
 *
 *
 * Arg:                 sli [UNKN ] Undocumented argument [SeqLookupInterface *]
 * Arg:                 mat [UNKN ] Undocumented argument [CompMat *]
 * Arg:            drop_off [UNKN ] Undocumented argument [int]
 * Arg:        score_cutoff [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [HSPScanInterface *]
 *
 */
# line 137 "hsp2hitscan.dy"
HSPScanInterface * new_twohit_HSPScanInterface(SeqLookupInterface * sli,CompMat * mat,int drop_off,int score_cutoff)
{
  HSPScanInterface * out;
  HSPScanPara * p;

  assert(sli);
  assert(mat);
  
  out = HSPScanInterface_alloc();

  p = HSPScanPara_alloc();
  p->sli = hard_link_SeqLookupInterface(sli);
  p->mat = hard_link_CompMat(mat);
  p->drop_off = drop_off;
  p->score_cutoff = score_cutoff;

  out->data = (void*)p;
  out->free_data = simple_HSPScan_free;


  out->scan_query = one_off_two_hit_HSPscan_query_direct;


  return out;
}




/* Function:  one_off_two_hit_HSPscan_query_direct(data,seq,para)
 *
 * Descrip:    two hit approach HSPscan 
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:         seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        para [UNKN ] Undocumented argument [HSPScanInterfacePara *]
 *
 * Return [UNKN ]  Undocumented return value [LinearHSPmanager *]
 *
 */
# line 169 "hsp2hitscan.dy"
LinearHSPmanager * one_off_two_hit_HSPscan_query_direct(void * data,Sequence * seq,HSPScanInterfacePara * para)
{
  HSPmanager * hspm;
  LinearHSPmanager * out;
  HSPScanPara * p = (HSPScanPara *)data;
  char * std_aa = "ACDEFGHIKLMNPQRSTVWY";

  SeedHitManager * shm;

  char newseq[5];
  int seqno[5];
  int base[5];
  int start_base;
  ArraySeqHead * head;
  int no;
  int aa;
  int score;

  SeqLookupResultInterface * slri;
  SeqLookupClientInterface * slci;
  SeqLookupResultStruct * res = NULL;

  static struct rusage use;
  int i,k,j;


  fprintf(stderr,"Got into two hit system!\n");

  shm = new_SeedHitManager();

  assert(seq != NULL);
  assert(p   != NULL);
  assert(para != NULL);

  slci = (*p->sli->get_client)(p->sli->data);


  getrusage(RUSAGE_SELF,&use);


  fprintf(stderr,"START %d.%03du %d.%03ds \n",
	  use.ru_utime.tv_sec,
	  use.ru_utime.tv_sec/1000,
	  use.ru_stime.tv_sec,
	  use.ru_stime.tv_sec/1000
	  );


  for(i=0;i<seq->len-5;i++) {

    
    if( (*slci->is_populated)(slci->data,seq_number_aa_5mer(seq->seq+i)) ) {

      slri = (*slci->lookup)(slci->data,seq_number_aa_5mer(seq->seq+i));

      res = NULL;
      for(;(*slri->is_more)(slri->data);) {    
	res = (*slri->next)(slri->data,res);
	add_SeedHit_SeedHitManager(shm,res->seq,i,res->pos,5);

      }
      free_SeqLookupResultInterface(slri);
    }

    
    for(score=0,j=0;j<5;j++) {

      seqno[j] = base[j]*(toupper(seq->seq[i+j]-'A'));

    }

    for(j=0;j<5;j++) {
      for(aa=0;aa<20;aa++) {
	if( seq->seq[i+j] == std_aa[aa] ) {
	  continue;
	}

	seqno[j] = base[j]*(std_aa[aa]-'A');

	no= seqno[0]+seqno[1]+seqno[2]+seqno[3]+seqno[4];


	  if( (*slci->is_populated)(slci->data,seq_number_aa_5mer(newseq)) ) {

	    slri = (*slci->lookup)(slci->data,seq_number_aa_5mer(newseq));
	    res = NULL;
	    for(;(*slri->is_more)(slri->data);) {
	      res = (*slri->next)(slri->data,res);
	      add_SeedHit_SeedHitManager(shm,res->seq,i,res->pos,1);

	    }
	    free_SeqLookupResultInterface(slri);
	  }


	seqno[j]  = base[j]*(toupper(seq->seq[i+j]-'A'));
	newseq[j] = seq->seq[i+j];
      }
      
    }  
  }

  getrusage(RUSAGE_SELF,&use);

  fprintf(stderr,"END OF SEED %d.%03du %d.%03ds \n",
	  use.ru_utime.tv_sec,
	  use.ru_utime.tv_sec/1000,
	  use.ru_stime.tv_sec,
	  use.ru_stime.tv_sec/1000
	  );



  hspm = populate_HSP_from_SeedHitManager(seq,shm,para,p);

  getrusage(RUSAGE_SELF,&use);

  fprintf(stderr,"POPULATION %d.%03du %d.%03ds \n",
	  use.ru_utime.tv_sec,
	  use.ru_utime.tv_sec/1000,
	  use.ru_stime.tv_sec,
	  use.ru_stime.tv_sec/1000
	  );


  delete_SeedHitManager(shm);

  if( para->use_protein_heuristic == TRUE ) {
    out = new_LinearHSPmanager_heuristic_max(hspm,para->max_results);
  } else {
    out = new_LinearHSPmanager_flat(hspm);
  }

  fprintf(stdout,"LINEARISED %d.%03du %d.%03ds \n",
	  use.ru_utime.tv_sec,
	  use.ru_utime.tv_sec/1000,
	  use.ru_stime.tv_sec,
	  use.ru_stime.tv_sec/1000
	  );

  
  return out;

}

/* Function:  compare_SeedHitSet_score(one,two)
 *
 * Descrip:    internal for score sorting
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [SeedHitSet *]
 * Arg:        two [UNKN ] Undocumented argument [SeedHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 317 "hsp2hitscan.dy"
int compare_SeedHitSet_score(SeedHitSet * one,SeedHitSet * two)
{
  return two->score - one->score;
}

/* Function:  compare_SeedHit_diagonal(one,two)
 *
 * Descrip:    internal for diagonal sorting
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [SeedHit *]
 * Arg:        two [UNKN ] Undocumented argument [SeedHit *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 325 "hsp2hitscan.dy"
int compare_SeedHit_diagonal(SeedHit * one,SeedHit * two)
{
  return one->diagonal - two->diagonal;
}

/* Function:  sort_SeedHitManager_by_score(shm)
 *
 * Descrip:    Sorts hit managers by score (highest first)
 *
 *
 * Arg:        shm [UNKN ] Undocumented argument [SeedHitManager *]
 *
 */
# line 333 "hsp2hitscan.dy"
void sort_SeedHitManager_by_score(SeedHitManager * shm)
{
  sort_SeedHitManager(shm,compare_SeedHitSet_score);
}

/* Function:  sort_SeedHitSet_by_diagonal(sh)
 *
 * Descrip:    Sorts SeedHit by diagonal (lowest first)
 *
 *
 * Arg:        sh [UNKN ] Undocumented argument [SeedHitSet *]
 *
 */
# line 341 "hsp2hitscan.dy"
void sort_SeedHitSet_by_diagonal(SeedHitSet * sh)
{
  sort_SeedHitSet(sh,compare_SeedHit_diagonal);
}


/* Function:  new_SeedHitManager(void)
 *
 * Descrip:    Makes an empty SeedHit Manager
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
# line 350 "hsp2hitscan.dy"
SeedHitManager * new_SeedHitManager(void)
{
  SeedHitManager * out;

  out = SeedHitManager_alloc_std();
  out->hash =  g_hash_table_new(g_direct_hash,g_direct_equal);

  return out;
}


/* Function:  delete_SeedHitManager(shm)
 *
 * Descrip:    Free SeedHitManager
 *
 *
 *
 * Arg:        shm [UNKN ] Undocumented argument [SeedHitManager *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
# line 365 "hsp2hitscan.dy"
SeedHitManager * delete_SeedHitManager(SeedHitManager * shm)
{
  g_hash_table_destroy(shm->hash);

  return free_SeedHitManager(shm);
}
# line 419 "hsp2hitscan.c"
/* Function:  hard_link_SeedHit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeedHit *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHit *]
 *
 */
SeedHit * hard_link_SeedHit(SeedHit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeedHit object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeedHit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHit *]
 *
 */
SeedHit * SeedHit_alloc(void) 
{
    SeedHit * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeedHit *) ckalloc (sizeof(SeedHit))) == NULL)  {  
      warn("SeedHit_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->diagonal = 0;   
    out->posi = 0;   
    out->posj = 0;   
    out->used = 0;   


    return out;  
}    


/* Function:  free_SeedHit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeedHit *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHit *]
 *
 */
SeedHit * free_SeedHit(SeedHit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SeedHit obj. Should be trappable");   
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


/* Function:  swap_SeedHitSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SeedHitSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SeedHit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SeedHitSet(SeedHit ** list,int i,int j)  
{
    SeedHit * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SeedHitSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SeedHitSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SeedHit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SeedHitSet(SeedHit ** list,int left,int right,int (*comp)(SeedHit * ,SeedHit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SeedHitSet(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SeedHitSet (list,++last,i); 
      }  
    swap_SeedHitSet (list,left,last);    
    qsort_SeedHitSet(list,left,last-1,comp); 
    qsort_SeedHitSet(list,last+1,right,comp);    
}    


/* Function:  sort_SeedHitSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SeedHitSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SeedHitSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SeedHitSet(SeedHitSet * obj,int (*comp)(SeedHit *, SeedHit *)) 
{
    qsort_SeedHitSet(obj->sh,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_SeedHitSet(obj,len)
 *
 * Descrip:    Really an internal function for add_SeedHitSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeedHitSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SeedHitSet(SeedHitSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SeedHitSet called with no need"); 
      return TRUE;   
      }  


    if( (obj->sh = (SeedHit ** ) ckrealloc (obj->sh,sizeof(SeedHit *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_SeedHitSet, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SeedHitSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeedHitSet *]
 * Arg:        add [OWNER] Object to add to the list [SeedHit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SeedHitSet(SeedHitSet * obj,SeedHit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SeedHitSet(obj,obj->len + SeedHitSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->sh[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_SeedHitSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SeedHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SeedHitSet(SeedHitSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->sh[i] != NULL)    {  
        free_SeedHit(obj->sh[i]);    
        obj->sh[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SeedHitSet_alloc_std(void)
 *
 * Descrip:    Equivalent to SeedHitSet_alloc_len(SeedHitSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHitSet *]
 *
 */
SeedHitSet * SeedHitSet_alloc_std(void) 
{
    return SeedHitSet_alloc_len(SeedHitSetLISTLENGTH);   
}    


/* Function:  SeedHitSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitSet *]
 *
 */
SeedHitSet * SeedHitSet_alloc_len(int len) 
{
    SeedHitSet * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SeedHitSet_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->sh = (SeedHit ** ) ckcalloc (len,sizeof(SeedHit *))) == NULL)   {  
      warn("Warning, ckcalloc failed in SeedHitSet_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SeedHitSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeedHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitSet *]
 *
 */
SeedHitSet * hard_link_SeedHitSet(SeedHitSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeedHitSet object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeedHitSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHitSet *]
 *
 */
SeedHitSet * SeedHitSet_alloc(void) 
{
    SeedHitSet * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeedHitSet *) ckalloc (sizeof(SeedHitSet))) == NULL)    {  
      warn("SeedHitSet_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->sh = NULL;  
    out->len = out->maxlen = 0;  
    out->score = 0;  


    return out;  
}    


/* Function:  free_SeedHitSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeedHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitSet *]
 *
 */
SeedHitSet * free_SeedHitSet(SeedHitSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SeedHitSet obj. Should be trappable");    
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
    if( obj->sh != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->sh[i] != NULL)  
          free_SeedHit(obj->sh[i]);  
        }  
      ckfree(obj->sh);   
      }  
    /* obj->target is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_SeedHitManager(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SeedHitManager
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SeedHitSet **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SeedHitManager(SeedHitSet ** list,int i,int j)  
{
    SeedHitSet * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SeedHitManager(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SeedHitManager which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SeedHitSet **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SeedHitManager(SeedHitSet ** list,int left,int right,int (*comp)(SeedHitSet * ,SeedHitSet * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SeedHitManager(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SeedHitManager (list,++last,i); 
      }  
    swap_SeedHitManager (list,left,last);    
    qsort_SeedHitManager(list,left,last-1,comp); 
    qsort_SeedHitManager(list,last+1,right,comp);    
}    


/* Function:  sort_SeedHitManager(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SeedHitManager
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SeedHitManager *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SeedHitManager(SeedHitManager * obj,int (*comp)(SeedHitSet *, SeedHitSet *)) 
{
    qsort_SeedHitManager(obj->set,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_SeedHitManager(obj,len)
 *
 * Descrip:    Really an internal function for add_SeedHitManager
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeedHitManager *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SeedHitManager(SeedHitManager * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SeedHitManager called with no need"); 
      return TRUE;   
      }  


    if( (obj->set = (SeedHitSet ** ) ckrealloc (obj->set,sizeof(SeedHitSet *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_SeedHitManager, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SeedHitManager(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SeedHitManager *]
 * Arg:        add [OWNER] Object to add to the list [SeedHitSet *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SeedHitManager(SeedHitManager * obj,SeedHitSet * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SeedHitManager(obj,obj->len + SeedHitManagerLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->set[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_SeedHitManager(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SeedHitManager *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SeedHitManager(SeedHitManager * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->set[i] != NULL)   {  
        free_SeedHitSet(obj->set[i]);    
        obj->set[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SeedHitManager_alloc_std(void)
 *
 * Descrip:    Equivalent to SeedHitManager_alloc_len(SeedHitManagerLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * SeedHitManager_alloc_std(void) 
{
    return SeedHitManager_alloc_len(SeedHitManagerLISTLENGTH);   
}    


/* Function:  SeedHitManager_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * SeedHitManager_alloc_len(int len) 
{
    SeedHitManager * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SeedHitManager_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->set = (SeedHitSet ** ) ckcalloc (len,sizeof(SeedHitSet *))) == NULL)    {  
      warn("Warning, ckcalloc failed in SeedHitManager_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SeedHitManager(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SeedHitManager *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * hard_link_SeedHitManager(SeedHitManager * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SeedHitManager object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SeedHitManager_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * SeedHitManager_alloc(void) 
{
    SeedHitManager * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SeedHitManager *) ckalloc (sizeof(SeedHitManager))) == NULL)    {  
      warn("SeedHitManager_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->set = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_SeedHitManager(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SeedHitManager *]
 *
 * Return [UNKN ]  Undocumented return value [SeedHitManager *]
 *
 */
SeedHitManager * free_SeedHitManager(SeedHitManager * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SeedHitManager obj. Should be trappable");    
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
    /* obj->hash is linked in */ 
    if( obj->set != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->set[i] != NULL) 
          free_SeedHitSet(obj->set[i]);  
        }  
      ckfree(obj->set);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

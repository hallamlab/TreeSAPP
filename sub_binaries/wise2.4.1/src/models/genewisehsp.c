#ifdef _cplusplus
extern "C" {
#endif
#include "genewisehsp.h"

# line 36 "genewisehsp.dy"
GeneWiseRunPara * new_GeneWiseRunPara_from_argv(int * argc,char ** argv)
{
  GeneWiseRunPara * out;

  out = GeneWiseRunPara_alloc();
  out->use_hsp = FALSE;
  out->edge_query = 10;
  out->edge_target = 3000;
  out->splice_spread = 5;

  strip_out_boolean_def_argument(argc,argv,"gwhsp",&out->use_hsp);
  strip_out_boolean_def_argument(argc,argv,"gwdebug",&out->debug);
  strip_out_integer_argument(argc,argv,"gw_edgequery",&out->edge_query);
  strip_out_integer_argument(argc,argv,"gw_edgetarget",&out->edge_target);
  strip_out_integer_argument(argc,argv,"gw_splicespread",&out->splice_spread);

  return out;
}

# line 55 "genewisehsp.dy"
void show_help_GeneWiseRunPara(FILE * ofp)
{
  fprintf(ofp,"Genewise protein running heuristics\n");
  fprintf(ofp,"   -[no]gwhsp        use heuristics for proteins [FALSE currently]\n");
  fprintf(ofp,"   -gw_edgequery     at start/end, amount of protein area to expand [10]\n");
  fprintf(ofp,"   -gw_edgetarget    at start/end, amount of DNA area to expand [3000]\n");
  fprintf(ofp,"   -gw_splicespread  spread around splice sites in codons [5]\n");
  fprintf(ofp,"   -gwdebug          print out debugging of heuristics on stdout\n");

}
    


# line 68 "genewisehsp.dy"
DPEnvelope * DPEnvelope_from_protein_gen(Sequence * prot,Sequence * dna,CompMat * mat,CodonTable * ct,GeneWiseRunPara *p)
{
  SeqLookupInterface * sli;
  HSPScanInterface * hsi;
  HSPScanInterfacePara para;
  Sequence * trans;
  int i;
  int j;
  int frame;
  char * temp_seq;
  int can_use = 0;
  
  LinearHSPmanager * lm;
  GeneWiseHSPmanager * gwh;
  GeneWiseHSPmanager * final;

  DPEnvelope * out;
  DPUnit * dpunit;

  int prev_query;
  int prev_target;

  SeqLookupLoadPara loadpara;


  info("Using HSP based heuristic for bounds calculation. Result may not be accurate (in particular for low similarity)");

  loadpara.tile_freq = 1;

  sli = new_ghash_SeqLookupInterface();
  (*sli->add_seq)(sli->data,prot,&loadpara);

  hsi = new_one_off_HSPScanInterface(sli,mat,15,5);

  para.min_score= 10;
  para.max_results = 750;
  para.use_protein_heuristic = FALSE;
 

  gwh = GeneWiseHSPmanager_alloc_std();

  for(frame=0;frame<3;frame++) {
    temp_seq = calloc(1+dna->len/3,sizeof(char));
    for(j=0,i=frame;i+3 < dna->len;i+=3,j++) {
      temp_seq[j] = aminoacid_from_seq(ct,dna->seq+i);
    }
    temp_seq[j]= '\0';
    trans = Sequence_alloc();
    trans->name = stringalloc("temp_seq");
    trans->seq  = temp_seq;
    trans->len  = strlen(temp_seq);

    if( p->debug ) {
      info("starting scan in frame %d\n",frame);
    }
    
    lm = (*hsi->scan_query)(hsi->data,trans,&para);
    if( p->debug ) {
      info("Retrieved %d hits in frame %d\n",lm->len,frame);
    }

    if( lm->len > 0 ) {
      add_GeneWiseHSPmanager_HSPset(gwh,lm->set[0],frame);
    }

    free_LinearHSPmanager(lm);
  }

  /* if this is empty, do something! */

  if( gwh->len == 0 ) {
    info("For genewise %s vs %s, no HSPs generated, heuristic failed",prot->name,dna->name);
    return NULL;
  }


 
  /* sort by score, descend list */

  sort_GeneWiseHSPmanager(gwh,compare_GeneWiseHSP_score);
  
  if( p->debug == TRUE ) {
    for(i=0;i<gwh->len;i++) {
      fprintf(stdout,"GWHSP [Before] %d,%d to %d,%d %d\n",
	      gwh->hsp[i]->query_start,
	      gwh->hsp[i]->query_end,
	      gwh->hsp[i]->target_start,
	      gwh->hsp[i]->target_end,
	      gwh->hsp[i]->score);
    }
  }
  
  final = GeneWiseHSPmanager_alloc_std();

  for(i=0;i<gwh->len;i++) {
    can_use = 1;
    for(j=0;j<final->len;j++) {
      if( consistent_GeneWiseHSP(final->hsp[j],gwh->hsp[i]) == 0 ) {
	can_use = 0;
	break;
      }
    }

    if( can_use == 1 ) {
      if( p->debug == TRUE ) {
	fprintf(stdout,"GWHSP accepting %d,%d to %d,%d\n",
	      gwh->hsp[i]->query_start,
	      gwh->hsp[i]->query_end,
	      gwh->hsp[i]->target_start,
	      gwh->hsp[i]->target_end
		);
      }
      add_GeneWiseHSPmanager(final,hard_link_GeneWiseHSP(gwh->hsp[i]));
    } else {
      if( p->debug == TRUE ) {
	fprintf(stdout,"GWHSP rejecting %d,%d to %d,%d\n",
	      gwh->hsp[i]->query_start,
	      gwh->hsp[i]->query_end,
	      gwh->hsp[i]->target_start,
	      gwh->hsp[i]->target_end
		);
      }
    }
  }


  out = DPEnvelope_alloc_std();

  sort_GeneWiseHSPmanager(final,compare_GeneWiseHSP_start);
  
  prev_query = final->hsp[0]->query_start - p->edge_query;
  prev_target = final->hsp[0]->target_start - p->edge_target;

  if( prev_query < 0 ) {
    prev_query = 0;
  }
  
  if( prev_target < 0 ) {
    prev_target = 0;
  }

  for(i=0;i<final->len;i++) {
    /* rectangle from previous exon to this one */
    dpunit = DPUnit_alloc();
    dpunit->type = DPENV_RECT;
    dpunit->starti = prev_query;
    dpunit->startj = prev_target;
    dpunit->height = final->hsp[i]->query_start - dpunit->starti + p->splice_spread;
    dpunit->length = final->hsp[i]->target_start - dpunit->startj + p->splice_spread*3;
    
    if( dpunit->height <= 0 || dpunit->length <= 0 ) {
      if( p->debug ) {
	fprintf(stdout,"On position %d , error with the jigging %d plays %d vs %d plays %d",i,
	      final->hsp[i]->query_start,dpunit->starti,final->hsp[i]->target_start,dpunit->startj);
      }
      if( dpunit->height <= 0 ) {
	dpunit->height = 1;
      } 
      if( dpunit->length <= 0 ) {
	dpunit->length = 1;
      }
      
    }

    if( p->debug ) {
      fprintf(stdout,"GWHSP Bridging to %d,%d to %d,%d from %dth element with %d,%d coordinate\n",
	      dpunit->starti,
	      dpunit->startj,
	      dpunit->height,
	      dpunit->length,i,final->hsp[i]->query_start,
	      final->hsp[i]->target_start
	      );
    }

    add_DPEnvelope(out,dpunit);

    /* rectange for this exon */

    dpunit = DPUnit_alloc();
    dpunit->type = DPENV_RECT;

    dpunit->starti = final->hsp[i]->query_start  - p->splice_spread;
    dpunit->startj = final->hsp[i]->target_start - p->splice_spread*3;;
    dpunit->height = final->hsp[i]->query_end - dpunit->starti + p->splice_spread;
    dpunit->length = final->hsp[i]->target_end - dpunit->startj + p->splice_spread*3;

    if( dpunit->height <= 0 || dpunit->length <= 0 ) {
      fatal("On exon position %d, error with the jigging of the positions",i);
    }

    if( p->debug ) {
      fprintf(stdout,"GWHSP Exon: %d,%d to %d,%d from %dth element with %d,%d coordinate\n",
	      dpunit->starti,
	      dpunit->startj,
	      dpunit->height,
	      dpunit->length,i,final->hsp[i]->query_start,
	      final->hsp[i]->target_start
	      );
    }
   
    add_DPEnvelope(out,dpunit);
    
    prev_query = final->hsp[i]->query_end - p->splice_spread;
    prev_target = final->hsp[i]->target_end - p->splice_spread*3;
  }

  free_GeneWiseHSPmanager(gwh);
  free_GeneWiseHSPmanager(final);
  free_HSPScanInterface(hsi);


  if( p->debug ) {
    fprintf(stdout,"GWHSP - exited with %d DP units\n",out->len);
    fflush(stdout);
  }

  

  return out;
}



# line 291 "genewisehsp.dy"
int consistent_GeneWiseHSP(GeneWiseHSP * true,GeneWiseHSP * proposed)
{
  int query_centre;
  int target_centre;

  assert(true);
  assert(proposed);
  /* is this left or right of the true HSP */
  query_centre = (proposed->query_start + proposed->query_end) / 2;
  target_centre = (proposed->target_start + proposed->target_end) / 2;
  
  /* overlap criteria first */
  if( query_centre >= true->query_start && query_centre <= true->query_end ) {
    return 0;
  }
  
  /* left or right */
  if( query_centre > true->query_end ) {
    if( target_centre < true->target_end ) {
      /* no - inconsistent */
      return 0;
    }
  } else {
    if( target_centre > true->target_start ) {
      return 0;
    }
  }


  /* overlap */

  if( target_centre >= true->target_start && target_centre <= true->target_end ) {
    return 0;
  }


  return 1;

}

# line 331 "genewisehsp.dy"
int compare_GeneWiseHSP_start(GeneWiseHSP * one,GeneWiseHSP * two)
{
  return one->query_start - two->query_start;
}

/* Function:  compare_GeneWiseHSP_score(one,two)
 *
 * Descrip:    internal for sort by score
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [GeneWiseHSP *]
 * Arg:        two [UNKN ] Undocumented argument [GeneWiseHSP *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 340 "genewisehsp.dy"
int compare_GeneWiseHSP_score(GeneWiseHSP * one,GeneWiseHSP * two)
{
  return two->score - one->score;
}


# line 346 "genewisehsp.dy"
void add_GeneWiseHSPmanager_HSPset(GeneWiseHSPmanager * gwh,HSPset * set,int frame)
{
  int i;
  GeneWiseHSP * h;


  for(i=0;i<set->len;i++) {
    h = GeneWiseHSP_alloc();
    h->query_start = set->hsp[i]->target_start;
    h->query_end   = set->hsp[i]->target_start+set->hsp[i]->length;
    h->target_start = set->hsp[i]->query_start*3 + frame;
    h->target_end  = (set->hsp[i]->query_start+set->hsp[i]->length)*3 + frame;
    h->score = set->hsp[i]->score;
    add_GeneWiseHSPmanager(gwh,h);
  }
}

# line 344 "genewisehsp.c"
/* Function:  hard_link_GeneWiseRunPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseRunPara *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseRunPara *]
 *
 */
GeneWiseRunPara * hard_link_GeneWiseRunPara(GeneWiseRunPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneWiseRunPara object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneWiseRunPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseRunPara *]
 *
 */
GeneWiseRunPara * GeneWiseRunPara_alloc(void) 
{
    GeneWiseRunPara * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWiseRunPara *) ckalloc (sizeof(GeneWiseRunPara))) == NULL)  {  
      warn("GeneWiseRunPara_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->use_hsp = TRUE; 
    out->edge_query = 10;    
    out->edge_target = 3000; 
    out->splice_spread = 4;  
    out->debug = FALSE;  


    return out;  
}    


/* Function:  free_GeneWiseRunPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseRunPara *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseRunPara *]
 *
 */
GeneWiseRunPara * free_GeneWiseRunPara(GeneWiseRunPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneWiseRunPara obj. Should be trappable");   
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


/* Function:  hard_link_GeneWiseHSP(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseHSP *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSP *]
 *
 */
GeneWiseHSP * hard_link_GeneWiseHSP(GeneWiseHSP * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneWiseHSP object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneWiseHSP_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSP *]
 *
 */
GeneWiseHSP * GeneWiseHSP_alloc(void) 
{
    GeneWiseHSP * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWiseHSP *) ckalloc (sizeof(GeneWiseHSP))) == NULL)  {  
      warn("GeneWiseHSP_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->query_start = 0;    
    out->query_end = 0;  
    out->target_start = 0;   
    out->target_end = 0; 
    out->score = 0;  
    out->frame = 0;  


    return out;  
}    


/* Function:  free_GeneWiseHSP(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseHSP *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSP *]
 *
 */
GeneWiseHSP * free_GeneWiseHSP(GeneWiseHSP * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneWiseHSP obj. Should be trappable");   
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


/* Function:  swap_GeneWiseHSPmanager(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GeneWiseHSPmanager
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GeneWiseHSP **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GeneWiseHSPmanager(GeneWiseHSP ** list,int i,int j)  
{
    GeneWiseHSP * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GeneWiseHSPmanager(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GeneWiseHSPmanager which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GeneWiseHSP **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GeneWiseHSPmanager(GeneWiseHSP ** list,int left,int right,int (*comp)(GeneWiseHSP * ,GeneWiseHSP * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GeneWiseHSPmanager(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GeneWiseHSPmanager (list,++last,i); 
      }  
    swap_GeneWiseHSPmanager (list,left,last);    
    qsort_GeneWiseHSPmanager(list,left,last-1,comp); 
    qsort_GeneWiseHSPmanager(list,last+1,right,comp);    
}    


/* Function:  sort_GeneWiseHSPmanager(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GeneWiseHSPmanager
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GeneWiseHSPmanager *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GeneWiseHSPmanager(GeneWiseHSPmanager * obj,int (*comp)(GeneWiseHSP *, GeneWiseHSP *)) 
{
    qsort_GeneWiseHSPmanager(obj->hsp,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_GeneWiseHSPmanager(obj,len)
 *
 * Descrip:    Really an internal function for add_GeneWiseHSPmanager
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneWiseHSPmanager *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GeneWiseHSPmanager(GeneWiseHSPmanager * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GeneWiseHSPmanager called with no need"); 
      return TRUE;   
      }  


    if( (obj->hsp = (GeneWiseHSP ** ) ckrealloc (obj->hsp,sizeof(GeneWiseHSP *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GeneWiseHSPmanager, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GeneWiseHSPmanager(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GeneWiseHSPmanager *]
 * Arg:        add [OWNER] Object to add to the list [GeneWiseHSP *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GeneWiseHSPmanager(GeneWiseHSPmanager * obj,GeneWiseHSP * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GeneWiseHSPmanager(obj,obj->len + GeneWiseHSPmanagerLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->hsp[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_GeneWiseHSPmanager(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GeneWiseHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GeneWiseHSPmanager(GeneWiseHSPmanager * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->hsp[i] != NULL)   {  
        free_GeneWiseHSP(obj->hsp[i]);   
        obj->hsp[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GeneWiseHSPmanager_alloc_std(void)
 *
 * Descrip:    Equivalent to GeneWiseHSPmanager_alloc_len(GeneWiseHSPmanagerLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSPmanager *]
 *
 */
GeneWiseHSPmanager * GeneWiseHSPmanager_alloc_std(void) 
{
    return GeneWiseHSPmanager_alloc_len(GeneWiseHSPmanagerLISTLENGTH);   
}    


/* Function:  GeneWiseHSPmanager_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSPmanager *]
 *
 */
GeneWiseHSPmanager * GeneWiseHSPmanager_alloc_len(int len) 
{
    GeneWiseHSPmanager * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GeneWiseHSPmanager_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->hsp = (GeneWiseHSP ** ) ckcalloc (len,sizeof(GeneWiseHSP *))) == NULL)  {  
      warn("Warning, ckcalloc failed in GeneWiseHSPmanager_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GeneWiseHSPmanager(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneWiseHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSPmanager *]
 *
 */
GeneWiseHSPmanager * hard_link_GeneWiseHSPmanager(GeneWiseHSPmanager * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneWiseHSPmanager object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneWiseHSPmanager_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSPmanager *]
 *
 */
GeneWiseHSPmanager * GeneWiseHSPmanager_alloc(void) 
{
    GeneWiseHSPmanager * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneWiseHSPmanager *) ckalloc (sizeof(GeneWiseHSPmanager))) == NULL)    {  
      warn("GeneWiseHSPmanager_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->hsp = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_GeneWiseHSPmanager(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneWiseHSPmanager *]
 *
 * Return [UNKN ]  Undocumented return value [GeneWiseHSPmanager *]
 *
 */
GeneWiseHSPmanager * free_GeneWiseHSPmanager(GeneWiseHSPmanager * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneWiseHSPmanager obj. Should be trappable");    
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
    if( obj->hsp != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->hsp[i] != NULL) 
          free_GeneWiseHSP(obj->hsp[i]); 
        }  
      ckfree(obj->hsp);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

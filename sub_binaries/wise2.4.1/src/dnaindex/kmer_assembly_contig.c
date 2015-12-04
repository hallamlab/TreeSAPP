#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_assembly_contig.h"

# line 34 "kmer_assembly_contig.dy"
KmerAssemblyContig * find_strict_mirrored_KmerAssemblyContig(KmerAssemblyContigSet * kcs,KmerAssemblyContig * c)
{
  int i;
  kmer_t start_rev;
  kmer_t end_rev;

  assert(kcs != NULL);
  assert(kcs->kai != NULL);
  assert(c != NULL);

  start_rev = reverse_complement_dna_number(c->start->number,kcs->kai->kii->kmer_size);
  end_rev   = reverse_complement_dna_number(c->end->number,kcs->kai->kii->kmer_size);

  for(i=0;i < kcs->len;i++) {
    if( kcs->contig[i]->mirror != NULL || kcs->contig[i]->is_mirror == 1) {
      continue;
    }
    /*    fprintf(stderr,"Looking at %ld,%ld vs %ld,%ld\n",kcs->contig[i]->start->number,start_rev,kcs->contig[i]->end->number,end_rev);*/

    if( kcs->contig[i]->start->number == end_rev && kcs->contig[i]->end->number == start_rev ) {
      /* should check labels */
      warn("Note to Ewan: Should check labels");
      return kcs->contig[i];
    }
  }

  return NULL;
}

# line 63 "kmer_assembly_contig.dy"
KmerAssemblyContigSet * KmerAssemblyContigSet_from_KmerAssemblyIndex(KmerAssemblyIndex * kai)
{
  KmerAssemblyContigSet * out;
  KmerAssemblyContigSet * final;
  kmer_t kmer;
  KmerAssemblyNode * node;
  KmerAssemblyContig * contig;
  boolean is_left_end;
  int i;

  assert(kai != NULL);
  assert(kai->kii != NULL );
  assert(kai->kii->next_filled_kmer != NULL);


  out = KmerAssemblyContigSet_alloc_std();
  out->kai = kai;

  kmer = -1;
  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer);

  fprintf(stderr,"CONTIG BUILD\n");

  for(;kmer != -1;  kmer = (*kai->kii->next_filled_kmer)(kai->kii->handle,kmer)) {
    node = (*kai->kii->retrieve_by_kmer)(kai->kii->handle,kmer);
    for(;node != NULL; node = node->node_chain ) {
      if( node->next_len == 0 ) {
	continue;
      }


      /* test to see if it is start; no back link or all TANGLED links */
      is_left_end = 0;
      if( node->prev_len == 0 ) {
	is_left_end = 1;
	/* this is the other main case, test before we do loop */
      } else if( node->prev_len == 1 && !(node->prev[0]->state & KMER_ASSEMBLY_PREV_TANGLED) ) {
	is_left_end = 0;
      } else {
	for(i=0;i<node->prev_len;i++) {
	  fprintf(stderr,"Node position %d, state %d\n",i,node->prev[i]->state);
	  if( !(node->prev[i]->state & KMER_ASSEMBLY_NEXT_TANGLED) ) {
	    is_left_end = 0;
	    break;
	  }
	}
	if( i >= node->prev_len ) {
	  is_left_end = 1;
	}
      }


      if( is_left_end == 0 && node->prev_len > 1 ) {
	fprintf(stderr,"Very weird: multi-prev but not left end, number is %d, state of 0 %d\n",node->prev_len,node->prev[0]->state);
      }
	     
      if( is_left_end &&  node->next_len == 1) {
	contig = new_KmerAssemblyContig(node);
	add_KmerAssemblyContigSet(out,contig);
      }
    }
  }

  return out;
}


# line 130 "kmer_assembly_contig.dy"
KmerAssemblyContig * new_KmerAssemblyContig(KmerAssemblyNode * node)
{
  KmerAssemblyContig * out;
  KmerAssemblyLink * link;
  int len;

  assert(node != NULL);

  if( node->next_len != 1 ) {
    warn("Cannot build a KmerAssemblyContig from a node with %d out going links",node->next_len);
  }

  out = KmerAssemblyContig_alloc();
  out->start = node;
  out->max_depth = 0;

  if( node->prev_len > 1 ) {
    out->clean_start = 0;
  } else {
    out->clean_start = 1;
  }

  len  = 0;
  for(link = node->next[0];link != NULL && link->next->next_len == 1; link = link->next->next[0],len++) {
    if( link->state & KMER_ASSEMBLY_PREV_TANGLED ) {
      break;
    }

    if( link->prev->prev_len > 1 ) {
      fprintf(stderr,"Assembly is still tangled!");
      break;
      /* bad end found */
    }

    if( link->sequence_label_len > out->max_depth ) {
      out->max_depth = link->sequence_label_len;
    }
  }

  out->end = link->next;

  if( link->next->next_len != 0 ) {
    out->clean_end = 0;
  } else {
    out->clean_end = 1;
  }

  out->len = len;
  
  return out;
}


# line 183 "kmer_assembly_contig.dy"
Assembly * Assembly_from_KmerAssemblyIndex(KmerAssemblyIndex * kai,KmerAssemblyContigPara * p)
{
  Assembly * out;
  KmerAssemblyContigSet * kacs;
  KmerAssemblyContigSet * final;

  KmerAssemblyContig * mirror;
  int i;

  kacs = KmerAssemblyContigSet_from_KmerAssemblyIndex(kai);

  final = KmerAssemblyContigSet_alloc_std();

  for(i=0;i<kacs->len;i++) {
    if( kacs->contig[i]->is_mirror == 1 ) {
      continue;
    }

    if( kacs->contig[i]->mirror == NULL ) {
      if( (mirror = find_strict_mirrored_KmerAssemblyContig(kacs,kacs->contig[i])) != NULL ) {
	kacs->contig[i]->mirror = hard_link_KmerAssemblyContig(mirror);
	mirror->is_mirror = 1;
	add_KmerAssemblyContigSet(final,hard_link_KmerAssemblyContig(kacs->contig[i]));
      } else {
	warn("Unable to mirror contig, adding anyway");
	add_KmerAssemblyContigSet(final,hard_link_KmerAssemblyContig(kacs->contig[i]));
      }
    }

  }
  out = Assembly_from_KmerAssemblyContigSet(final,p);

  free_KmerAssemblyContigSet(kacs);
  free_KmerAssemblyContigSet(final);

  return out;
}

# line 221 "kmer_assembly_contig.dy"
Assembly * Assembly_from_KmerAssemblyContigSet(KmerAssemblyContigSet * kacs,KmerAssemblyContigPara * p)
{
  Assembly * out;
  AssemblyContig * ac;
  int i;

  assert(kacs != NULL);
  assert(p != NULL);

  out = Assembly_alloc_std();

  for(i=0;i<kacs->len;i++) {
    if( kacs->contig[i]->len > p->minimum_len && kacs->contig[i]->max_depth >= p->minimum_depth ) {
      ac = AssemblyContig_from_KmerAssemblyContig(kacs->contig[i]);
      add_Assembly(out,ac);
    }
  }

  return out;
}

# line 242 "kmer_assembly_contig.dy"
AssemblyContig * AssemblyContig_from_KmerAssemblyContig(KmerAssemblyContig * kac)
{
  AssemblyContig * out;
  KmerAssemblyLink * link;
  Sequence * con;
  int i;
  char buffer[512];


  con = Sequence_alloc();
  con->seq = calloc(kac->len+1,sizeof(char));
  con->len = kac->len;
  con->maxlen = con->len;
  sprintf(buffer,"contig_%ld",(long int)kac);

  con->name = stringalloc(buffer);

  out = AssemblyContig_alloc_std();

  fprintf(stderr,"Starting contig fetching...\n");

  for(i=0,link = kac->start->next[0] ;link != NULL && link->next->next_len == 1;link = link->next->next[0],i++) {
    con->seq[i] = link->base;
    /* NOT DEALING WITH READ PLACEMENT YET */
    if( link->next == kac->end ) {
      break;
    }
  }

  out->clean_start = kac->clean_start;
  out->clean_end   = kac->clean_end;
  out->max_depth   = kac->max_depth;
  

  con->seq[i] = '\0';
  con->len = strlen(con->seq);
  
  out->consensus = con;

  return out;
}



# line 263 "kmer_assembly_contig.c"
/* Function:  hard_link_KmerAssemblyContig(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerAssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContig *]
 *
 */
KmerAssemblyContig * hard_link_KmerAssemblyContig(KmerAssemblyContig * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a KmerAssemblyContig object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  KmerAssemblyContig_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContig *]
 *
 */
KmerAssemblyContig * KmerAssemblyContig_alloc(void) 
{
    KmerAssemblyContig * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(KmerAssemblyContig *) ckalloc (sizeof(KmerAssemblyContig))) == NULL)    {  
      warn("KmerAssemblyContig_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->len = 0;    
    out->clean_start = FALSE;    
    out->clean_end = FALSE;  
    out->max_depth = 0;  
    out->mirror = NULL;  
    out->is_mirror = 0;  


    return out;  
}    


/* Function:  free_KmerAssemblyContig(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerAssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContig *]
 *
 */
KmerAssemblyContig * free_KmerAssemblyContig(KmerAssemblyContig * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a KmerAssemblyContig obj. Should be trappable");    
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
    /* obj->start is linked in */ 
    /* obj->end is linked in */ 
    if( obj->mirror != NULL) 
      free_KmerAssemblyContig(obj->mirror);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_KmerAssemblyContigSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_KmerAssemblyContigSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [KmerAssemblyContig **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_KmerAssemblyContigSet(KmerAssemblyContig ** list,int i,int j)  
{
    KmerAssemblyContig * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_KmerAssemblyContigSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_KmerAssemblyContigSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [KmerAssemblyContig **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_KmerAssemblyContigSet(KmerAssemblyContig ** list,int left,int right,int (*comp)(KmerAssemblyContig * ,KmerAssemblyContig * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_KmerAssemblyContigSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_KmerAssemblyContigSet (list,++last,i);  
      }  
    swap_KmerAssemblyContigSet (list,left,last); 
    qsort_KmerAssemblyContigSet(list,left,last-1,comp);  
    qsort_KmerAssemblyContigSet(list,last+1,right,comp); 
}    


/* Function:  sort_KmerAssemblyContigSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_KmerAssemblyContigSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [KmerAssemblyContigSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_KmerAssemblyContigSet(KmerAssemblyContigSet * obj,int (*comp)(KmerAssemblyContig *, KmerAssemblyContig *)) 
{
    qsort_KmerAssemblyContigSet(obj->contig,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_KmerAssemblyContigSet(obj,len)
 *
 * Descrip:    Really an internal function for add_KmerAssemblyContigSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [KmerAssemblyContigSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_KmerAssemblyContigSet(KmerAssemblyContigSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_KmerAssemblyContigSet called with no need");  
      return TRUE;   
      }  


    if( (obj->contig = (KmerAssemblyContig ** ) ckrealloc (obj->contig,sizeof(KmerAssemblyContig *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_KmerAssemblyContigSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_KmerAssemblyContigSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [KmerAssemblyContigSet *]
 * Arg:        add [OWNER] Object to add to the list [KmerAssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_KmerAssemblyContigSet(KmerAssemblyContigSet * obj,KmerAssemblyContig * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_KmerAssemblyContigSet(obj,obj->len + KmerAssemblyContigSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->contig[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_KmerAssemblyContigSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [KmerAssemblyContigSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_KmerAssemblyContigSet(KmerAssemblyContigSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->contig[i] != NULL)    {  
        free_KmerAssemblyContig(obj->contig[i]); 
        obj->contig[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  KmerAssemblyContigSet_alloc_std(void)
 *
 * Descrip:    Equivalent to KmerAssemblyContigSet_alloc_len(KmerAssemblyContigSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigSet *]
 *
 */
KmerAssemblyContigSet * KmerAssemblyContigSet_alloc_std(void) 
{
    return KmerAssemblyContigSet_alloc_len(KmerAssemblyContigSetLISTLENGTH); 
}    


/* Function:  KmerAssemblyContigSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigSet *]
 *
 */
KmerAssemblyContigSet * KmerAssemblyContigSet_alloc_len(int len) 
{
    KmerAssemblyContigSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = KmerAssemblyContigSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->contig = (KmerAssemblyContig ** ) ckcalloc (len,sizeof(KmerAssemblyContig *))) == NULL) {  
      warn("Warning, ckcalloc failed in KmerAssemblyContigSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_KmerAssemblyContigSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerAssemblyContigSet *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigSet *]
 *
 */
KmerAssemblyContigSet * hard_link_KmerAssemblyContigSet(KmerAssemblyContigSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a KmerAssemblyContigSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  KmerAssemblyContigSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigSet *]
 *
 */
KmerAssemblyContigSet * KmerAssemblyContigSet_alloc(void) 
{
    KmerAssemblyContigSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(KmerAssemblyContigSet *) ckalloc (sizeof(KmerAssemblyContigSet))) == NULL)  {  
      warn("KmerAssemblyContigSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->contig = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_KmerAssemblyContigSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerAssemblyContigSet *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigSet *]
 *
 */
KmerAssemblyContigSet * free_KmerAssemblyContigSet(KmerAssemblyContigSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a KmerAssemblyContigSet obj. Should be trappable"); 
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
    if( obj->contig != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->contig[i] != NULL)  
          free_KmerAssemblyContig(obj->contig[i]);   
        }  
      ckfree(obj->contig);   
      }  
    /* obj->kai is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_KmerAssemblyContigPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [KmerAssemblyContigPara *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigPara *]
 *
 */
KmerAssemblyContigPara * hard_link_KmerAssemblyContigPara(KmerAssemblyContigPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a KmerAssemblyContigPara object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  KmerAssemblyContigPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigPara *]
 *
 */
KmerAssemblyContigPara * KmerAssemblyContigPara_alloc(void) 
{
    KmerAssemblyContigPara * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(KmerAssemblyContigPara *) ckalloc (sizeof(KmerAssemblyContigPara))) == NULL)    {  
      warn("KmerAssemblyContigPara_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->minimum_len = 50;   
    out->minimum_depth = 1;  


    return out;  
}    


/* Function:  free_KmerAssemblyContigPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [KmerAssemblyContigPara *]
 *
 * Return [UNKN ]  Undocumented return value [KmerAssemblyContigPara *]
 *
 */
KmerAssemblyContigPara * free_KmerAssemblyContigPara(KmerAssemblyContigPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a KmerAssemblyContigPara obj. Should be trappable");    
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

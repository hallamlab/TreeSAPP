#ifdef _cplusplus
extern "C" {
#endif
#include "assembly.h"



# line 80 "assembly.dy"
void show_help_AssemblyOutputPara(FILE * ofp)
{
  fprintf(ofp,"Assembly Output options\n");
  fprintf(ofp,"  -as_min_size [500] min size of contig to show\n");

}

# line 87 "assembly.dy"
AssemblyOutputPara * new_AssemblyOutputPara_from_argv(int * argc,char ** argv)
{
  AssemblyOutputPara * out;

  out = AssemblyOutputPara_alloc();

  out->min_size = 500;

  strip_out_integer_argument(argc,argv,"as_min_size",&out->min_size);

  return out;
}

# line 100 "assembly.dy"
AssemblyOpaqueTypeSet * homopolymer_AssemblyOpaqueTypeSet(void)
{
  AssemblyOpaqueTypeSet * out;
  AssemblyOpaqueType * t;
  char * base = "ATGC";
  int i;

  out = AssemblyOpaqueTypeSet_alloc_std();

  for(i=0;i<4;i++) {
    t = new_homopolymer_AssemblyOpaqueType(i,base[i]);
    add_AssemblyOpaqueTypeSet(out,t);
  }

  return out;
}


# line 118 "assembly.dy"
int annotate_AssemblyOpaqueFeatures(AssemblyOpaqueTypeSet * aots,AssemblySequence * aseq,int kmer_size)
{
  int count = 0;
  int i,j;
  int k;
  AssemblyOpaqueFeature * aof;
  AssemblyOpaqueFeature * prev;

  for(i=0;i<aots->len;i++) {
    prev = NULL;
    for(j=0;j<aseq->len;) {
      if( aseq->seq->seq[j] == aots->type[i]->base ) {
	for(k=0;k<kmer_size;k++) {
	  if( aseq->seq->seq[j+k] != aots->type[i]->base ) {
	    break;
	  }
	}

	if( k < kmer_size ) {
	  j += k;
	  continue;
	}

	/* else, is a good position, extend to the end of this run */
	for(;k+j<aseq->len;k++) {
	  if( aseq->seq->seq[j+k] != aots->type[i]->base ) {
	    break;
	  }
	}

	/* if the start is within kmer of prev, then simply extend prev */

	if( prev != NULL && prev->start+prev->length+kmer_size >= j ) {
	  prev->length = j+k - prev->start;
	} else {
	  /* new feature */
	  aof = AssemblyOpaqueFeature_alloc();
	  aof->start = j;
	  aof->length = k;
	  aof->type = aots->type[i];
	  add_opq_AssemblySequence(aseq,aof);
	  prev = aof;
	  count++;
	}
      } else {
	j++;
      }
    }
  }

  return count;
}


# line 172 "assembly.dy"
void show_AssemblySequence(AssemblySequence * aseq,FILE * ofp)
{
  int i;
  assert(aseq!=NULL);

  for(i=0;i<aseq->opq_len;i++) {
    fprintf(ofp,"Opaque type %c from %d to %d\n",aseq->opaque[i]->type->base,aseq->opaque[i]->start,aseq->opaque[i]->start+aseq->opaque[i]->length);
  }
  write_fasta_Sequence(aseq->seq,ofp);
  fprintf(ofp,"//\n");

}

# line 185 "assembly.dy"
AssemblyOpaqueType * new_homopolymer_AssemblyOpaqueType(int int_type,char base)
{
  AssemblyOpaqueType * out;

  out = AssemblyOpaqueType_alloc();
  out->base = base;
  out->int_type = int_type;

  return out;
}

# line 196 "assembly.dy"
AssemblySequence * read_plain_fasta_AssemblySequence(FILE * ifp,int report_log,FILE * report)
{
  AssemblySequence * out;
  Sequence * seq;

  seq = read_large_dna_Sequence(ifp,report_log,report);

  if( seq == NULL ) {
    return NULL;
  }

  out = AssemblySequence_alloc_std();

  out->seq = seq;

  return out;
}

# line 214 "assembly.dy"
AssemblySequence * mirrored_AssemblySequence(AssemblySequence * aseq)
{
  AssemblySequence * out;

  out = AssemblySequence_alloc();
  out->seq = reverse_complement_Sequence(aseq->seq);
  out->mirror_seq = 1;

  out->mirror   = aseq;
  aseq->mirror  = out;
  
  return out;
}

# line 228 "assembly.dy"
void dump_contigs_as_fasta_Assembly(Assembly * assembly,AssemblyOutputPara * aop,FILE * ofp)
{
  int i;

  for(i=0;i<assembly->len;i++) {
    AssemblyContig * c = assembly->contig[i];
    if( aop->min_size > c->consensus->len ) {
      continue;
    }

    fprintf(ofp,">%s max_depth=%d clean_start=%d clean_end=%d\n",c->consensus->name,c->max_depth,c->clean_start,c->clean_end);
    show_line(c->consensus->seq,60,ofp);
  }

}


# line 181 "assembly.c"
/* Function:  hard_link_AssemblySequenceEdit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblySequenceEdit *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceEdit *]
 *
 */
AssemblySequenceEdit * hard_link_AssemblySequenceEdit(AssemblySequenceEdit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AssemblySequenceEdit object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AssemblySequenceEdit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceEdit *]
 *
 */
AssemblySequenceEdit * AssemblySequenceEdit_alloc(void) 
{
    AssemblySequenceEdit * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AssemblySequenceEdit *) ckalloc (sizeof(AssemblySequenceEdit))) == NULL)    {  
      warn("AssemblySequenceEdit_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 'u'; 
    out->base = 'u'; 


    return out;  
}    


/* Function:  free_AssemblySequenceEdit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblySequenceEdit *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequenceEdit *]
 *
 */
AssemblySequenceEdit * free_AssemblySequenceEdit(AssemblySequenceEdit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AssemblySequenceEdit obj. Should be trappable");  
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


/* Function:  hard_link_AssemblyOpaqueType(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyOpaqueType *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueType *]
 *
 */
AssemblyOpaqueType * hard_link_AssemblyOpaqueType(AssemblyOpaqueType * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AssemblyOpaqueType object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AssemblyOpaqueType_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueType *]
 *
 */
AssemblyOpaqueType * AssemblyOpaqueType_alloc(void) 
{
    AssemblyOpaqueType * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AssemblyOpaqueType *) ckalloc (sizeof(AssemblyOpaqueType))) == NULL)    {  
      warn("AssemblyOpaqueType_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->int_type = 0;   
    out->base = 'u'; 


    return out;  
}    


/* Function:  free_AssemblyOpaqueType(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyOpaqueType *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueType *]
 *
 */
AssemblyOpaqueType * free_AssemblyOpaqueType(AssemblyOpaqueType * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AssemblyOpaqueType obj. Should be trappable");    
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


/* Function:  swap_AssemblyOpaqueTypeSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AssemblyOpaqueTypeSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AssemblyOpaqueType **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AssemblyOpaqueTypeSet(AssemblyOpaqueType ** list,int i,int j)  
{
    AssemblyOpaqueType * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AssemblyOpaqueTypeSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AssemblyOpaqueTypeSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AssemblyOpaqueType **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AssemblyOpaqueTypeSet(AssemblyOpaqueType ** list,int left,int right,int (*comp)(AssemblyOpaqueType * ,AssemblyOpaqueType * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AssemblyOpaqueTypeSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AssemblyOpaqueTypeSet (list,++last,i);  
      }  
    swap_AssemblyOpaqueTypeSet (list,left,last); 
    qsort_AssemblyOpaqueTypeSet(list,left,last-1,comp);  
    qsort_AssemblyOpaqueTypeSet(list,last+1,right,comp); 
}    


/* Function:  sort_AssemblyOpaqueTypeSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AssemblyOpaqueTypeSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AssemblyOpaqueTypeSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj,int (*comp)(AssemblyOpaqueType *, AssemblyOpaqueType *)) 
{
    qsort_AssemblyOpaqueTypeSet(obj->type,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_AssemblyOpaqueTypeSet(obj,len)
 *
 * Descrip:    Really an internal function for add_AssemblyOpaqueTypeSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblyOpaqueTypeSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AssemblyOpaqueTypeSet called with no need");  
      return TRUE;   
      }  


    if( (obj->type = (AssemblyOpaqueType ** ) ckrealloc (obj->type,sizeof(AssemblyOpaqueType *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_AssemblyOpaqueTypeSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_AssemblyOpaqueTypeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblyOpaqueTypeSet *]
 * Arg:        add [OWNER] Object to add to the list [AssemblyOpaqueType *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj,AssemblyOpaqueType * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AssemblyOpaqueTypeSet(obj,obj->len + AssemblyOpaqueTypeSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->type[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_AssemblyOpaqueTypeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AssemblyOpaqueTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->type[i] != NULL)  {  
        free_AssemblyOpaqueType(obj->type[i]);   
        obj->type[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  AssemblyOpaqueTypeSet_alloc_std(void)
 *
 * Descrip:    Equivalent to AssemblyOpaqueTypeSet_alloc_len(AssemblyOpaqueTypeSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueTypeSet *]
 *
 */
AssemblyOpaqueTypeSet * AssemblyOpaqueTypeSet_alloc_std(void) 
{
    return AssemblyOpaqueTypeSet_alloc_len(AssemblyOpaqueTypeSetLISTLENGTH); 
}    


/* Function:  AssemblyOpaqueTypeSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueTypeSet *]
 *
 */
AssemblyOpaqueTypeSet * AssemblyOpaqueTypeSet_alloc_len(int len) 
{
    AssemblyOpaqueTypeSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AssemblyOpaqueTypeSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->type = (AssemblyOpaqueType ** ) ckcalloc (len,sizeof(AssemblyOpaqueType *))) == NULL)   {  
      warn("Warning, ckcalloc failed in AssemblyOpaqueTypeSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_AssemblyOpaqueTypeSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyOpaqueTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueTypeSet *]
 *
 */
AssemblyOpaqueTypeSet * hard_link_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AssemblyOpaqueTypeSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AssemblyOpaqueTypeSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueTypeSet *]
 *
 */
AssemblyOpaqueTypeSet * AssemblyOpaqueTypeSet_alloc(void) 
{
    AssemblyOpaqueTypeSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AssemblyOpaqueTypeSet *) ckalloc (sizeof(AssemblyOpaqueTypeSet))) == NULL)  {  
      warn("AssemblyOpaqueTypeSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_AssemblyOpaqueTypeSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyOpaqueTypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueTypeSet *]
 *
 */
AssemblyOpaqueTypeSet * free_AssemblyOpaqueTypeSet(AssemblyOpaqueTypeSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AssemblyOpaqueTypeSet obj. Should be trappable"); 
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
    if( obj->type != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->type[i] != NULL)    
          free_AssemblyOpaqueType(obj->type[i]); 
        }  
      ckfree(obj->type); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_AssemblyOpaqueFeature(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyOpaqueFeature *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueFeature *]
 *
 */
AssemblyOpaqueFeature * hard_link_AssemblyOpaqueFeature(AssemblyOpaqueFeature * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AssemblyOpaqueFeature object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AssemblyOpaqueFeature_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueFeature *]
 *
 */
AssemblyOpaqueFeature * AssemblyOpaqueFeature_alloc(void) 
{
    AssemblyOpaqueFeature * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AssemblyOpaqueFeature *) ckalloc (sizeof(AssemblyOpaqueFeature))) == NULL)  {  
      warn("AssemblyOpaqueFeature_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->length = 0; 


    return out;  
}    


/* Function:  free_AssemblyOpaqueFeature(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyOpaqueFeature *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOpaqueFeature *]
 *
 */
AssemblyOpaqueFeature * free_AssemblyOpaqueFeature(AssemblyOpaqueFeature * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AssemblyOpaqueFeature obj. Should be trappable"); 
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
    /* obj->type is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_AssemblySequence(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AssemblySequence
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AssemblySequenceEdit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AssemblySequence(AssemblySequenceEdit ** list,int i,int j)  
{
    AssemblySequenceEdit * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AssemblySequence(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AssemblySequence which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AssemblySequenceEdit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AssemblySequence(AssemblySequenceEdit ** list,int left,int right,int (*comp)(AssemblySequenceEdit * ,AssemblySequenceEdit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AssemblySequence(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AssemblySequence (list,++last,i);   
      }  
    swap_AssemblySequence (list,left,last);  
    qsort_AssemblySequence(list,left,last-1,comp);   
    qsort_AssemblySequence(list,last+1,right,comp);  
}    


/* Function:  sort_AssemblySequence(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AssemblySequence
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AssemblySequence *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AssemblySequence(AssemblySequence * obj,int (*comp)(AssemblySequenceEdit *, AssemblySequenceEdit *)) 
{
    qsort_AssemblySequence(obj->edit,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_AssemblySequence(obj,len)
 *
 * Descrip:    Really an internal function for add_AssemblySequence
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblySequence *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AssemblySequence(AssemblySequence * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AssemblySequence called with no need");   
      return TRUE;   
      }  


    if( (obj->edit = (AssemblySequenceEdit ** ) ckrealloc (obj->edit,sizeof(AssemblySequenceEdit *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_AssemblySequence, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_AssemblySequence(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblySequence *]
 * Arg:        add [OWNER] Object to add to the list [AssemblySequenceEdit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_AssemblySequence(AssemblySequence * obj,AssemblySequenceEdit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AssemblySequence(obj,obj->len + AssemblySequenceLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->edit[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_AssemblySequence(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AssemblySequence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_AssemblySequence(AssemblySequence * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->edit[i] != NULL)  {  
        free_AssemblySequenceEdit(obj->edit[i]); 
        obj->edit[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_opq_AssemblySequence(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_opq_AssemblySequence
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AssemblyOpaqueFeature **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_opq_AssemblySequence(AssemblyOpaqueFeature ** list,int i,int j)  
{
    AssemblyOpaqueFeature * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_opq_AssemblySequence(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_opq_AssemblySequence which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AssemblyOpaqueFeature **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_opq_AssemblySequence(AssemblyOpaqueFeature ** list,int left,int right,int (*comp)(AssemblyOpaqueFeature * ,AssemblyOpaqueFeature * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_opq_AssemblySequence(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_opq_AssemblySequence (list,++last,i);   
      }  
    swap_opq_AssemblySequence (list,left,last);  
    qsort_opq_AssemblySequence(list,left,last-1,comp);   
    qsort_opq_AssemblySequence(list,last+1,right,comp);  
}    


/* Function:  sort_opq_AssemblySequence(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_opq_AssemblySequence
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AssemblySequence *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_opq_AssemblySequence(AssemblySequence * obj,int (*comp)(AssemblyOpaqueFeature *, AssemblyOpaqueFeature *)) 
{
    qsort_opq_AssemblySequence(obj->opaque,0,obj->opq_len-1,comp);   
    return;  
}    


/* Function:  expand_opq_AssemblySequence(obj,len)
 *
 * Descrip:    Really an internal function for add_opq_AssemblySequence
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblySequence *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_opq_AssemblySequence(AssemblySequence * obj,int len) 
{


    if( obj->opq_maxlen > obj->opq_len )     {  
      warn("expand_AssemblySequenceopq_ called with no need");   
      return TRUE;   
      }  


    if( (obj->opaque = (AssemblyOpaqueFeature ** ) ckrealloc (obj->opaque,sizeof(AssemblyOpaqueFeature *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_AssemblySequence, returning FALSE"); 
      return FALSE;  
      }  
    obj->opq_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_opq_AssemblySequence(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblySequence *]
 * Arg:        add [OWNER] Object to add to the list [AssemblyOpaqueFeature *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_opq_AssemblySequence(AssemblySequence * obj,AssemblyOpaqueFeature * add) 
{
    if( obj->opq_len >= obj->opq_maxlen) {  
      if( expand_opq_AssemblySequence(obj,obj->opq_len + AssemblySequenceLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->opaque[obj->opq_len++]=add; 
    return TRUE; 
}    


/* Function:  flush_opq_AssemblySequence(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AssemblySequence *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_opq_AssemblySequence(AssemblySequence * obj) 
{
    int i;   


    for(i=0;i<obj->opq_len;i++)  { /*for i over list length*/ 
      if( obj->opaque[i] != NULL)    {  
        free_AssemblyOpaqueFeature(obj->opaque[i]);  
        obj->opaque[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->opq_len = 0;    
    return i;    
}    


/* Function:  AssemblySequence_alloc_std(void)
 *
 * Descrip:    Equivalent to AssemblySequence_alloc_len(AssemblySequenceLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequence *]
 *
 */
AssemblySequence * AssemblySequence_alloc_std(void) 
{
    return AssemblySequence_alloc_len(AssemblySequenceLISTLENGTH);   
}    


/* Function:  AssemblySequence_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequence *]
 *
 */
AssemblySequence * AssemblySequence_alloc_len(int len) 
{
    AssemblySequence * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AssemblySequence_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->edit = (AssemblySequenceEdit ** ) ckcalloc (len,sizeof(AssemblySequenceEdit *))) == NULL)   {  
      warn("Warning, ckcalloc failed in AssemblySequence_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->opaque = (AssemblyOpaqueFeature ** ) ckcalloc (len,sizeof(AssemblyOpaqueFeature *))) == NULL)   {  
      warn("Warning, ckcalloc failed in AssemblySequence_alloc_len");    
      return NULL;   
      }  
    out->opq_len = 0;    
    out->opq_maxlen = len;   


    return out;  
}    


/* Function:  hard_link_AssemblySequence(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblySequence *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequence *]
 *
 */
AssemblySequence * hard_link_AssemblySequence(AssemblySequence * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AssemblySequence object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AssemblySequence_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequence *]
 *
 */
AssemblySequence * AssemblySequence_alloc(void) 
{
    AssemblySequence * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AssemblySequence *) ckalloc (sizeof(AssemblySequence))) == NULL)    {  
      warn("AssemblySequence_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seq = NULL; 
    out->quality = NULL; 
    out->edit = NULL;    
    out->len = out->maxlen = 0;  
    out->state = 'u';    
    out->mirror_seq = 0; 
    out->opaque = NULL;  
    out->opq_len = out->opq_maxlen = 0;  
    out->orig = NULL;    
    out->abrev_repeat = NULL;    


    return out;  
}    


/* Function:  free_AssemblySequence(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblySequence *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequence *]
 *
 */
AssemblySequence * free_AssemblySequence(AssemblySequence * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AssemblySequence obj. Should be trappable");  
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
    if( obj->seq != NULL)    
      free_Sequence(obj->seq);   
    if( obj->quality != NULL)    
      ckfree(obj->quality);  
    if( obj->edit != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->edit[i] != NULL)    
          free_AssemblySequenceEdit(obj->edit[i]);   
        }  
      ckfree(obj->edit); 
      }  
    /* obj->mirror is linked in */ 
    /* obj->pair is linked in */ 
    if( obj->opaque != NULL) {  
      for(i=0;i<obj->opq_len;i++)    {  
        if( obj->opaque[i] != NULL)  
          free_AssemblyOpaqueFeature(obj->opaque[i]);    
        }  
      ckfree(obj->opaque);   
      }  
    if( obj->orig != NULL)   
      free_Sequence(obj->orig);  
    if( obj->abrev_repeat != NULL)   
      ckfree(obj->abrev_repeat);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_AssemblySequencePlacement(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblySequencePlacement *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequencePlacement *]
 *
 */
AssemblySequencePlacement * hard_link_AssemblySequencePlacement(AssemblySequencePlacement * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AssemblySequencePlacement object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AssemblySequencePlacement_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequencePlacement *]
 *
 */
AssemblySequencePlacement * AssemblySequencePlacement_alloc(void) 
{
    AssemblySequencePlacement * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AssemblySequencePlacement *) ckalloc (sizeof(AssemblySequencePlacement))) == NULL)  {  
      warn("AssemblySequencePlacement_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->aseq = NULL;    
    out->contig_start = 0;   
    out->contig_end = 0; 
    out->seq_start = 0;  
    out->seq_end = 0;    
    out->ungapped = FALSE;   


    return out;  
}    


/* Function:  free_AssemblySequencePlacement(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblySequencePlacement *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblySequencePlacement *]
 *
 */
AssemblySequencePlacement * free_AssemblySequencePlacement(AssemblySequencePlacement * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AssemblySequencePlacement obj. Should be trappable"); 
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
    if( obj->aseq != NULL)   
      free_AssemblySequence(obj->aseq);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_AssemblyContig(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AssemblyContig
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AssemblySequencePlacement **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AssemblyContig(AssemblySequencePlacement ** list,int i,int j)  
{
    AssemblySequencePlacement * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AssemblyContig(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AssemblyContig which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AssemblySequencePlacement **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AssemblyContig(AssemblySequencePlacement ** list,int left,int right,int (*comp)(AssemblySequencePlacement * ,AssemblySequencePlacement * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AssemblyContig(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AssemblyContig (list,++last,i); 
      }  
    swap_AssemblyContig (list,left,last);    
    qsort_AssemblyContig(list,left,last-1,comp); 
    qsort_AssemblyContig(list,last+1,right,comp);    
}    


/* Function:  sort_AssemblyContig(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AssemblyContig
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AssemblyContig *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AssemblyContig(AssemblyContig * obj,int (*comp)(AssemblySequencePlacement *, AssemblySequencePlacement *)) 
{
    qsort_AssemblyContig(obj->reads,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_AssemblyContig(obj,len)
 *
 * Descrip:    Really an internal function for add_AssemblyContig
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblyContig *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AssemblyContig(AssemblyContig * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AssemblyContig called with no need"); 
      return TRUE;   
      }  


    if( (obj->reads = (AssemblySequencePlacement ** ) ckrealloc (obj->reads,sizeof(AssemblySequencePlacement *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_AssemblyContig, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_AssemblyContig(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AssemblyContig *]
 * Arg:        add [OWNER] Object to add to the list [AssemblySequencePlacement *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_AssemblyContig(AssemblyContig * obj,AssemblySequencePlacement * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AssemblyContig(obj,obj->len + AssemblyContigLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->reads[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_AssemblyContig(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_AssemblyContig(AssemblyContig * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->reads[i] != NULL) {  
        free_AssemblySequencePlacement(obj->reads[i]);   
        obj->reads[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  AssemblyContig_alloc_std(void)
 *
 * Descrip:    Equivalent to AssemblyContig_alloc_len(AssemblyContigLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyContig *]
 *
 */
AssemblyContig * AssemblyContig_alloc_std(void) 
{
    return AssemblyContig_alloc_len(AssemblyContigLISTLENGTH);   
}    


/* Function:  AssemblyContig_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyContig *]
 *
 */
AssemblyContig * AssemblyContig_alloc_len(int len) 
{
    AssemblyContig * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AssemblyContig_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->reads = (AssemblySequencePlacement ** ) ckcalloc (len,sizeof(AssemblySequencePlacement *))) == NULL)    {  
      warn("Warning, ckcalloc failed in AssemblyContig_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_AssemblyContig(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyContig *]
 *
 */
AssemblyContig * hard_link_AssemblyContig(AssemblyContig * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AssemblyContig object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AssemblyContig_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyContig *]
 *
 */
AssemblyContig * AssemblyContig_alloc(void) 
{
    AssemblyContig * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AssemblyContig *) ckalloc (sizeof(AssemblyContig))) == NULL)    {  
      warn("AssemblyContig_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->consensus = NULL;   
    out->reads = NULL;   
    out->len = out->maxlen = 0;  
    out->clean_start = FALSE;    
    out->clean_end = FALSE;  
    out->max_depth = 0;  


    return out;  
}    


/* Function:  free_AssemblyContig(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyContig *]
 *
 */
AssemblyContig * free_AssemblyContig(AssemblyContig * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AssemblyContig obj. Should be trappable");    
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
    if( obj->consensus != NULL)  
      free_Sequence(obj->consensus);     
    if( obj->reads != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->reads[i] != NULL)   
          free_AssemblySequencePlacement(obj->reads[i]); 
        }  
      ckfree(obj->reads);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_Assembly(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Assembly
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AssemblyContig **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Assembly(AssemblyContig ** list,int i,int j)  
{
    AssemblyContig * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Assembly(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Assembly which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AssemblyContig **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Assembly(AssemblyContig ** list,int left,int right,int (*comp)(AssemblyContig * ,AssemblyContig * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Assembly(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Assembly (list,++last,i);   
      }  
    swap_Assembly (list,left,last);  
    qsort_Assembly(list,left,last-1,comp);   
    qsort_Assembly(list,last+1,right,comp);  
}    


/* Function:  sort_Assembly(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Assembly
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Assembly *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Assembly(Assembly * obj,int (*comp)(AssemblyContig *, AssemblyContig *)) 
{
    qsort_Assembly(obj->contig,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_Assembly(obj,len)
 *
 * Descrip:    Really an internal function for add_Assembly
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Assembly *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Assembly(Assembly * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Assembly called with no need");   
      return TRUE;   
      }  


    if( (obj->contig = (AssemblyContig ** ) ckrealloc (obj->contig,sizeof(AssemblyContig *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Assembly, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Assembly(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Assembly *]
 * Arg:        add [OWNER] Object to add to the list [AssemblyContig *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Assembly(Assembly * obj,AssemblyContig * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Assembly(obj,obj->len + AssemblyLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->contig[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_Assembly(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Assembly *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Assembly(Assembly * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->contig[i] != NULL)    {  
        free_AssemblyContig(obj->contig[i]); 
        obj->contig[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_chaff_Assembly(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_chaff_Assembly
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AssemblySequence **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_chaff_Assembly(AssemblySequence ** list,int i,int j)  
{
    AssemblySequence * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_chaff_Assembly(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_chaff_Assembly which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AssemblySequence **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_chaff_Assembly(AssemblySequence ** list,int left,int right,int (*comp)(AssemblySequence * ,AssemblySequence * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_chaff_Assembly(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_chaff_Assembly (list,++last,i); 
      }  
    swap_chaff_Assembly (list,left,last);    
    qsort_chaff_Assembly(list,left,last-1,comp); 
    qsort_chaff_Assembly(list,last+1,right,comp);    
}    


/* Function:  sort_chaff_Assembly(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_chaff_Assembly
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Assembly *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_chaff_Assembly(Assembly * obj,int (*comp)(AssemblySequence *, AssemblySequence *)) 
{
    qsort_chaff_Assembly(obj->chaff,0,obj->chaff_len-1,comp);    
    return;  
}    


/* Function:  expand_chaff_Assembly(obj,len)
 *
 * Descrip:    Really an internal function for add_chaff_Assembly
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Assembly *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_chaff_Assembly(Assembly * obj,int len) 
{


    if( obj->chaff_maxlen > obj->chaff_len )     {  
      warn("expand_Assemblychaff_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->chaff = (AssemblySequence ** ) ckrealloc (obj->chaff,sizeof(AssemblySequence *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Assembly, returning FALSE"); 
      return FALSE;  
      }  
    obj->chaff_maxlen = len; 
    return TRUE; 
}    


/* Function:  add_chaff_Assembly(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Assembly *]
 * Arg:        add [OWNER] Object to add to the list [AssemblySequence *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_chaff_Assembly(Assembly * obj,AssemblySequence * add) 
{
    if( obj->chaff_len >= obj->chaff_maxlen) {  
      if( expand_chaff_Assembly(obj,obj->chaff_len + AssemblyLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->chaff[obj->chaff_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_chaff_Assembly(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Assembly *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_chaff_Assembly(Assembly * obj) 
{
    int i;   


    for(i=0;i<obj->chaff_len;i++)    { /*for i over list length*/ 
      if( obj->chaff[i] != NULL) {  
        free_AssemblySequence(obj->chaff[i]);    
        obj->chaff[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->chaff_len = 0;  
    return i;    
}    


/* Function:  Assembly_alloc_std(void)
 *
 * Descrip:    Equivalent to Assembly_alloc_len(AssemblyLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Assembly *]
 *
 */
Assembly * Assembly_alloc_std(void) 
{
    return Assembly_alloc_len(AssemblyLISTLENGTH);   
}    


/* Function:  Assembly_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Assembly *]
 *
 */
Assembly * Assembly_alloc_len(int len) 
{
    Assembly * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Assembly_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->contig = (AssemblyContig ** ) ckcalloc (len,sizeof(AssemblyContig *))) == NULL) {  
      warn("Warning, ckcalloc failed in Assembly_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->chaff = (AssemblySequence ** ) ckcalloc (len,sizeof(AssemblySequence *))) == NULL)  {  
      warn("Warning, ckcalloc failed in Assembly_alloc_len");    
      return NULL;   
      }  
    out->chaff_len = 0;  
    out->chaff_maxlen = len; 


    return out;  
}    


/* Function:  hard_link_Assembly(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Assembly *]
 *
 * Return [UNKN ]  Undocumented return value [Assembly *]
 *
 */
Assembly * hard_link_Assembly(Assembly * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Assembly object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Assembly_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Assembly *]
 *
 */
Assembly * Assembly_alloc(void) 
{
    Assembly * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Assembly *) ckalloc (sizeof(Assembly))) == NULL)    {  
      warn("Assembly_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->contig = NULL;  
    out->len = out->maxlen = 0;  
    out->chaff = NULL;   
    out->chaff_len = out->chaff_maxlen = 0;  


    return out;  
}    


/* Function:  free_Assembly(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Assembly *]
 *
 * Return [UNKN ]  Undocumented return value [Assembly *]
 *
 */
Assembly * free_Assembly(Assembly * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Assembly obj. Should be trappable");  
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
          free_AssemblyContig(obj->contig[i]);   
        }  
      ckfree(obj->contig);   
      }  
    if( obj->chaff != NULL)  {  
      for(i=0;i<obj->chaff_len;i++)  {  
        if( obj->chaff[i] != NULL)   
          free_AssemblySequence(obj->chaff[i]);  
        }  
      ckfree(obj->chaff);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_AssemblyOutputPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AssemblyOutputPara *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOutputPara *]
 *
 */
AssemblyOutputPara * hard_link_AssemblyOutputPara(AssemblyOutputPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AssemblyOutputPara object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AssemblyOutputPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOutputPara *]
 *
 */
AssemblyOutputPara * AssemblyOutputPara_alloc(void) 
{
    AssemblyOutputPara * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AssemblyOutputPara *) ckalloc (sizeof(AssemblyOutputPara))) == NULL)    {  
      warn("AssemblyOutputPara_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->min_size = 0;   


    return out;  
}    


/* Function:  free_AssemblyOutputPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AssemblyOutputPara *]
 *
 * Return [UNKN ]  Undocumented return value [AssemblyOutputPara *]
 *
 */
AssemblyOutputPara * free_AssemblyOutputPara(AssemblyOutputPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AssemblyOutputPara obj. Should be trappable");    
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

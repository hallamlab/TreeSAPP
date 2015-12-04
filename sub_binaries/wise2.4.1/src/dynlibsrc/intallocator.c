#ifdef _cplusplus
extern "C" {
#endif
#include "intallocator.h"


/* Function:  new_IntAllocatorSet(max_size)
 *
 * Descrip:    Makes a new IntAllocatorSet up to a certain size
 *
 *
 * Arg:        max_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocatorSet *]
 *
 */
# line 43 "intallocator.dy"
IntAllocatorSet * new_IntAllocatorSet(int max_size)
{
  IntAllocatorSet * out;


  out = IntAllocatorSet_alloc();

  out->allocator_set = calloc(max_size,sizeof(IntAllocator*));
  out->max_size = max_size;
  
  return out;

}

/* Function:  realloc_intarray_IntAllocatorSet(ias,current,old_size,new_size)
 *
 * Descrip:    reallocates a piece of memory
 *
 *
 * Arg:             ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 * Arg:         current [UNKN ] Undocumented argument [int *]
 * Arg:        old_size [UNKN ] Undocumented argument [int]
 * Arg:        new_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int *]
 *
 */
# line 60 "intallocator.dy"
int * realloc_intarray_IntAllocatorSet(IntAllocatorSet * ias,int * current,int old_size,int new_size)
{
  int i;
  int * new_a;

  assert(ias);
  assert(new_size > old_size);

  new_a = alloc_intarray_IntAllocatorSet(ias,new_size);
  assert(new_a);

  for(i=0;i<old_size;i++) 
    new_a[i] = current[i];

  free_intarray_IntAllocatorSet(ias,current,old_size);

  return new_a;

}

/* Function:  free_intarray_IntAllocatorSet(ias,array,size)
 *
 * Descrip:    Frees a piece of memory in a IntAllocatorSet
 *
 *
 * Arg:          ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 * Arg:        array [UNKN ] Undocumented argument [int *]
 * Arg:         size [UNKN ] Undocumented argument [int]
 *
 */
# line 83 "intallocator.dy"
void free_intarray_IntAllocatorSet(IntAllocatorSet * ias,int * array,int size)
{
  assert(ias);
  assert(array);

  if( size > ias->max_size ) {
    return free(array);
  }

  assert(ias->allocator_set[size]);

  return free_intarray_IntAllocator(ias->allocator_set[size],array);
}


/* Function:  alloc_intarray_IntAllocatorSet(ias,size)
 *
 * Descrip:    Allocates a new piece of memory in a IntAllocatorSet
 *
 *
 * Arg:         ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 * Arg:        size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int *]
 *
 */
# line 101 "intallocator.dy"
int * alloc_intarray_IntAllocatorSet(IntAllocatorSet * ias,int size)
{
  assert(ias);
  if( size > ias->max_size ) {
    return (int*)calloc(size,sizeof(int));
  }

  if( ias->allocator_set[size] == NULL ) {
    ias->allocator_set[size] = new_IntAllocator(size);
  }

  return alloc_intarray_IntAllocator(ias->allocator_set[size]);

}


/* Function:  new_IntAllocator(size)
 *
 * Descrip:    Makes a new int allocator
 *
 *
 * Arg:        size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocator *]
 *
 */
# line 120 "intallocator.dy"
IntAllocator * new_IntAllocator(int size)
{
  IntAllocator * out;

#ifdef IntAllocator_PARANOIA 
  warn("IntAllocator in paranoia mode. Will perform horrendously slowly due to cycle checking");
#endif

  out = IntAllocator_alloc();
  out->size = size;
  out->start_of_free = NULL;
  out->allocated_blocks = NULL;
  out->max_allocated_blocks = 0;
  out->current_allocated_block = 0;

  return out;
}

#ifdef IntAllocator_PARANOIA

/* Function:  is_acyclic_IntAllocator(ia)
 *
 * Descrip:    Detect cycle
 *
 *
 * Arg:        ia [UNKN ] Undocumented argument [IntAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 143 "intallocator.dy"
boolean is_acyclic_IntAllocator(IntAllocator * ia)
{
  GHashTable * gh;
  IntAllocatorHeader * h;
  int count;

  int dummy;

  assert(ia);

  gh = g_hash_table_new(g_direct_hash,g_direct_equal);

  count = 0;
  for(h = ia->start_of_free;h != NULL;h = h->s.next) {
    if( g_hash_table_lookup(gh,(gconstpointer)h) != NULL ) {
      warn("Found cycle at memory position %d, count %d",h,count);
      return FALSE;
    } else {
      g_hash_table_insert(gh,(gpointer)h,&dummy);
      count++;
    }
  }

  g_hash_table_destroy(gh);

  return TRUE;

}

#endif


/* Function:  show_allocator_status_IntAllocatorSet(ias,ofp)
 *
 * Descrip:    Show status of intallocator set
 *
 *
 * Arg:        ias [UNKN ] Undocumented argument [IntAllocatorSet *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 178 "intallocator.dy"
void show_allocator_status_IntAllocatorSet(IntAllocatorSet * ias,FILE * ofp)
{
  int i;
  int count;
  IntAllocatorHeader * h;
  IntAllocator * ia;
  int mem;
  long total = 0;



  for(i=0;i<ias->max_size;i++) {
    if( ias->allocator_set[i] == NULL ) {
      fprintf(ofp,"[%4d] No allocator\n",i);
    } else {
      ia = ias->allocator_set[i];
      for(h = ia->start_of_free,count =0;h != NULL;h = h->s.next) {
	count++;
      }
      mem = ia->current_allocated_block * (sizeof(IntAllocatorHeader)+(sizeof(int)*ia->size)) * IntAllocator_BLOCKSIZE;
      total += mem;
      fprintf(ofp,"[%4d] %d allocated, %d free, total bytes %d\n",i,ia->current_allocated_block*IntAllocator_BLOCKSIZE,count,mem);
      
    }
  }

  fprintf(ofp,"In total, %ld bytes allocated\n",total);

}

/* Function:  show_allocator_status_IntAllocator(ia,ofp)
 *
 * Descrip:    Shows allocator status
 *
 *
 * Arg:         ia [UNKN ] Undocumented argument [IntAllocator *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 211 "intallocator.dy"
void show_allocator_status_IntAllocator(IntAllocator * ia,FILE * ofp)
{
  int count = 0;
  IntAllocatorHeader * h;

  fprintf(ofp,"%d blocks allocated, using %ld bytes\n",ia->current_allocated_block,ia->current_allocated_block * (sizeof(IntAllocatorHeader)+(sizeof(int)*ia->size)) * IntAllocator_BLOCKSIZE);

  for(h = ia->start_of_free;h != NULL;h = h->s.next) {
    count++;
  }

  fprintf(ofp,"%d units free, %.2f %% occupancy\n",count,(count*100.0)/(ia->current_allocated_block * IntAllocator_BLOCKSIZE));
}

/* Function:  free_intarray_IntAllocator(ia,array)
 *
 * Descrip:    returns an integer back to the pool. NOTE:
 *             This integer * must have come from the pool otherwise
 *             there is going to be a disaster...
 *
 *
 * Arg:           ia [UNKN ] Undocumented argument [IntAllocator *]
 * Arg:        array [UNKN ] Undocumented argument [int *]
 *
 */
# line 230 "intallocator.dy"
void free_intarray_IntAllocator(IntAllocator * ia,int * array)
{
  char * runner;
  IntAllocatorHeader * h;

  assert(ia);
  assert(array);

  runner = (char*) array;
  runner = runner - sizeof(IntAllocatorHeader);

  h = (IntAllocatorHeader*) runner;

  h->s.next = ia->start_of_free;
  ia->start_of_free = h;

#ifdef IntAllocator_PARANOIA
  if( is_acyclic_IntAllocator(ia) == FALSE ) {
    fatal("cycle detected on freeing block for %d position",h);
  }
#endif

  return;

}

/* Function:  alloc_intarray_IntAllocator(ia)
 *
 * Descrip:    returns an integer * from this allocator
 *
 *
 * Arg:        ia [UNKN ] Undocumented argument [IntAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [int *]
 *
 */
# line 259 "intallocator.dy"
int * alloc_intarray_IntAllocator(IntAllocator * ia)
{
  IntAllocatorHeader * h;
  char * runner;

  assert(ia);

  if( ia->start_of_free == NULL ) {
    ia->start_of_free = allocate_new_block_IntAllocator(ia);
  }

#ifdef IntAllocator_PARANOIA
  if( is_acyclic_IntAllocator(ia) == FALSE ) {
    fatal("cycle detected on allocating new block");
  }
#endif


  h = ia->start_of_free;
  ia->start_of_free = h->s.next;

#ifdef IntAllocator_PARANOIA
  if( is_acyclic_IntAllocator(ia) == FALSE ) {
    fatal("cycle detected on returning new block");
  }
#endif

  runner = (char*) h;
  runner = runner + sizeof(IntAllocatorHeader);

  return (int*) runner;
}

/* Function:  allocate_new_block_IntAllocator(ia)
 *
 * Descrip:    internal function to allocate and segment
 *             a block read for use, storing the memory
 *             and segmenting it correctly. Returned pointer
 *             is the first header block
 *
 *
 * Arg:        ia [UNKN ] Undocumented argument [IntAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocatorHeader *]
 *
 */
# line 298 "intallocator.dy"
IntAllocatorHeader * allocate_new_block_IntAllocator(IntAllocator * ia)
{
  int i;
  int step_size;
  char * new_block;
  IntAllocatorHeader * h;
  char * runner;

  assert(ia);

  step_size = sizeof(IntAllocatorHeader)+(sizeof(int)*ia->size);
  new_block = calloc(IntAllocator_BLOCKSIZE,step_size);

  add_new_block_to_memory_handlers_IA(ia,new_block);

  for(i=0;i<IntAllocator_BLOCKSIZE;i++) {
    runner = new_block+(i*step_size); /* in bytes first */
    h = (IntAllocatorHeader*) runner;

    if( i+1 < IntAllocator_BLOCKSIZE ) {
      runner = new_block+((i+1)*step_size);
      h->s.next = (IntAllocatorHeader*) runner;
    } else {
      /* last one */
      h->s.next = NULL;
    }
  }


  return (IntAllocatorHeader*) new_block;
}

/* Function:  add_new_block_to_memory_handlers_IA(ia,new_block)
 *
 * Descrip:    internal function to ensure new block is added, with growth of block array
 *             if needed
 *
 *
 * Arg:               ia [UNKN ] Undocumented argument [IntAllocator *]
 * Arg:        new_block [UNKN ] Undocumented argument [void *]
 *
 */
# line 334 "intallocator.dy"
void add_new_block_to_memory_handlers_IA(IntAllocator * ia,void * new_block)
{
  assert(ia);
  assert(new_block);

  if( ia->allocated_blocks == NULL ) {
    ia->allocated_blocks = calloc(IntAllocator_MEMORY_BLOCK_SIZE,sizeof(void*));
    ia->max_allocated_blocks = IntAllocator_MEMORY_BLOCK_SIZE;
    ia->current_allocated_block = 0;
  } else if( ia->current_allocated_block >= ia->max_allocated_blocks ) {
    ia->allocated_blocks = realloc(ia->allocated_blocks,(ia->max_allocated_blocks+IntAllocator_MEMORY_BLOCK_SIZE)*sizeof(void*));
    ia->max_allocated_blocks = (ia->max_allocated_blocks+IntAllocator_MEMORY_BLOCK_SIZE);
  }
      
  ia->allocated_blocks[ia->current_allocated_block++] = new_block;

  return;
}

# line 391 "intallocator.c"
/* Function:  hard_link_IntAllocator(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [IntAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocator *]
 *
 */
IntAllocator * hard_link_IntAllocator(IntAllocator * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a IntAllocator object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  IntAllocator_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [IntAllocator *]
 *
 */
IntAllocator * IntAllocator_alloc(void) 
{
    IntAllocator * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(IntAllocator *) ckalloc (sizeof(IntAllocator))) == NULL)    {  
      warn("IntAllocator_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->size = 0;   
    out->max_allocated_blocks = 0;   
    out->current_allocated_block = 0;    


    return out;  
}    


/* Function:  free_IntAllocator(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [IntAllocator *]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocator *]
 *
 */
IntAllocator * free_IntAllocator(IntAllocator * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a IntAllocator obj. Should be trappable");  
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
    /* obj->start_of_free is linked in */ 
    /* obj->allocated_blocks is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_IntAllocatorSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [IntAllocatorSet *]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocatorSet *]
 *
 */
IntAllocatorSet * hard_link_IntAllocatorSet(IntAllocatorSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a IntAllocatorSet object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  IntAllocatorSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [IntAllocatorSet *]
 *
 */
IntAllocatorSet * IntAllocatorSet_alloc(void) 
{
    IntAllocatorSet * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(IntAllocatorSet *) ckalloc (sizeof(IntAllocatorSet))) == NULL)  {  
      warn("IntAllocatorSet_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->allocator_set = NULL;   
    out->max_size = 0;   


    return out;  
}    


/* Function:  free_IntAllocatorSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [IntAllocatorSet *]
 *
 * Return [UNKN ]  Undocumented return value [IntAllocatorSet *]
 *
 */
IntAllocatorSet * free_IntAllocatorSet(IntAllocatorSet * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a IntAllocatorSet obj. Should be trappable");   
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
    if( obj->allocator_set != NULL)  
      ckfree(obj->allocator_set);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

#ifdef _cplusplus
extern "C" {
#endif
#include "staticseq.h"


/* Function:  new_StaticSeqHolder(void)
 *
 * Descrip:    makes a new StaticSeqHolder
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StaticSeqHolder *]
 *
 */
# line 25 "staticseq.dy"
StaticSeqHolder * new_StaticSeqHolder(void)
{
  StaticSeqHolder * out;

  out = StaticSeqHolder_alloc();
  assert(out);
  out->gstring_chunk = g_string_chunk_new(1024*1024);

  return out;
}

/* Function:  free_GStringChunk(gs)
 *
 * Descrip:    for registering glib thingy for freeing
 *
 *
 * Arg:        gs [UNKN ] Undocumented argument [GStringChunk *]
 *
 * Return [UNKN ]  Undocumented return value [GStringChunk *]
 *
 */
# line 39 "staticseq.dy"
GStringChunk * free_GStringChunk(GStringChunk * gs)
{
  assert(gs);

  g_string_chunk_free(gs);
  return NULL;
}

/* Function:  new_Sequence_StaticSeqHolder(ssh,seq)
 *
 * Descrip:    Making a new sequence from a staticseq holder - consumes
 *             the sequence object (actually recycling the shell of it)
 *
 *
 * Arg:        ssh [UNKN ] Undocumented argument [StaticSeqHolder *]
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence *]
 *
 */
# line 51 "staticseq.dy"
Sequence * new_Sequence_StaticSeqHolder(StaticSeqHolder * ssh,Sequence * seq)
{
  char * str;
  
  assert(seq);
  assert(ssh);
  assert(ssh->gstring_chunk);
  
  str = g_string_chunk_insert(ssh->gstring_chunk,seq->seq);
  ckfree(seq->seq);
  seq->seq = str;

  seq->dynamite_hard_link++;

  return seq;
}

# line 70 "staticseq.c"
/* Function:  hard_link_StaticSeqHolder(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StaticSeqHolder *]
 *
 * Return [UNKN ]  Undocumented return value [StaticSeqHolder *]
 *
 */
StaticSeqHolder * hard_link_StaticSeqHolder(StaticSeqHolder * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a StaticSeqHolder object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  StaticSeqHolder_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StaticSeqHolder *]
 *
 */
StaticSeqHolder * StaticSeqHolder_alloc(void) 
{
    StaticSeqHolder * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(StaticSeqHolder *) ckalloc (sizeof(StaticSeqHolder))) == NULL)  {  
      warn("StaticSeqHolder_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->gstring_chunk = NULL;   


    return out;  
}    


/* Function:  free_StaticSeqHolder(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StaticSeqHolder *]
 *
 * Return [UNKN ]  Undocumented return value [StaticSeqHolder *]
 *
 */
StaticSeqHolder * free_StaticSeqHolder(StaticSeqHolder * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a StaticSeqHolder obj. Should be trappable");   
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
    if( obj->gstring_chunk != NULL)  
      free_GStringChunk(obj->gstring_chunk);     


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

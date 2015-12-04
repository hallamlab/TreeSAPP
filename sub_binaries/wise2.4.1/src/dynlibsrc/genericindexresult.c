#ifdef _cplusplus
extern "C" {
#endif
#include "genericindexresult.h"


/* Function:  next_interface_GenericIndexResult(data,prev)
 *
 * Descrip:    For interface, returns next position in SeqLookup Results
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 * Arg:        prev [UNKN ] Undocumented argument [SeqLookupResultStruct *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupResultStruct *]
 *
 */
# line 23 "genericindexresult.dy"
SeqLookupResultStruct * next_interface_GenericIndexResult(void * data,SeqLookupResultStruct * prev) 
{
  GenericIndexResult * gir = (GenericIndexResult*)data;

  return &(gir->result[gir->current_pos++]);
}

/* Function:  is_more_interface_GenericIndexResult(data)
 *
 * Descrip:    For interface, indicates whether there is more stuff to find or not
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 33 "genericindexresult.dy"
boolean is_more_interface_GenericIndexResult(void * data)
{
  GenericIndexResult * gir = (GenericIndexResult*)data;

  if( gir->current_pos < gir->len ) {
    return TRUE;
  } 
  return FALSE;
}

/* Function:  free_noop_GenericIndexResult(data)
 *
 * Descrip:    Frees the data
 *
 *
 * Arg:        data [UNKN ] Undocumented argument [void *]
 *
 */
# line 46 "genericindexresult.dy"
void free_noop_GenericIndexResult(void * data)
{
  return;
}

/* Function:  add_GenericIndexResult(gir,seq,pos)
 *
 * Descrip:    Adds another result to a IndexResult
 *
 *
 * Arg:        gir [UNKN ] Undocumented argument [GenericIndexResult *]
 * Arg:        seq [UNKN ] Undocumented argument [Sequence *]
 * Arg:        pos [UNKN ] Undocumented argument [int]
 *
 */
# line 54 "genericindexresult.dy"
void add_GenericIndexResult(GenericIndexResult * gir,Sequence * seq,int pos)
{
  assert(gir);
  assert(seq);

  if( gir->len >= gir->max_len ) {
    gir->result = realloc(gir->result,2*gir->max_len*sizeof(SeqLookupResultStruct));
    gir->max_len = 2* gir->max_len;
  }
  
  gir->result[gir->len].seq = seq;
  gir->result[gir->len].pos = pos;
  
  gir->len++;

}


/* Function:  free_GenericIndexResult(p)
 *
 * Descrip:    Frees GenericIndexResults - overrides dynamite default
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [GenericIndexResult *]
 *
 * Return [UNKN ]  Undocumented return value [GenericIndexResult *]
 *
 */
# line 76 "genericindexresult.dy"
GenericIndexResult * free_GenericIndexResult(GenericIndexResult * p)
{
  assert(p);
  free(p->result);
  free(p);

  return NULL;
}



# line 101 "genericindexresult.c"
/* Function:  hard_link_GenericIndexResult(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenericIndexResult *]
 *
 * Return [UNKN ]  Undocumented return value [GenericIndexResult *]
 *
 */
GenericIndexResult * hard_link_GenericIndexResult(GenericIndexResult * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenericIndexResult object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenericIndexResult_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenericIndexResult *]
 *
 */
GenericIndexResult * GenericIndexResult_alloc(void) 
{
    GenericIndexResult * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenericIndexResult *) ckalloc (sizeof(GenericIndexResult))) == NULL)    {  
      warn("GenericIndexResult_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->result = NULL;  
    out->len = 0;    
    out->max_len = 0;    
    out->current_pos = 0;    


    return out;  
}    



#ifdef _cplusplus
}
#endif

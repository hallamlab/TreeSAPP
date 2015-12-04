#ifdef _cplusplus
extern "C" {
#endif
#include "proteinindexcons.h"

/* Function:  new_SeqLookupInterface_from_ProteinIndexConstructor(pic)
 *
 * Descrip:    Makes a new SeqLookupInterface from ProteinIndexConstructor
 *
 *
 * Arg:        pic [UNKN ] Undocumented argument [ProteinIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [SeqLookupInterface *]
 *
 */
# line 34 "proteinindexcons.dy"
SeqLookupInterface * new_SeqLookupInterface_from_ProteinIndexConstructor(ProteinIndexConstructor * pic)
{
  switch(pic->type) {

    case ProteinIndexConstructor_Array  :
    return new_ArraySeq_SeqLookupInterface(SEQLOOKUP_5AA_SIZE,1000);
    
    case  ProteinIndexConstructor_Hash :
    return new_ghash_SeqLookupInterface();
    
    case ProteinIndexConstructor_Stream : 
    return new_ProteinStreamedIndex_SeqLookupInterface(pic->waypost);
    
    case ProteinIndexConstructor_Shadow :
    return new_ShadowSequenceIndex_SeqLookupInterface(pic->shadowlength,pic->has_maxlen,pic->max_seqlen,pic->shadow_error);
    
    default:
    fatal("Cannot process type %d as protein index constructor",pic->type);
  }
  return NULL;
}

/* Function:  show_help_ProteinIndexConstructor(ofp)
 *
 * Descrip:    provides help for protein index constructor
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 59 "proteinindexcons.dy"
void show_help_ProteinIndexConstructor(FILE * ofp)
{
  fprintf(ofp,"Protein Index construction options\n");
  fprintf(ofp,"   -pitype [array/hash/stream/shadow] - default array\n");
  fprintf(ofp,"   -piwaypost [number]  - waypost for streamed cases, default 3\n");
  fprintf(ofp,"   -pishadow [number]   - shadow length for shadow cases, default 15\n");
  fprintf(ofp,"   -pishadow_err [number] - errors per 100 identities tolerated, 3\n");
  fprintf(ofp,"   -piseqmax            - indexes can assumme maximum length of seq\n");
  fprintf(ofp,"   -piseqmax_len [number] - assummed max sequnce length, default 1000\n"); 

  return;
}


/* Function:  new_ProteinIndexConstructor_from_argv(argc,argv)
 *
 * Descrip:    Provides a ProteinIndexConstructor argument from argv
 *
 *
 * Arg:        argc [UNKN ] Undocumented argument [int *]
 * Arg:        argv [UNKN ] Undocumented argument [char **]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIndexConstructor *]
 *
 */
# line 76 "proteinindexcons.dy"
ProteinIndexConstructor * new_ProteinIndexConstructor_from_argv(int * argc,char ** argv)
{
  ProteinIndexConstructor * out;
  char * temp;

  out = ProteinIndexConstructor_alloc();

  temp = strip_out_assigned_argument(argc,argv,"pitype");
  if( temp != NULL ) {
    if( strcmp(temp,"array") == 0 ) {
      out->type = ProteinIndexConstructor_Array;
    } else if ( strcmp(temp,"hash") == 0 ) {
      out->type = ProteinIndexConstructor_Hash;
    } else if ( strcmp(temp,"stream") == 0 ) {
      out->type = ProteinIndexConstructor_Stream;
    } else if ( strcmp(temp,"shadow") == 0 ) {
      out->type = ProteinIndexConstructor_Shadow;
    } else {
      fatal("Could not interpret %s as a protein index type",temp);
    }
  }

  strip_out_integer_argument(argc,argv,"piwaypost",&out->waypost);
  strip_out_integer_argument(argc,argv,"pishadow",&out->shadowlength);

  strip_out_boolean_def_argument(argc,argv,"piseqmax",&out->waypost);
  strip_out_integer_argument(argc,argv,"piseqmax_len",&out->shadowlength);


  return out;

}

# line 100 "proteinindexcons.c"
/* Function:  hard_link_ProteinIndexConstructor(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ProteinIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIndexConstructor *]
 *
 */
ProteinIndexConstructor * hard_link_ProteinIndexConstructor(ProteinIndexConstructor * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ProteinIndexConstructor object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ProteinIndexConstructor_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinIndexConstructor *]
 *
 */
ProteinIndexConstructor * ProteinIndexConstructor_alloc(void) 
{
    ProteinIndexConstructor * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ProteinIndexConstructor *) ckalloc (sizeof(ProteinIndexConstructor))) == NULL)  {  
      warn("ProteinIndexConstructor_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = ProteinIndexConstructor_Array;   
    out->waypost = 3;    
    out->shadowlength = 15;  
    out->has_maxlen = 0; 
    out->max_seqlen = 1000;  
    out->shadow_error = 3;   


    return out;  
}    


/* Function:  free_ProteinIndexConstructor(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ProteinIndexConstructor *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinIndexConstructor *]
 *
 */
ProteinIndexConstructor * free_ProteinIndexConstructor(ProteinIndexConstructor * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ProteinIndexConstructor obj. Should be trappable");   
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

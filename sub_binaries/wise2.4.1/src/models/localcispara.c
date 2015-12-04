#ifdef _cplusplus
extern "C" {
#endif
#include "localcispara.h"

# line 41 "localcispara.dy"
LocalCisHitScore * standard_LocalCisHitScore(NMaskType nmask)
{
  LocalCisHitProb *p;
  LocalCisHitScore * out;

  p = standard_LocalCisHitProb(nmask);
  out = LocalCisHitScore_from_LocalCisHitProb(p);

  free_LocalCisHitProb(p);

  return out;
}

# line 54 "localcispara.dy"
LocalCisHitProb * standard_LocalCisHitProb(NMaskType nmask)
{
  LocalCisHitProb * out;

  out = LocalCisHitProb_alloc();

  out->comp65 = DnaProbMatrix_from_match(0.65,nmask);  
  flat_null_DnaProbMatrix(out->comp65);  

  out->comp75 = DnaProbMatrix_from_match(0.75,nmask);  
  flat_null_DnaProbMatrix(out->comp75);  

  out->comp85 = DnaProbMatrix_from_match(0.85,nmask);  
  flat_null_DnaProbMatrix(out->comp85);  

  out->comp95 = DnaProbMatrix_from_match(0.95,nmask);  
  flat_null_DnaProbMatrix(out->comp95);  

  out->g = LCH_GAP / LCH_UNMATCHED_PEN;
  out->u = 1.0;
  out->v = (1.0- LCH_UNMATCHED_PEN)/(LCH_UNMATCHED_PEN);
  out->s = (1.0 - (LCH_GAP + LCH_GAP + LCH_BLOCKOPEN))/(LCH_UNMATCHED_PEN * LCH_UNMATCHED_PEN);
  out->b = (LCH_BLOCKOPEN/LCH_UNMATCHED_PEN);

  return out;
}


# line 82 "localcispara.dy"
LocalCisHitScore * LocalCisHitScore_from_LocalCisHitProb(LocalCisHitProb * lchp)
{
  LocalCisHitScore * lchs;

  lchs = LocalCisHitScore_alloc();

  lchs->comp65 = DnaMatrix_from_DnaProbMatrix(lchp->comp65);
  lchs->comp75 = DnaMatrix_from_DnaProbMatrix(lchp->comp75);
  lchs->comp85 = DnaMatrix_from_DnaProbMatrix(lchp->comp85);
  lchs->comp95 = DnaMatrix_from_DnaProbMatrix(lchp->comp95);
  lchs->g = Probability2Score(lchp->g);
  lchs->u = Probability2Score(lchp->u);
  lchs->v = Probability2Score(lchp->v);
  lchs->s = Probability2Score(lchp->s);
  lchs->b = Probability2Score(lchp->b);

  return lchs;
}


# line 69 "localcispara.c"
/* Function:  hard_link_LocalCisHitProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocalCisHitProb *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitProb *]
 *
 */
LocalCisHitProb * hard_link_LocalCisHitProb(LocalCisHitProb * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LocalCisHitProb object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LocalCisHitProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitProb *]
 *
 */
LocalCisHitProb * LocalCisHitProb_alloc(void) 
{
    LocalCisHitProb * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LocalCisHitProb *) ckalloc (sizeof(LocalCisHitProb))) == NULL)  {  
      warn("LocalCisHitProb_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->comp65 = NULL;  
    out->comp75 = NULL;  
    out->comp85 = NULL;  
    out->comp95 = NULL;  
    out->g = 0.0;    
    out->u = 0.0;    
    out->v = 0.0;    
    out->s = 0.0;    
    out->b = 0.0;    


    return out;  
}    


/* Function:  free_LocalCisHitProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalCisHitProb *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitProb *]
 *
 */
LocalCisHitProb * free_LocalCisHitProb(LocalCisHitProb * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LocalCisHitProb obj. Should be trappable");   
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
    if( obj->comp65 != NULL) 
      free_DnaProbMatrix(obj->comp65);   
    if( obj->comp75 != NULL) 
      free_DnaProbMatrix(obj->comp75);   
    if( obj->comp85 != NULL) 
      free_DnaProbMatrix(obj->comp85);   
    if( obj->comp95 != NULL) 
      free_DnaProbMatrix(obj->comp95);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_LocalCisHitScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocalCisHitScore *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitScore *]
 *
 */
LocalCisHitScore * hard_link_LocalCisHitScore(LocalCisHitScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LocalCisHitScore object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LocalCisHitScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitScore *]
 *
 */
LocalCisHitScore * LocalCisHitScore_alloc(void) 
{
    LocalCisHitScore * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LocalCisHitScore *) ckalloc (sizeof(LocalCisHitScore))) == NULL)    {  
      warn("LocalCisHitScore_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->comp65 = NULL;  
    out->comp75 = NULL;  
    out->comp85 = NULL;  
    out->comp95 = NULL;  
    out->g = 0;  
    out->u = 0;  
    out->v = 0;  
    out->s = 0;  
    out->b = 0;  


    return out;  
}    


/* Function:  free_LocalCisHitScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalCisHitScore *]
 *
 * Return [UNKN ]  Undocumented return value [LocalCisHitScore *]
 *
 */
LocalCisHitScore * free_LocalCisHitScore(LocalCisHitScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LocalCisHitScore obj. Should be trappable");  
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
    if( obj->comp65 != NULL) 
      free_DnaMatrix(obj->comp65);   
    if( obj->comp75 != NULL) 
      free_DnaMatrix(obj->comp75);   
    if( obj->comp85 != NULL) 
      free_DnaMatrix(obj->comp85);   
    if( obj->comp95 != NULL) 
      free_DnaMatrix(obj->comp95);   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

#ifdef _cplusplus
extern "C" {
#endif
#include "cdparser.h"

/* Function:  removed_probability_from_cds_cdna(cdp)
 *
 * Descrip:    Makes a convienient sum over all the transition
 *             probabilities
 *
 *
 * Arg:        cdp [UNKN ] Undocumented argument [cDNAParser *]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
# line 52 "cdparser.dy"
Probability removed_probability_from_cds_cdna(cDNAParser * cdp)
{
  return cdp->trans[PCD_INSERT_2_BASE] +
    cdp->trans[PCD_INSERT_1_BASE] +
      cdp->trans[PCD_DELETE_2_BASE] +
	cdp->trans[PCD_DELETE_1_BASE];
}
 
/* Function:  cDNAParserScore_from_cDNAParser(cdp)
 *
 * Descrip:    Makes a new Score object from its probability
 *             counterpart
 *
 *
 * Arg:        cdp [UNKN ] Undocumented argument [cDNAParser *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParserScore *]
 *
 */
# line 64 "cdparser.dy"
cDNAParserScore * cDNAParserScore_from_cDNAParser(cDNAParser * cdp)
{
  cDNAParserScore * out;

  out = cDNAParserScore_alloc();

  Probability2Score_move(cdp->trans,out->trans,PCD_PARSER_TRANS_LEN);

  return out;
}

/* Function:  flat_cDNAParser(p)
 *
 * Descrip:    Makes a flat (ie, indels of 1 or 2 == p)
 *             cDNA parser. This means that insertions
 *             and deletions of both 1 or 2 bases are
 *             all parameterised at the same probability
 *
 *
 *
 * Arg:        p [READ ] probability of an indel [Probability]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParser *]
 *
 */
# line 84 "cdparser.dy"
cDNAParser * flat_cDNAParser(Probability p)
{
  cDNAParser * out;

  out = cDNAParser_alloc();

  out->trans[PCD_INSERT_2_BASE] = p;
  out->trans[PCD_INSERT_1_BASE] = p;
  out->trans[PCD_DELETE_2_BASE] = p;
  out->trans[PCD_DELETE_1_BASE] = p;
    
  return out;
}


# line 72 "cdparser.c"
/* Function:  hard_link_cDNAParser(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [cDNAParser *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParser *]
 *
 */
cDNAParser * hard_link_cDNAParser(cDNAParser * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a cDNAParser object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  cDNAParser_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [cDNAParser *]
 *
 */
cDNAParser * cDNAParser_alloc(void) 
{
    cDNAParser * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(cDNAParser *) ckalloc (sizeof(cDNAParser))) == NULL)    {  
      warn("cDNAParser_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* trans[PCD_PARSER_TRANS_LEN] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_cDNAParser(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [cDNAParser *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParser *]
 *
 */
cDNAParser * free_cDNAParser(cDNAParser * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a cDNAParser obj. Should be trappable");    
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


/* Function:  hard_link_cDNAParserScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParserScore *]
 *
 */
cDNAParserScore * hard_link_cDNAParserScore(cDNAParserScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a cDNAParserScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  cDNAParserScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [cDNAParserScore *]
 *
 */
cDNAParserScore * cDNAParserScore_alloc(void) 
{
    cDNAParserScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(cDNAParserScore *) ckalloc (sizeof(cDNAParserScore))) == NULL)  {  
      warn("cDNAParserScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* trans[PCD_PARSER_TRANS_LEN] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_cDNAParserScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [cDNAParserScore *]
 *
 * Return [UNKN ]  Undocumented return value [cDNAParserScore *]
 *
 */
cDNAParserScore * free_cDNAParserScore(cDNAParserScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a cDNAParserScore obj. Should be trappable");   
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

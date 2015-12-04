#ifdef _cplusplus
extern "C" {
#endif
#include "geneparser4.h"


# line 34 "geneparser4.dy"
GeneParser4 * std_GeneParser4(double indel,double intron2cds)
{
  int i;
  GeneParser4 * out;

  out = GeneParser4_alloc();

  out->transition[GP4_INTRON2CDS] = intron2cds;

  out->transition[GP4_INTRON2INTRON]    = 1.0;
  out->transition[GP4_DELETE_1_BASE]    = indel;
  out->transition[GP4_DELETE_2_BASE]    = indel;
  out->transition[GP4_INSERT_1_BASE]    = indel;
  out->transition[GP4_INSERT_2_BASE]    = indel;
  out->transition[GP4_LOOP2MODEL]       = 1.0;
  out->transition[GP4_LOOP2LOOP]        = 0;


  for(i=0;i<5;i++)
    out->intron[i] = 1.0;

  return out;
}



# line 60 "geneparser4.dy"
GeneParser4Score * GeneParser4Score_from_GeneParser21Score(GeneParser21Score * gp21s)
{
  int i;
  GeneParser4Score * out;

  out = GeneParser4Score_alloc();

  out->transition[GP4_INTRON2CDS] = gp21s->transition[GP21_CENTRAL2PY] + gp21s->transition[GP21_PY2SPACER] + gp21s->transition[GP21_SPACER2CDS];

  out->transition[GP4_INTRON2INTRON]    = gp21s->transition[GP21_CENTRAL2CENTRAL];
  out->transition[GP4_DELETE_1_BASE]    = gp21s->transition[GP21_DELETE_1_BASE];
  out->transition[GP4_DELETE_2_BASE]    = gp21s->transition[GP21_DELETE_2_BASE];
  out->transition[GP4_INSERT_1_BASE]    = gp21s->transition[GP21_INSERT_1_BASE];
  out->transition[GP4_INSERT_2_BASE]    = gp21s->transition[GP21_INSERT_2_BASE];
  out->transition[GP4_LOOP2MODEL]       = gp21s->transition[GP21_RND2MODEL];
  /*  out->transition[GP4_LOOP2LOOP]        = gp21s->transition[GP21_RND2RND]; */
  out->transition[GP4_LOOP2LOOP]        = 0;
  /*  fprintf(stderr,"Loop score is %d\n",out->transition[GP4_LOOP2LOOP]); */


  for(i=0;i<5;i++)
    out->intron[i] = gp21s->central[i];

  return out;
}
  
  

# line 88 "geneparser4.dy"
GeneParser4Score * GeneParser4Score_from_GeneParser4(GeneParser4 * gp4)
{
  GeneParser4Score * out;

  out = GeneParser4Score_alloc();

  Probability2Score_move(gp4->transition,out->transition,GP4_TRANSITION_LEN);
  Probability2Score_move(gp4->intron,out->intron,5);

  return out;
}


# line 76 "geneparser4.c"
/* Function:  hard_link_GeneParser4(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneParser4 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4 *]
 *
 */
GeneParser4 * hard_link_GeneParser4(GeneParser4 * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneParser4 object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneParser4_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4 *]
 *
 */
GeneParser4 * GeneParser4_alloc(void) 
{
    GeneParser4 * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneParser4 *) ckalloc (sizeof(GeneParser4))) == NULL)  {  
      warn("GeneParser4_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* transition[GP4_TRANSITION_LEN] is an array: no default possible */ 
    /* intron[5] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_GeneParser4(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneParser4 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4 *]
 *
 */
GeneParser4 * free_GeneParser4(GeneParser4 * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneParser4 obj. Should be trappable");   
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


/* Function:  hard_link_GeneParser4Score(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneParser4Score *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4Score *]
 *
 */
GeneParser4Score * hard_link_GeneParser4Score(GeneParser4Score * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneParser4Score object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneParser4Score_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4Score *]
 *
 */
GeneParser4Score * GeneParser4Score_alloc(void) 
{
    GeneParser4Score * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneParser4Score *) ckalloc (sizeof(GeneParser4Score))) == NULL)    {  
      warn("GeneParser4Score_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* transition[GP4_TRANSITION_LEN] is an array: no default possible */ 
    /* intron[5] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_GeneParser4Score(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneParser4Score *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser4Score *]
 *
 */
GeneParser4Score * free_GeneParser4Score(GeneParser4Score * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneParser4Score obj. Should be trappable");  
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

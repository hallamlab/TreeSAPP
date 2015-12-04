#ifdef _cplusplus
extern "C" {
#endif
#include "geneparser21.h"



# line 64 "geneparser21.dy"
RandomModelDNA * fudged_mixed_RandomModelDNA_from_GeneParser21(GeneParser21 * gp,RandomModelDNA * rnd)
{
  RandomModelDNA * out;
  int i;

  out = RandomModelDNA_alloc();

#define PROBMAX(this,that) (this > that ? this : that)


  for(i=0;i<5;i++) {
    out->base[i] = PROBMAX(PROBMAX(gp->central[i],gp->py[i]),PROBMAX(gp->spacer[i],rnd->base[i]))/rnd->base[i];
  }

#undef PROBMAX

  return out;
}

  
/***
  Adds the error probabilities

***/

# line 89 "geneparser21.dy"
void add_flat_error_probabilities_GeneParser21(GeneParser21 * gp21,Probability error)
{
  add_error_probabilities_GeneParser21(gp21,error,error,error,error);
}

# line 94 "geneparser21.dy"
void add_error_probabilities_GeneParser21(GeneParser21 * gp21,Probability insert_1,Probability insert_2,Probability delete_1,Probability delete_2)
{
  gp21->transition[GP21_INSERT_1_BASE] = insert_1;
  gp21->transition[GP21_INSERT_2_BASE] = insert_2;
  gp21->transition[GP21_DELETE_1_BASE] = delete_1;
  gp21->transition[GP21_DELETE_2_BASE] = delete_2;
}



  
/****
  Makes the GeneParser21 structure, which essentially has
  transition probabilities for geneparsing

  this sets everything to zero and then 
  only allocates gene parsing probs. Frameshifting error is not
  delt with here.

  ****/


# line 116 "geneparser21.dy"
GeneParser21 * GeneParser21_from_GeneFrequency21_cds(GeneFrequency21 * gf,Probability rnd_loop,Probability cds_loop,Probability rnd_to_model,Probability link_loop,Probability link_to_model)
{
  GeneParser21 * out;
  register int i;
  double total;


  out = GeneParser21_alloc();

  set_Probability_array(out->transition,0.0,GENEPARSER21_TRANSITION_LEN);
  set_Probability_array(out->central,0.0,GENEPARSER21_EMISSION_LEN);
  set_Probability_array(out->spacer,0.0,GENEPARSER21_EMISSION_LEN);
  set_Probability_array(out->py,0.0,GENEPARSER21_EMISSION_LEN);


  total = sum_Probability_array(gf->central,4);  


  for(i=0;i<4;i++) {
    out->central[i] = gf->central[i]/total;
  }
  out->central[4] = 1.0;

  total = sum_Probability_array(gf->py,4);  
  for(i=0;i<4;i++) {
    out->py[i] = gf->py[i]/total;
  }
  out->py[4] = 1.0;

  total = sum_Probability_array(gf->spacer,4);  
  for(i=0;i<4;i++) {
    out->spacer[i] = gf->spacer[i]/total;
  }
  out->spacer[4] = 1.0;

  out->transition[GP21_CENTRAL2CENTRAL] = gf->transition[GF21_CENTRAL_STAY];
  out->transition[GP21_CENTRAL2PY]      = 1-  gf->transition[GF21_CENTRAL_STAY];
  out->transition[GP21_PY2PY]           = gf->transition[GF21_PY_STAY];
  out->transition[GP21_PY2SPACER]       = 1 - gf->transition[GF21_PY_STAY] * (1-gf->transition[GF21_NO_SPACER]);
  out->transition[GP21_PY2CDS]       = 1 - gf->transition[GF21_PY_STAY] * (gf->transition[GF21_NO_SPACER]);
  out->transition[GP21_SPACER2SPACER]   = gf->transition[GF21_SPACER_STAY];
  out->transition[GP21_SPACER2CDS]      = 1 - gf->transition[GF21_SPACER_STAY];
  

  out->transition[GP21_CDS2CDS] = cds_loop;
  out->transition[GP21_CDS2RND] = (1-cds_loop);
  out->transition[GP21_RND2RND] = rnd_loop;
  /*  fprintf(stderr,"Score is %f\n",out->transition[GP21_RND2RND]); */
  out->transition[GP21_RND2CDS] = (1-rnd_loop-rnd_to_model);
  out->transition[GP21_RND2MODEL] = rnd_to_model;
  out->transition[GP21_LINK2MODEL] = link_to_model;
  out->transition[GP21_LINK2LINK] = link_loop;
  out->transition[GP21_LINK2RND] = (1- link_loop - link_to_model) ;

  return out;
}

# line 173 "geneparser21.dy"
RandomModelDNA * RandomModelDNA_from_central_GeneParser21(GeneParser21 *gp21)
{
  RandomModelDNA * rmd;

  rmd = RandomModelDNA_alloc();
  if( rmd == NULL ) 
    return NULL;

  rmd->base[BASE_A] = gp21->central[BASE_A];
  rmd->base[BASE_T] = gp21->central[BASE_T];
  rmd->base[BASE_G] = gp21->central[BASE_G];
  rmd->base[BASE_C] = gp21->central[BASE_C];
  rmd->base[BASE_N] = gp21->central[BASE_N];

  return rmd;
}

# line 190 "geneparser21.dy"
void show_GeneParser21(GeneParser21 * gp21,FILE * ofp)
{
  fprintf(ofp,"Central emissions\n");
  show_Probability_array(gp21->central,GENEPARSER21_EMISSION_LEN,ofp);
  fprintf(ofp,"\nPyrimidine emissions\n");
  show_Probability_array(gp21->py,GENEPARSER21_EMISSION_LEN,ofp);
  fprintf(ofp,"\nSpacer emissions\n");
  show_Probability_array(gp21->spacer,GENEPARSER21_EMISSION_LEN,ofp);
  fprintf(ofp,"\nTransitions\n");
  show_Probability_array(gp21->transition,GENEPARSER21_TRANSITION_LEN,ofp);

}


# line 204 "geneparser21.dy"
GeneParser21 * std_GeneParser21(void)
{
  GeneParser21 * out;


  out = GeneParser21_alloc();

  out->central[BASE_T] = 0.25;
  out->central[BASE_G] = 0.25;
  out->central[BASE_C] = 0.25;
  out->central[BASE_A] = 0.25;

  out->py[BASE_T] = 0.4;
  out->py[BASE_C] = 0.4;
  out->py[BASE_G] = 0.1;
  out->py[BASE_A] = 0.1;

  out->spacer[BASE_T] = 0.4;
  out->spacer[BASE_G] = 0.2;
  out->spacer[BASE_C] = 0.2;
  out->spacer[BASE_A] = 0.2;
   
  out->central[BASE_N] = out->py[BASE_N] = out->spacer[BASE_N] = 1.0;
  

  out->transition[GP21_CDS2CENTRAL]     = 0.0;
  out->transition[GP21_CENTRAL2CENTRAL] = 0.998;
  out->transition[GP21_CENTRAL2PY]      = 0.002;
  out->transition[GP21_PY2PY]           = 0.94;
  out->transition[GP21_PY2CDS]          = 0.01;
  out->transition[GP21_PY2SPACER]       = 0.05;
  out->transition[GP21_SPACER2SPACER]   = 0.9;
  out->transition[GP21_SPACER2CDS]      = 0.1;
  out->transition[GP21_INSERT_1_BASE] = 0.001;
  out->transition[GP21_INSERT_2_BASE] = 0.001;
  out->transition[GP21_DELETE_1_BASE] = 0.002;
  out->transition[GP21_DELETE_2_BASE] = 0.003;

  return out;
}


# line 246 "geneparser21.dy"
Probability removed_probability_from_cds(GeneParser21 * gp21)
{

  Probability ret;

  ret =  gp21->transition[GP21_INSERT_1_BASE] +
    gp21->transition[GP21_INSERT_2_BASE] +
      gp21->transition[GP21_DELETE_2_BASE] +
	gp21->transition[GP21_DELETE_1_BASE];

  return ret;
}

# line 259 "geneparser21.dy"
void GeneParser21_fold_in_RandomModelDNA(GeneParser21 * gp21,RandomModelDNA * rmd)
{
  (void)Probability_array_divide(gp21->central,gp21->central,rmd->base,5);
  (void)Probability_array_divide(gp21->py,gp21->py,rmd->base,5);
  (void)Probability_array_divide(gp21->spacer,gp21->spacer,rmd->base,5);
	  
}

# line 267 "geneparser21.dy"
GeneParser21Score * GeneParser21Score_from_GeneParser21(GeneParser21 * gp21)
{
  GeneParser21Score * out;

  out = GeneParser21Score_alloc();

  if( out == NULL )
    return NULL;

  Probability2Score_move(gp21->transition,out->transition,GENEPARSER21_TRANSITION_LEN);
  Probability2Score_move(gp21->central,   out->central,   GENEPARSER21_EMISSION_LEN);
  Probability2Score_move(gp21->py,        out->py,        GENEPARSER21_EMISSION_LEN);
  Probability2Score_move(gp21->spacer,    out->spacer,    GENEPARSER21_EMISSION_LEN);

  return out;
}



# line 239 "geneparser21.c"
/* Function:  hard_link_GeneParser21(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneParser21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21 *]
 *
 */
GeneParser21 * hard_link_GeneParser21(GeneParser21 * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneParser21 object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneParser21_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21 *]
 *
 */
GeneParser21 * GeneParser21_alloc(void) 
{
    GeneParser21 * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneParser21 *) ckalloc (sizeof(GeneParser21))) == NULL)    {  
      warn("GeneParser21_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* transition[GENEPARSER21_TRANSITION_LEN] is an array: no default possible */ 
    /* central[GENEPARSER21_EMISSION_LEN] is an array: no default possible */ 
    /* py[GENEPARSER21_EMISSION_LEN] is an array: no default possible */ 
    /* spacer[GENEPARSER21_EMISSION_LEN] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_GeneParser21(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneParser21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21 *]
 *
 */
GeneParser21 * free_GeneParser21(GeneParser21 * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneParser21 obj. Should be trappable");  
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


/* Function:  hard_link_GeneParser21Score(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GeneParser21Score *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21Score *]
 *
 */
GeneParser21Score * hard_link_GeneParser21Score(GeneParser21Score * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GeneParser21Score object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GeneParser21Score_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21Score *]
 *
 */
GeneParser21Score * GeneParser21Score_alloc(void) 
{
    GeneParser21Score * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GeneParser21Score *) ckalloc (sizeof(GeneParser21Score))) == NULL)  {  
      warn("GeneParser21Score_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* transition[GENEPARSER21_TRANSITION_LEN] is an array: no default possible */ 
    /* central[GENEPARSER21_EMISSION_LEN] is an array: no default possible */ 
    /* py[GENEPARSER21_EMISSION_LEN] is an array: no default possible */ 
    /* spacer[GENEPARSER21_EMISSION_LEN] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_GeneParser21Score(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneParser21Score *]
 *
 * Return [UNKN ]  Undocumented return value [GeneParser21Score *]
 *
 */
GeneParser21Score * free_GeneParser21Score(GeneParser21Score * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GeneParser21Score obj. Should be trappable"); 
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

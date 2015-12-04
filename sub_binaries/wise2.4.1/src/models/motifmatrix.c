#ifdef _cplusplus
extern "C" {
#endif
#include "motifmatrix.h"


# line 46 "motifmatrix.dy"
MotifMatrixScore * MotifMatrixScore_from_MotifMatrixPara(MotifMatrixPara * mmp)
{
  MotifMatrixScore * out;
  DnaProbMatrix * dmp;

  out = MotifMatrixScore_alloc();

  dmp = DnaProbMatrix_from_match(mmp->comp_in_match,NMaskType_BANNED);  
  assert(dmp);
  flat_null_DnaProbMatrix(dmp);  

  out->comp_in_motif = DnaMatrix_from_DnaProbMatrix(dmp);
  free_DnaProbMatrix(dmp);

  dmp = DnaProbMatrix_from_match(mmp->comp_out_match,NMaskType_BANNED);  
  assert(dmp);
  flat_null_DnaProbMatrix(dmp);  

  out->comp_out_motif = DnaMatrix_from_DnaProbMatrix(dmp);
  free_DnaProbMatrix(dmp);


  dmp = DnaProbMatrix_from_match(mmp->comp_spacer,NMaskType_BANNED);  
  assert(dmp);
  flat_null_DnaProbMatrix(dmp);  

  out->comp_spacer = DnaMatrix_from_DnaProbMatrix(dmp);
  free_DnaProbMatrix(dmp);

  out->region_in       = Probability2Score(mmp->region_in);
  out->motif_indel     = Probability2Score(mmp->motif_indel);
  out->cons_indel      = Probability2Score(mmp->cons_indel);
  out->spacer_indel    = Probability2Score(mmp->spacer_indel);
  out->spacer_to_cons  = Probability2Score(mmp->spacer_to_cons);
  out->spacer_to_motif = Probability2Score(mmp->spacer_to_motif);
  out->spacer_duration = Probability2Score(mmp->spacer_duration);
  out->motif_duration  = Probability2Score(mmp->motif_duration);
  out->cons_duration   = Probability2Score(mmp->cons_duration);


  return out;

}

# line 90 "motifmatrix.dy"
MotifMatrixPara * new_MotifMatrixPara_from_argv(int * argc,char ** argv)
{
  MotifMatrixPara * out;

  out = MotifMatrixPara_alloc();

  out->comp_in_match  = 0.9;
  out->comp_out_match = 0.75;
  out->comp_spacer    = 0.35;
  out->region_in      = 0.0001;
  out->motif_indel    = 0.00001;
  out->cons_indel     = 0.025;
  out->spacer_indel   = 0.1;
  out->spacer_duration = 0.2;
  out->motif_duration = 0.9;
  out->cons_duration  = 0.7;
  out->spacer_to_motif = 0.05;
  out->spacer_to_cons  = 0.000001;

  strip_out_float_argument(argc,argv,"mm_motif",&out->comp_in_match);
  strip_out_float_argument(argc,argv,"mm_cons",&out->comp_out_match);
  strip_out_float_argument(argc,argv,"mm_spacer",&out->comp_spacer);
  strip_out_float_argument(argc,argv,"mm_motif_indel",&out->motif_indel);
  strip_out_float_argument(argc,argv,"mm_cons_indel",&out->cons_indel);
  strip_out_float_argument(argc,argv,"mm_spacer_indel",&out->spacer_indel);
  strip_out_float_argument(argc,argv,"mm_switch_motif",&out->spacer_to_motif);
  strip_out_float_argument(argc,argv,"mm_switch_cons",&out->spacer_to_cons);

  return out;
}

# line 121 "motifmatrix.dy"
void show_help_MotifMatrixPara(FILE * ofp)
{
  fprintf(ofp,"Motif Matrix matching paramters\n");
  fprintf(ofp,"  -mm_motif [0.9]  Probability of a match in a motif\n");
  fprintf(ofp,"  -mm_cons  [0.75] Probability of a match in a non-motif conserved\n");
  fprintf(ofp,"  -mm_spacer[0.35] Probability of a match in a spacer\n");
  fprintf(ofp,"  -mm_motif_indel [0.00001] indel inside a motif\n");
  fprintf(ofp,"  -mm_cons_indel  [0.025]   indel inside a conserved region\n");
  fprintf(ofp,"  -mm_spacer_indel [0.1]    indel inside a spacer\n");
  fprintf(ofp,"  -mm_switch_motif [0.05]    cost of switching to motif match\n");
  fprintf(ofp,"  -mm_switch_cons  [0.000001]  cost of switching to conserved match\n");

}

# line 135 "motifmatrix.dy"
MotifConsMatrix * new_MotifConsMatrix(TransFactorMatchSet * one,int starti,int endi,TransFactorMatchSet * two,int startj,int endj)
{
  int i;
  int j;
  int l;
  int k;
  int z;
  MotifConsMatrix * out;
  int motif_len;

  assert(one != NULL);
  assert(two != NULL);

  out = MotifConsMatrix_alloc_matrix(endi-starti,endj-startj);

  for(i=0;i< (endi-starti);i++) {
    for(j=0;j< (endj-startj) ;j++) {
      out->mat[i][j] = 0;
    }
  }

  for(l=0;l<one->len;l++) {
    if( one->match[l]->start < starti || one->match[l]->end > endi ) {
      continue;
    }
    for(k=0;k<two->len;k++) {
      if( two->match[k]->start < startj || two->match[k]->end > endj ) {
	continue;
      }
      
      if( two->match[k]->factor != one->match[l]->factor ) {
	continue;
      }

      motif_len = one->match[l]->end - one->match[l]->start;

      for(z=0;z<motif_len;z++) {
	out->mat[one->match[l]->start-starti+z][two->match[k]->start-startj+z] = 1;
      }      

    }
  }


  return out;
}
# line 145 "motifmatrix.c"
/* Function:  hard_link_MotifMatrixPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MotifMatrixPara *]
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixPara *]
 *
 */
MotifMatrixPara * hard_link_MotifMatrixPara(MotifMatrixPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MotifMatrixPara object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MotifMatrixPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixPara *]
 *
 */
MotifMatrixPara * MotifMatrixPara_alloc(void) 
{
    MotifMatrixPara * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MotifMatrixPara *) ckalloc (sizeof(MotifMatrixPara))) == NULL)  {  
      warn("MotifMatrixPara_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->comp_in_match = 0.0;    
    out->comp_out_match = 0.0;   
    out->comp_spacer = 0.0;  
    out->region_in = 0.0;    
    out->motif_indel = 0.0;  
    out->cons_indel = 0.0;   
    out->spacer_indel = 0.0; 
    out->spacer_to_cons = 0.0;   
    out->spacer_to_motif = 0.0;  
    out->spacer_duration = 0.0;  
    out->motif_duration = 0.0;   
    out->cons_duration = 0.0;    


    return out;  
}    


/* Function:  free_MotifMatrixPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MotifMatrixPara *]
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixPara *]
 *
 */
MotifMatrixPara * free_MotifMatrixPara(MotifMatrixPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MotifMatrixPara obj. Should be trappable");   
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


/* Function:  hard_link_MotifMatrixScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MotifMatrixScore *]
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixScore *]
 *
 */
MotifMatrixScore * hard_link_MotifMatrixScore(MotifMatrixScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MotifMatrixScore object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MotifMatrixScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixScore *]
 *
 */
MotifMatrixScore * MotifMatrixScore_alloc(void) 
{
    MotifMatrixScore * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MotifMatrixScore *) ckalloc (sizeof(MotifMatrixScore))) == NULL)    {  
      warn("MotifMatrixScore_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->comp_in_motif = NULL;   
    out->comp_out_motif = NULL;  
    out->comp_spacer = NULL; 
    out->region_in = 0;  
    out->motif_indel = 0;    
    out->cons_indel = 0; 
    out->spacer_indel = 0;   
    out->spacer_to_cons = 0; 
    out->spacer_to_motif = 0;    
    out->spacer_duration = 0;    
    out->motif_duration = 0; 
    out->cons_duration = 0;  


    return out;  
}    


/* Function:  free_MotifMatrixScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MotifMatrixScore *]
 *
 * Return [UNKN ]  Undocumented return value [MotifMatrixScore *]
 *
 */
MotifMatrixScore * free_MotifMatrixScore(MotifMatrixScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MotifMatrixScore obj. Should be trappable");  
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
    if( obj->comp_in_motif != NULL)  
      free_DnaMatrix(obj->comp_in_motif);    
    if( obj->comp_out_motif != NULL) 
      free_DnaMatrix(obj->comp_out_motif);   
    if( obj->comp_spacer != NULL)    
      free_DnaMatrix(obj->comp_spacer);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  MotifConsMatrix_alloc_matrix(leni,lenj)
 *
 * Descrip:    Allocates structure and matrix
 *
 *
 * Arg:        leni [UNKN ] Length of first dimension of matrix [int]
 * Arg:        lenj [UNKN ] Length of second dimension of matrix [int]
 *
 * Return [UNKN ]  Undocumented return value [MotifConsMatrix *]
 *
 */
MotifConsMatrix * MotifConsMatrix_alloc_matrix(int leni,int lenj) 
{
    MotifConsMatrix * out;  /* out is exported          */ 
    register int i; /* for stepping down matrix */ 
    register int j; /* for stepping across matrix */ 


    /* Call alloc function, return NULL if NULL */ 
    if((out = MotifConsMatrix_alloc()) == NULL)  
      return NULL;   


    /* Allocate memory for mat  */ 
    if((out->mat = (char **) ckcalloc (leni,sizeof(char *))) == NULL)    {  
      warn("Memory allocation problem in matrix for MotifConsMatrix mat, first pointer set");    
      ckfree(out);   
      return NULL;   
      }  


    /* Add NULL to all matrix pointers so free can be called */ 
    for(i=0;i<leni;i++)  
      out->mat[i] = NULL;    


    /* Allocate each matrix row */ 
    for(i=0;i<leni;i++)  {  
      out->mat[i] = (char *) ckcalloc (lenj,sizeof(char ));  
      if( out->mat[i] == NULL)   {  
        warn("Failed alloc on %d, calling free and returning NULL",i);   
        free_MotifConsMatrix(out);   
        return NULL;     
        }  
      }  


    for(i=0;i<leni;i++)  {  
      for(j=0;j<lenj;j++)    
        out->mat[i][j] = NULL;   
      }  


    out->leni=out->maxleni=leni;     
    out->lenj=out->maxlenj=lenj;     


    return out;  
}    


/* Function:  expand_MotifConsMatrix(obj,leni,lenj)
 *
 * Descrip:    Expands matrix. Rarely used
 *
 *
 * Arg:         obj [UNKN ] Undocumented argument [MotifConsMatrix *]
 * Arg:        leni [UNKN ] Undocumented argument [int]
 * Arg:        lenj [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_MotifConsMatrix(MotifConsMatrix * obj,int leni,int lenj) 
{
    int i;   
    int actualj;     


    if( obj == NULL)     {  
      warn("Trying to expand a MotifConsMatrix but is NULL!");   
      return FALSE;  
      }  


    if( leni <= obj->maxleni && lenj <= obj->maxlenj)    
      return TRUE;   


    if( obj->maxleni < leni )    {  
      if( (obj->mat=(char **) ckrealloc (obj->mat,sizeof(char *)*leni)) == NULL) 
        return FALSE;    
      obj->maxleni=obj->leni=leni;   
      }  
    if( lenj > obj->maxlenj )    
      actualj = lenj;    
    else actualj = obj->maxlenj; 
    for(i=0;i<obj->leni;i++) {  
      if((obj->mat[i] = (char *) realloc (obj->mat[i],sizeof(char ) * actualj)) == NULL) 
        return FALSE;    
      }  
    return TRUE;     
}    


/* Function:  hard_link_MotifConsMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MotifConsMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [MotifConsMatrix *]
 *
 */
MotifConsMatrix * hard_link_MotifConsMatrix(MotifConsMatrix * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MotifConsMatrix object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MotifConsMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MotifConsMatrix *]
 *
 */
MotifConsMatrix * MotifConsMatrix_alloc(void) 
{
    MotifConsMatrix * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MotifConsMatrix *) ckalloc (sizeof(MotifConsMatrix))) == NULL)  {  
      warn("MotifConsMatrix_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->mat = NULL; 
    out->leni=out->maxleni=0;    
    out->lenj=out->maxlenj=0;    


    return out;  
}    


/* Function:  free_MotifConsMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MotifConsMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [MotifConsMatrix *]
 *
 */
MotifConsMatrix * free_MotifConsMatrix(MotifConsMatrix * obj) 
{
    int return_early = 0;    
    int i;   
    int j;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MotifConsMatrix obj. Should be trappable");   
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
    if( obj->mat != NULL)    {  
      for(i=0;i<obj->leni;i++)   {  
        if( obj->mat[i] != NULL) 
          ckfree(obj->mat[i]);   
        }  
      ckfree(obj->mat);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

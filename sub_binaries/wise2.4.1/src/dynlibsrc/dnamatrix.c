#ifdef _cplusplus
extern "C" {
#endif
#include "dnamatrix.h"



/* Function:  new_CompMat_from_DnaMatrix_flat(dm)
 *
 * Descrip:    Builds a CompMat mapping of a DnaMatrix
 *
 *
 * Arg:        dm [UNKN ] Undocumented argument [DnaMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [CompMat *]
 *
 */
# line 44 "dnamatrix.dy"
CompMat * new_CompMat_from_DnaMatrix_flat(DnaMatrix * dm)
{
  int i,j;

  CompMat * mat;

  mat = CompMat_alloc();

  for(i=0;i<26;i++) {
    for(j=0;j<26;j++) {
      mat->comp[i][j] = NEGI;
    }
  }

  for(i=0;i<5;i++) {
    for(j=0;j<5;j++) {
	mat->comp[char_from_base(i)-'A'][char_from_base(j)-'A'] = dm->score[i][j];
    }
  }

  return mat;
}


/* Function:  DnaProbMatrix_from_match(match,nmask_type)
 *
 * Descrip:    Makes a probability matrix from simple match/mismatch 
 *             probabilities.
 *
 *
 *
 * Arg:             match [UNKN ] Undocumented argument [Probability]
 * Arg:        nmask_type [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProbMatrix *]
 *
 */
# line 73 "dnamatrix.dy"
DnaProbMatrix * DnaProbMatrix_from_match(Probability match,int nmask_type)
{
  int i,j;
  DnaProbMatrix * out;
  Probability factor;

  switch (nmask_type ) {
  case NMaskType_BASE :
    factor = ((1.0 - match)/4.0);
    break;
  case NMaskType_VARIABLE :
  case NMaskType_EXCLUDED :
  case NMaskType_BANNED :
    factor = ((1.0 - match)/3.0);
    break;
  default :
    warn("No valid mask type. Ugh!");
    return NULL;
  }

    

  out = DnaProbMatrix_alloc();
  
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      if( i == j ) {
	out->prob[i][j]  = match;
      } else {
	out->prob[i][j] = factor;
      }
    }
  }

  for(i=0;i<5;i++) {
    switch (nmask_type ) {
    case NMaskType_BASE :
      if( i == BASE_N ) {
	out->prob[i][i]  = match;
      } else {
	out->prob[BASE_N][i] = out->prob[i][BASE_N] = factor;
      }
      break;
    case NMaskType_VARIABLE :
      if( i == BASE_N ) {
	out->prob[i][i]  = 0.25;
      } else {
	out->prob[BASE_N][i] = out->prob[i][BASE_N] = 0.25;
      }
      break;

    case NMaskType_EXCLUDED :
      if( i == BASE_N ) {
	out->prob[i][i]  = 1.0;
      } else {
	out->prob[BASE_N][i] = out->prob[i][BASE_N] = 0.0;
      }
      break;
    case NMaskType_BANNED :
      out->prob[BASE_N][i] = out->prob[i][BASE_N] = 0.0;      
      break;
    default :
      warn("No valid mask type. Ugh! Shouldn't be here. A BAD  BAD bug!!!");
    }

  }

  return out;

}

/* Function:  flat_null_DnaProbMatrix(dpm)
 *
 * Descrip:    makes a odds of dpm via a 0.25 factor 
 *             into each base.
 *
 *
 * Arg:        dpm [UNKN ] Undocumented argument [DnaProbMatrix *]
 *
 */
# line 148 "dnamatrix.dy"
void flat_null_DnaProbMatrix(DnaProbMatrix * dpm)
{
  int i,j;

  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      dpm->prob[i][j] = dpm->prob[i][j]/0.25;
    }
  }

  return;
}
  
    
/* Function:  DnaMatrix_from_DnaProbMatrix(dpm)
 *
 * Descrip:    Maps probabilities to scores
 *
 *
 * Arg:        dpm [UNKN ] Undocumented argument [DnaProbMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatrix *]
 *
 */
# line 165 "dnamatrix.dy"
DnaMatrix * DnaMatrix_from_DnaProbMatrix(DnaProbMatrix * dpm)
{
  int i,j;
  DnaMatrix * out;
  
  out = DnaMatrix_alloc();
  for(i=0;i<5;i++) 
    for(j=0;j<5;j++) 
      out->score[i][j] = Probability2Score(dpm->prob[i][j]);

  return out;
}

/* Function:  show_DnaMatrix(dcm,ofp)
 *
 * Descrip:    Simple view of DnaMatrix
 *
 *
 * Arg:        dcm [UNKN ] Undocumented argument [DnaMatrix *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 181 "dnamatrix.dy"
void show_DnaMatrix(DnaMatrix * dcm,FILE * ofp)
{
  int i,j;

  for(i=0;i<5;i++) {
    for(j=0;j<5;j++) {
      fprintf(ofp,"%c %c - %d\n",char_from_base(i),char_from_base(j),dcm->score[i][j]);
    }
  }
}

/* Function:  show_DnaProbMatrix(dpm,ofp)
 *
 * Descrip:    Simple view of DnaProbMatrix
 *
 *
 * Arg:        dpm [UNKN ] Undocumented argument [DnaProbMatrix *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 195 "dnamatrix.dy"
void show_DnaProbMatrix(DnaProbMatrix * dpm,FILE * ofp)
{
  int i,j;

  for(i=0;i<5;i++) {
    for(j=0;j<5;j++) {
      fprintf(ofp,"%c %c - %g\n",char_from_base(i),char_from_base(j),dpm->prob[i][j]);
    }
  }
}


/* Function:  fail_safe_DnaMatrix_access(dm,one,two)
 *
 * Descrip:    Run-time checked that one and two are ok to pass
 *             into dm as bases
 *
 *
 *
 * Arg:         dm [UNKN ] DnaMatrix to get score from [DnaMatrix *]
 * Arg:        one [UNKN ] base of one sequence [base]
 * Arg:        two [UNKN ] base of the other sequence [base]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
# line 216 "dnamatrix.dy"
Score fail_safe_DnaMatrix_access(DnaMatrix * dm,base one,base two)
{
  if( dm == NULL) {
    warn("Passing in a NULL dna matrix into fail_safe_DnaMatrix_access, can't get a score therefore");
    return 0;
  }

  if( one < 0 || one > 4 || two < 0 || two > 4 ) {
    warn("In fail safe DnaMatrix access, trying to access a position [%d][%d] where this is meant to be 0-4",one,two);
    return 0;
  }

  return dm->score[one][two];
}

/* Function:  identity_DnaMatrix(id_score,mismatch)
 *
 * Descrip:    makes an idenity matrix wth id_score on the leading
 *             diagonal and mismatch elsewhere.
 *
 *
 *
 * Arg:        id_score [UNKN ] score of idenities [Score]
 * Arg:        mismatch [UNKN ] score of mistmatches [Score]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatrix *]
 *
 */
# line 239 "dnamatrix.dy"
DnaMatrix * identity_DnaMatrix(Score id_score,Score mismatch)
{
  DnaMatrix * out;
  int i;
  int j;

  out = DnaMatrix_alloc();

  for(i=0;i<4;i++)
    for(j=i;j<4;j++) {
      if( i == j ) 
	out->score[i][j] = id_score;
      else 
	out->score[i][j] = out->score[j][i] = mismatch;
    }
  
  for(i=0;i<4;i++) {
    out->score[i][4] = out->score[4][i] = 0;
  }
  out->score[4][4] = 0;

  return out;
}

# line 270 "dnamatrix.c"
/* Function:  hard_link_DnaProbMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProbMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProbMatrix *]
 *
 */
DnaProbMatrix * hard_link_DnaProbMatrix(DnaProbMatrix * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaProbMatrix object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaProbMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProbMatrix *]
 *
 */
DnaProbMatrix * DnaProbMatrix_alloc(void) 
{
    DnaProbMatrix * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProbMatrix *) ckalloc (sizeof(DnaProbMatrix))) == NULL)  {  
      warn("DnaProbMatrix_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* prob[5][5] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_DnaProbMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProbMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProbMatrix *]
 *
 */
DnaProbMatrix * free_DnaProbMatrix(DnaProbMatrix * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProbMatrix obj. Should be trappable"); 
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


/* Function:  hard_link_DnaMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatrix *]
 *
 */
DnaMatrix * hard_link_DnaMatrix(DnaMatrix * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaMatrix object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaMatrix *]
 *
 */
DnaMatrix * DnaMatrix_alloc(void) 
{
    DnaMatrix * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaMatrix *) ckalloc (sizeof(DnaMatrix))) == NULL)  {  
      warn("DnaMatrix_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* score[5][5] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_DnaMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DnaMatrix *]
 *
 */
DnaMatrix * free_DnaMatrix(DnaMatrix * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaMatrix obj. Should be trappable"); 
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

#ifdef _cplusplus
extern "C" {
#endif
#include "codonmatrix.h"



/* Function:  naive_CodonMatrixScore_from_prob(ct,cm)
 *
 * Descrip:    Combines CodonMatrixScore_from_CodonMatrix and naive_CodonMatrix
 *
 *
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        cm [UNKN ] Undocumented argument [CompProb *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
# line 23 "codonmatrix.dy"
CodonMatrixScore * naive_CodonMatrixScore_from_prob(CodonTable * ct,CompProb * cm)
{
  CodonMatrixScore * out;
  CodonMatrix * com;

  com = naive_CodonMatrix(ct,cm);

  out = CodonMatrixScore_from_CodonMatrix(com);

  free_CodonMatrix(com);

  return out;
}


/* Function:  CodonMatrixScore_from_CodonMatrix(cm)
 *
 * Descrip:    Makes a CodonMatrixScore from a CodonMatrix
 *
 *
 * Arg:        cm [UNKN ] Undocumented argument [CodonMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
# line 41 "codonmatrix.dy"
CodonMatrixScore * CodonMatrixScore_from_CodonMatrix(CodonMatrix * cm)
{
  int i,j;
  CodonMatrixScore * out;

  out = CodonMatrixScore_alloc();

  for(i=0;i<125;i++) 
    for(j=0;j<125;j++) 
      out->score[i][j] = Probability2Score(cm->prob[i][j]);
  
  return out;
}

/* Function:  naive_CodonMatrix(ct,comp)
 *
 * Descrip:    Builds a probability matrix
 *               No codon bias
 *               No errors
 *             N codons score 1.0, stop codons probability 0.00001
 *
 *
 * Arg:          ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        comp [UNKN ] Undocumented argument [CompProb *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrix *]
 *
 */
# line 61 "codonmatrix.dy"
CodonMatrix * naive_CodonMatrix(CodonTable * ct,CompProb * comp)
{
  int i;
  int j;
  CodonMatrix * out;


  out = CodonMatrix_alloc();

  for(i=0;i<125;i++)
    for(j=i;j<125;j++) {

      if( has_random_bases(i) == TRUE || has_random_bases(j) == TRUE ) {
	out->prob[i][j] = out->prob[j][i] = 1.0;
      } else if ( is_stop_codon(i,ct) == TRUE || is_stop_codon(j,ct) == TRUE ) {	
	out->prob[i][j] = out->prob[j][i] = 0.00001;
      } else { 
	out->prob[i][j] = out->prob[j][i] = comp->comp[aminoacid_no_from_codon(ct,i)][aminoacid_no_from_codon(ct,j)];
      }
    }
    

  return out;
}


/* Function:  naive_CodonMatrixScore(ct,comp)
 *
 * Descrip:    Builds a codon matrix from CompMat which assummes:
 *               No codon Bias
 *               No errors
 *             N codons score 0, stop codons score ??
 *
 *
 * Arg:          ct [UNKN ] CodonTable for codon->aa mapping [CodonTable *]
 * Arg:        comp [UNKN ] Comparison matrix for the score of the individual access [CompMat *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
# line 96 "codonmatrix.dy"
CodonMatrixScore * naive_CodonMatrixScore(CodonTable * ct,CompMat * comp)
{
  int i;
  int j;
  CodonMatrixScore * out;


  out = CodonMatrixScore_alloc();

  for(i=0;i<125;i++)
    for(j=i;j<125;j++) {
      if( has_random_bases(i) == TRUE || has_random_bases(j) == TRUE ) 
	out->score[i][j] = out->score[j][i] = 0;
      else 
	out->score[i][j] = out->score[j][i] = fail_safe_CompMat_access(comp,aminoacid_no_from_codon(ct,i),aminoacid_no_from_codon(ct,j));
    }

  return out;
}


/* Function:  show_CodonMatrixScore(cms,ct,ofp)
 *
 * Descrip:    Shows a codonmatrix
 *
 *
 * Arg:        cms [UNKN ] Undocumented argument [CodonMatrixScore *]
 * Arg:         ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 120 "codonmatrix.dy"
void show_CodonMatrixScore(CodonMatrixScore * cms,CodonTable * ct,FILE * ofp)
{
  int i;
  int j;

  for(i=0;i<125;i++) 
    for(j=0;j<125;j++) {
      fprintf(ofp,"%5d %c :%5d %c Score %5d\n",i,aminoacid_from_codon(ct,i),j,aminoacid_from_codon(ct,j),cms->score[i][j]);
    }

}


# line 150 "codonmatrix.c"
/* Function:  hard_link_CodonMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CodonMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrix *]
 *
 */
CodonMatrix * hard_link_CodonMatrix(CodonMatrix * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CodonMatrix object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CodonMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrix *]
 *
 */
CodonMatrix * CodonMatrix_alloc(void) 
{
    CodonMatrix * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CodonMatrix *) ckalloc (sizeof(CodonMatrix))) == NULL)  {  
      warn("CodonMatrix_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* prob[125][125] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_CodonMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CodonMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrix *]
 *
 */
CodonMatrix * free_CodonMatrix(CodonMatrix * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CodonMatrix obj. Should be trappable");   
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


/* Function:  hard_link_CodonMatrixScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CodonMatrixScore *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
CodonMatrixScore * hard_link_CodonMatrixScore(CodonMatrixScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CodonMatrixScore object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CodonMatrixScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
CodonMatrixScore * CodonMatrixScore_alloc(void) 
{
    CodonMatrixScore * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CodonMatrixScore *) ckalloc (sizeof(CodonMatrixScore))) == NULL)    {  
      warn("CodonMatrixScore_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* score[125][125] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_CodonMatrixScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CodonMatrixScore *]
 *
 * Return [UNKN ]  Undocumented return value [CodonMatrixScore *]
 *
 */
CodonMatrixScore * free_CodonMatrixScore(CodonMatrixScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CodonMatrixScore obj. Should be trappable");  
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

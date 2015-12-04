#ifdef _cplusplus
extern "C" {
#endif
#include "pwmdna.h"


/* Function:  max_prob_pwmDNA(in)
 *
 * Descrip:    maximum prob of pwmDNA
 *
 *
 * Arg:        in [UNKN ] Undocumented argument [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 64 "pwmdna.dy"
double max_prob_pwmDNA(pwmDNA * in)
{
  int i,k;
  double temp_max;
  double score = 0.0;

  for(i=0;i<in->len;i++) {
    temp_max = in->pos[i]->emit[0];
    for(k=1;k<4;k++) {
      if( temp_max < in->pos[i]->emit[k] ) {
	temp_max = in->pos[i]->emit[k];
      }
    }
    score *= temp_max;
  }

  return score;
}


/* Function:  min_prob_pwmDNA(in)
 *
 * Descrip:    minimum prob of pwmDNA
 *
 *
 * Arg:        in [UNKN ] Undocumented argument [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [double]
 *
 */
# line 87 "pwmdna.dy"
double min_prob_pwmDNA(pwmDNA * in)
{
  int i,k;
  double temp_min;
  double score = 0.0;

  for(i=0;i<in->len;i++) {
    temp_min = in->pos[i]->emit[0];
    for(k=1;k<4;k++) {
      if( temp_min > in->pos[i]->emit[k] ) {
	temp_min = in->pos[i]->emit[k];
      }
    }
    score *= temp_min;
  }

  return score;
}

/* Function:  circular_permuted_pwmDNA(in,rotate_number)
 *
 * Descrip:    Provides a rotated pwm for randomisation
 *
 *
 * Arg:                   in [UNKN ] Undocumented argument [pwmDNA *]
 * Arg:        rotate_number [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
# line 109 "pwmdna.dy"
pwmDNA * circular_permuted_pwmDNA(pwmDNA * in,int rotate_number)
{
  int i;
  pwmDNA * out;
  pwmColProb * col;
  int j;
  int pos;

  assert(in != NULL);

  out = pwmDNA_alloc_len(in->len);
  for(i=0;i<in->len;i++) {
    pos = (i+rotate_number) % in->len;
    col = pwmColProb_alloc();
    add_pwmDNA(out,col);
    for(j=0;j<5;j++) {
      col->emit[j] = in->pos[pos]->emit[j];
    }
  }

  return out;

}

/* Function:  score_pwmDNAScore_Sequence(pds,s,pos)
 *
 * Descrip:    This gives back a Score from a particular sequence and
 *             position
 *
 *
 * Arg:        pds [UNKN ] Undocumented argument [pwmDNAScore *]
 * Arg:          s [UNKN ] Undocumented argument [Sequence *]
 * Arg:        pos [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
# line 137 "pwmdna.dy"
Score score_pwmDNAScore_Sequence(pwmDNAScore * pds,Sequence * s,int pos)
{
  int score;

  if( pds->len + pos > s->len ) {
    warn("For sequence %s, position %d is unable to be matched to pwmDNA of length %d",s->name,pos,pds->len);
  }

  score = score_pwmDNAScore_string(pds,s->seq+pos);
  fprintf(stderr,"Making score %d for %d\n",score,pos);
  return score;
}

/* Function:  score_pwmDNAScore_string(pds,str)
 *
 * Descrip:    This gives back a Score from a particular string
 *
 *
 * Arg:        pds [UNKN ] Undocumented argument [pwmDNAScore *]
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
# line 153 "pwmdna.dy"
Score score_pwmDNAScore_string(pwmDNAScore * pds,char * str)
{
  int i;
  Score total = 0.0;

  /*
  if( strlen(str) < pds->len ) {
    warn("String [%s] is shorter than the length of the pds. Should not be using it!",str);
    return NEGI;
  }
  */

  for(i=0;i<pds->len;i++)
    total += pds->pos[i]->emit[base_from_char(str[i])];

/*  printf("Score %c%c%c to score with %s  [%d]\n",str[0],str[1],str[2],pds->ref_align->seq[0]->seq,total);*/
 
  return total;
}

/* Function:  prob_pwmDNA_string(pds,str)
 *
 * Descrip:    This gives back a Probability from a particular string
 *
 *
 * Arg:        pds [UNKN ] Undocumented argument [pwmDNA *]
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [Probability]
 *
 */
# line 176 "pwmdna.dy"
Probability prob_pwmDNA_string(pwmDNA * pds,char * str)
{
  int i;
  Probability total = 1.0;


  for(i=0;i<pds->len;i++) {
    if( str[i] == '\0' ) {
      return 0.0;
    } else {
      total *= pds->pos[i]->emit[base_from_char(str[i])];
      /*      fprintf(stdout,"As %d, total %.6f with %c\n",i,total,str[i]);*/
    }
  }

  return total;
}

/* Function:  fold_randommodel_pwmDNA(pd,rmd)
 *
 * Descrip:    This folds in a randommodel into pwmDNA
 *
 *
 * Arg:         pd [UNKN ] Undocumented argument [pwmDNA *]
 * Arg:        rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 */
# line 197 "pwmdna.dy"
void fold_randommodel_pwmDNA(pwmDNA * pd,RandomModelDNA * rmd)
{
  int i;
  int j;

  assert(pd);
  assert(rmd);
  for(i=0;i<pd->len;i++)
    for(j=0;j<5;j++)
      if( rmd->base[j] < 0.00000000000000001 ) {
	warn("Zero base %d, skipping",j);
      } else {
	pd->pos[i]->emit[j] /= rmd->base[j];
      }

}

/* Function:  pwmDNA_from_SeqAlign(sa,simple_pseudocount)
 *
 * Descrip:    This function makes a single pwmDNA from
 *             a SeqAlign
 *
 *             FIXME: This DOES NOT handle ambiguity codes well
 *
 *
 * Arg:                        sa [UNKN ] Undocumented argument [SeqAlign *]
 * Arg:        simple_pseudocount [UNKN ] Undocumented argument [double]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
# line 220 "pwmdna.dy"
pwmDNA * pwmDNA_from_SeqAlign(SeqAlign * sa,double simple_pseudocount)
{
  pwmDNA * out;
  pwmColProb * col;
  ColumnCount * cc;
  int i;

  assert(sa);

  out = pwmDNA_alloc_std();
  out->ref_align = hard_link_SeqAlign(sa);

  for(i=0;i<sa->seq[0]->len;i++) {
    cc = ColumnCount_from_SeqAlign(sa,i);
    col = pwmColProb_from_ColumnCount(cc,simple_pseudocount);
    add_pwmDNA(out,col);
    free_ColumnCount(cc);
  }

  return out;
}


/* Function:  show_pwmDNA_col(pd,ofp)
 *
 * Descrip:    Shows a columns along the page
 *
 *
 * Arg:         pd [UNKN ] Undocumented argument [pwmDNA *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 246 "pwmdna.dy"
void show_pwmDNA_col(pwmDNA * pd,FILE * ofp)
{
  char base[] = "ATGC";
  int i;
  int j;

  for(j=0;j<4;j++) {
    fprintf(ofp," %c ",base[j]);
    for(i=0;i<pd->len;i++) {
      fprintf(ofp," %.4f ",pd->pos[i]->emit[base_from_char(base[j])]);
    }
    fprintf(ofp,"\n");
  }
}

  

/* Function:  pwmColProb_from_ColumnCount(cc,simple_pseudocount)
 *
 * Descrip:    This function makes a single pwmColProb from
 *             a ColumnCount, applying a simple pseudocount
 *             method
 *
 *             FIXME: This DOES NOT handle ambiguity codes well
 *
 *
 * Arg:                        cc [UNKN ] Undocumented argument [ColumnCount *]
 * Arg:        simple_pseudocount [UNKN ] Undocumented argument [double]
 *
 * Return [UNKN ]  Undocumented return value [pwmColProb *]
 *
 */
# line 271 "pwmdna.dy"
pwmColProb * pwmColProb_from_ColumnCount(ColumnCount * cc,double simple_pseudocount)
{
  pwmColProb * out;
  double total;
  char base[] = "ATGC";
  int i;
  out = pwmColProb_alloc();
  
  for(i=0,total = 0.0;i<4;i++)
    total += (cc->count[base[i]-'A'] + simple_pseudocount);
  
  for(i=0;i<4;i++)
    out->emit[base_from_char(base[i])] = (cc->count[base[i]-'A'] + simple_pseudocount) / total;
  out->emit[4] = 1.0;

  return out;
}
  
/* Function:  pwmDNAScore_from_pwmDNA_RandomModelDNA(pwm,rmd)
 *
 * Descrip:    This function makes score represention of a
 *             position weight matrix from a probability, 
 *             with a random model folded in on-the-fly
 *
 *
 * Arg:        pwm [UNKN ] Undocumented argument [pwmDNA *]
 * Arg:        rmd [UNKN ] Undocumented argument [RandomModelDNA *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
# line 294 "pwmdna.dy"
pwmDNAScore * pwmDNAScore_from_pwmDNA_RandomModelDNA(pwmDNA * pwm,RandomModelDNA * rmd)
{
  int i;
  pwmDNAScore * out;

  out = pwmDNAScore_alloc_std();
  out->ref_align = hard_link_SeqAlign(pwm->ref_align);
  for(i=0;i<pwm->len;i++) {
    add_pwmDNAScore(out,pwmColScore_from_pwmColProb_rmd(pwm->pos[i],rmd));
  }

  return out;
}

/* Function:  pwmDNAScore_from_pwmDNA(pwm)
 *
 * Descrip:    This function makes score represention of a
 *             position weight matrix from a probability
 *
 *
 * Arg:        pwm [UNKN ] Undocumented argument [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
# line 312 "pwmdna.dy"
pwmDNAScore * pwmDNAScore_from_pwmDNA(pwmDNA * pwm)
{
  int i;
  pwmDNAScore * out;

  out = pwmDNAScore_alloc_std();
  out->ref_align = hard_link_SeqAlign(pwm->ref_align);

  for(i=0;i<pwm->len;i++) {
    add_pwmDNAScore(out,pwmColScore_from_pwmColProb(pwm->pos[i]));
  }

  return out;
}

/* Function:  pwmColScore_from_pwmColProb(p)
 *
 * Descrip:    This function makes a score representation of a
 *             position weight matrix from a column representation
 *
 *
 * Arg:        p [UNKN ] Undocumented argument [pwmColProb *]
 *
 * Return [UNKN ]  Undocumented return value [pwmColScore *]
 *
 */
# line 332 "pwmdna.dy"
pwmColScore * pwmColScore_from_pwmColProb(pwmColProb * p)
{
  pwmColScore * out;

  out = pwmColScore_alloc();
  Probability2Score_move(p->emit,out->emit,5);

  return out;
}

/* Function:  pwmColScore_from_pwmColProb_rmd(p,rmd)
 *
 * Descrip:    This function makes a score representation of a
 *             position weight matrix from a column representation,
 *             with a RandomModel factored in
 *
 *
 * Arg:          p [UNKN ] Undocumented argument [pwmColProb *]
 * Arg:        rmd [UNKN ] Undocumented argument [RandomModelDNA  *]
 *
 * Return [UNKN ]  Undocumented return value [pwmColScore *]
 *
 */
# line 348 "pwmdna.dy"
pwmColScore * pwmColScore_from_pwmColProb_rmd(pwmColProb * p,RandomModelDNA  * rmd)
{
  pwmColScore * out;
  int i;

  out = pwmColScore_alloc();
  for(i=0;i<5;i++)
    out->emit[i] = Probability2Score(p->emit[i] / rmd->base[i]);

  return out;
}

# line 393 "pwmdna.c"
/* Function:  hard_link_pwmColScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [pwmColScore *]
 *
 * Return [UNKN ]  Undocumented return value [pwmColScore *]
 *
 */
pwmColScore * hard_link_pwmColScore(pwmColScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a pwmColScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  pwmColScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmColScore *]
 *
 */
pwmColScore * pwmColScore_alloc(void) 
{
    pwmColScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(pwmColScore *) ckalloc (sizeof(pwmColScore))) == NULL)  {  
      warn("pwmColScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* emit[5] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_pwmColScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [pwmColScore *]
 *
 * Return [UNKN ]  Undocumented return value [pwmColScore *]
 *
 */
pwmColScore * free_pwmColScore(pwmColScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a pwmColScore obj. Should be trappable");   
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


/* Function:  hard_link_pwmColProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [pwmColProb *]
 *
 * Return [UNKN ]  Undocumented return value [pwmColProb *]
 *
 */
pwmColProb * hard_link_pwmColProb(pwmColProb * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a pwmColProb object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  pwmColProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmColProb *]
 *
 */
pwmColProb * pwmColProb_alloc(void) 
{
    pwmColProb * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(pwmColProb *) ckalloc (sizeof(pwmColProb))) == NULL)    {  
      warn("pwmColProb_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* emit[5] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_pwmColProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [pwmColProb *]
 *
 * Return [UNKN ]  Undocumented return value [pwmColProb *]
 *
 */
pwmColProb * free_pwmColProb(pwmColProb * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a pwmColProb obj. Should be trappable");    
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


/* Function:  swap_pwmDNA(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_pwmDNA
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [pwmColProb **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_pwmDNA(pwmColProb ** list,int i,int j)  
{
    pwmColProb * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_pwmDNA(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_pwmDNA which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [pwmColProb **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_pwmDNA(pwmColProb ** list,int left,int right,int (*comp)(pwmColProb * ,pwmColProb * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_pwmDNA(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_pwmDNA (list,++last,i); 
      }  
    swap_pwmDNA (list,left,last);    
    qsort_pwmDNA(list,left,last-1,comp); 
    qsort_pwmDNA(list,last+1,right,comp);    
}    


/* Function:  sort_pwmDNA(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_pwmDNA
 *
 *
 * Arg:         obj [UNKN ] Object containing list [pwmDNA *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_pwmDNA(pwmDNA * obj,int (*comp)(pwmColProb *, pwmColProb *)) 
{
    qsort_pwmDNA(obj->pos,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_pwmDNA(obj,len)
 *
 * Descrip:    Really an internal function for add_pwmDNA
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [pwmDNA *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_pwmDNA(pwmDNA * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_pwmDNA called with no need"); 
      return TRUE;   
      }  


    if( (obj->pos = (pwmColProb ** ) ckrealloc (obj->pos,sizeof(pwmColProb *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_pwmDNA, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_pwmDNA(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [pwmDNA *]
 * Arg:        add [OWNER] Object to add to the list [pwmColProb *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_pwmDNA(pwmDNA * obj,pwmColProb * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_pwmDNA(obj,obj->len + pwmDNALISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->pos[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_pwmDNA(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_pwmDNA(pwmDNA * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->pos[i] != NULL)   {  
        free_pwmColProb(obj->pos[i]);    
        obj->pos[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  pwmDNA_alloc_std(void)
 *
 * Descrip:    Equivalent to pwmDNA_alloc_len(pwmDNALISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * pwmDNA_alloc_std(void) 
{
    return pwmDNA_alloc_len(pwmDNALISTLENGTH);   
}    


/* Function:  pwmDNA_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * pwmDNA_alloc_len(int len) 
{
    pwmDNA * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = pwmDNA_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->pos = (pwmColProb ** ) ckcalloc (len,sizeof(pwmColProb *))) == NULL)    {  
      warn("Warning, ckcalloc failed in pwmDNA_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_pwmDNA(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * hard_link_pwmDNA(pwmDNA * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a pwmDNA object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  pwmDNA_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * pwmDNA_alloc(void) 
{
    pwmDNA * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(pwmDNA *) ckalloc (sizeof(pwmDNA))) == NULL)    {  
      warn("pwmDNA_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->pos = NULL; 
    out->len = out->maxlen = 0;  
    out->ref_align = NULL;   


    return out;  
}    


/* Function:  free_pwmDNA(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [pwmDNA *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNA *]
 *
 */
pwmDNA * free_pwmDNA(pwmDNA * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a pwmDNA obj. Should be trappable");    
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
    if( obj->pos != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->pos[i] != NULL) 
          free_pwmColProb(obj->pos[i]);  
        }  
      ckfree(obj->pos);  
      }  
    if( obj->ref_align != NULL)  
      free_SeqAlign(obj->ref_align);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_pwmDNAScore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_pwmDNAScore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [pwmColScore **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_pwmDNAScore(pwmColScore ** list,int i,int j)  
{
    pwmColScore * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_pwmDNAScore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_pwmDNAScore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [pwmColScore **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_pwmDNAScore(pwmColScore ** list,int left,int right,int (*comp)(pwmColScore * ,pwmColScore * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_pwmDNAScore(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_pwmDNAScore (list,++last,i);    
      }  
    swap_pwmDNAScore (list,left,last);   
    qsort_pwmDNAScore(list,left,last-1,comp);    
    qsort_pwmDNAScore(list,last+1,right,comp);   
}    


/* Function:  sort_pwmDNAScore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_pwmDNAScore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [pwmDNAScore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_pwmDNAScore(pwmDNAScore * obj,int (*comp)(pwmColScore *, pwmColScore *)) 
{
    qsort_pwmDNAScore(obj->pos,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_pwmDNAScore(obj,len)
 *
 * Descrip:    Really an internal function for add_pwmDNAScore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [pwmDNAScore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_pwmDNAScore(pwmDNAScore * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_pwmDNAScore called with no need");    
      return TRUE;   
      }  


    if( (obj->pos = (pwmColScore ** ) ckrealloc (obj->pos,sizeof(pwmColScore *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_pwmDNAScore, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_pwmDNAScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [pwmDNAScore *]
 * Arg:        add [OWNER] Object to add to the list [pwmColScore *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_pwmDNAScore(pwmDNAScore * obj,pwmColScore * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_pwmDNAScore(obj,obj->len + pwmDNAScoreLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->pos[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_pwmDNAScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [pwmDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_pwmDNAScore(pwmDNAScore * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->pos[i] != NULL)   {  
        free_pwmColScore(obj->pos[i]);   
        obj->pos[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  pwmDNAScore_alloc_std(void)
 *
 * Descrip:    Equivalent to pwmDNAScore_alloc_len(pwmDNAScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * pwmDNAScore_alloc_std(void) 
{
    return pwmDNAScore_alloc_len(pwmDNAScoreLISTLENGTH); 
}    


/* Function:  pwmDNAScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * pwmDNAScore_alloc_len(int len) 
{
    pwmDNAScore * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = pwmDNAScore_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->pos = (pwmColScore ** ) ckcalloc (len,sizeof(pwmColScore *))) == NULL)  {  
      warn("Warning, ckcalloc failed in pwmDNAScore_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_pwmDNAScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [pwmDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * hard_link_pwmDNAScore(pwmDNAScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a pwmDNAScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  pwmDNAScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * pwmDNAScore_alloc(void) 
{
    pwmDNAScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(pwmDNAScore *) ckalloc (sizeof(pwmDNAScore))) == NULL)  {  
      warn("pwmDNAScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->pos = NULL; 
    out->len = out->maxlen = 0;  
    out->ref_align = NULL;   


    return out;  
}    


/* Function:  free_pwmDNAScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [pwmDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [pwmDNAScore *]
 *
 */
pwmDNAScore * free_pwmDNAScore(pwmDNAScore * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a pwmDNAScore obj. Should be trappable");   
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
    if( obj->pos != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->pos[i] != NULL) 
          free_pwmColScore(obj->pos[i]); 
        }  
      ckfree(obj->pos);  
      }  
    if( obj->ref_align != NULL)  
      free_SeqAlign(obj->ref_align);     


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

#ifdef _cplusplus
extern "C" {
#endif
#include "dnaprofile.h"



# line 46 "dnaprofile.dy"
TransFactor * new_TransFactor_from_DnaProfile(DnaProfile * dp,char * name)
{
  TransFactor * out;
  assert(dp != NULL);
  
  out = TransFactor_alloc();
  out->name = stringalloc(name);
  out->seed = hard_link_SeqAlign(dp->sa);
  
  return out;
}

# line 58 "dnaprofile.dy"
double bits_entropy_DnaProfile(DnaProfile * dnap,RandomModelDNA * background)
{
  double entropy = 0.0;
  double prob = 0;
  int i;
  int j;

  
  for(i=0;i<dnap->len;i++) {
    for(j=0;j<4;j++) {
      prob = dnap->col[i]->emit[j] * background->base[j];
      entropy +=  prob * (-(log(prob)/log(2.0)));
    }
  }
  
  return entropy;
}

# line 76 "dnaprofile.dy"
SeqAlign * merged_SeqAlign(DnaProfile * one,DnaProfile * two,AlnBlock * alb)
{
  SeqAlign * out;
  Sequence * temp;
  int i;
  int j;
  int len = 0;
  AlnColumn * alc;
  char * temp_seq;
  out = SeqAlign_alloc_std();

  for(alc=alb->start,len = 0;alc != NULL;alc = alc->next ) {
    len++;
  }

  temp_seq = calloc(len+1,sizeof(char));
  /* loop along the top sequence */
  
  for(i=0;i<one->sa->len;i++) {
    j = 0;
    for(alc=alb->start;alc!=NULL;alc=alc->next) {
      if( strcmp(alc->alu[1]->text_label,"END") == 0 ) {
	break;
      }
      if( strcmp(alc->alu[0]->text_label,"MATCH") == 0 ||
	  strcmp(alc->alu[0]->text_label,"INSERT") == 0 ) {
	temp_seq[j++] = one->sa->seq[i]->seq[alc->alu[0]->end];
      } else if( strcmp(alc->alu[0]->text_label,"SEQ_UNMATCHED") == 0) {
	temp_seq[j++] = tolower(one->sa->seq[i]->seq[alc->alu[0]->end]);
      } else {
	temp_seq[j++] = '-';
      }
    }
    temp_seq[j] = '\0';
    temp = Sequence_from_static_memory(one->sa->seq[i]->name,temp_seq);
    add_SeqAlign(out,temp);
  }

  for(i=0;i<two->sa->len;i++) {
    j = 0;
    for(alc=alb->start;alc!=NULL;alc=alc->next) {
      if( strcmp(alc->alu[1]->text_label,"END") == 0 ) {
	break;
      }
      if( strcmp(alc->alu[1]->text_label,"MATCH") == 0 ||
	  strcmp(alc->alu[1]->text_label,"INSERT") == 0 ) {
	temp_seq[j++] = two->sa->seq[i]->seq[alc->alu[1]->end];
      } else if( strcmp(alc->alu[1]->text_label,"SEQ_UNMATCHED") == 0) {
	temp_seq[j++] = tolower(two->sa->seq[i]->seq[alc->alu[1]->end]);
      } else {
	temp_seq[j++] = '-';
      }
    }
    temp_seq[j] = '\0';
    temp = Sequence_from_static_memory(two->sa->seq[i]->name,temp_seq);
    add_SeqAlign(out,temp);
  }
      
 
  return out;
}

# line 138 "dnaprofile.dy"
DnaProfileScore * DnaProfileScore_from_DnaProfile(DnaProfile * dp)
{
  int i,k;
  DnaProfileColScore * col;
  DnaProfileScore * out;

  out = DnaProfileScore_alloc_len(dp->len);


  for(i=0;i<dp->len;i++) {
    col = DnaProfileColScore_alloc();
    Probability2Score_move(dp->col[i]->emit,col->emit,4);

    Probability2Score_move(dp->col[i]->trans,col->trans,DnaProfile_TRANS_LENGTH);

    add_DnaProfileScore(out,col);
  }

  return out;
}

# line 159 "dnaprofile.dy"
void fold_RandomModel_DnaProfile(DnaProfile * dp,RandomModelDNA * rm)
{
  int i,k;
  for(i=0;i<dp->len;i++) {
    for(k=0;k<4;k++) {
      dp->col[i]->emit[k] /= rm->base[k];
    }
  }

  dp->folded_random = TRUE;
}

# line 171 "dnaprofile.dy"
DnaProfile * naive_DnaProfile_from_Sequence(Sequence * seq,double id,double m2i,double m2d,double i2i,double d2d)
{
  DnaProfile * out;
  DnaProfileCol * col;
  int i,k;
  SeqAlign * sa;

  sa = SeqAlign_alloc_len(1);
  add_SeqAlign(sa,hard_link_Sequence(seq));

  
  out= DnaProfile_alloc_len(seq->len);
  out->sa = sa;

  for(i=0;i<seq->len;i++) {
    col = DnaProfileCol_alloc();

    for(k=0;k<4;k++) {
      col->emit[k]= (1.0-id)/3.0;
    }

    if( base_from_char(seq->seq[i]) != BASE_N ) {
      col->emit[base_from_char(seq->seq[i])] = id;
    }


    col->trans[DnaProfile_M2I] = m2i;
    col->trans[DnaProfile_M2D] = m2d;
    col->trans[DnaProfile_I2I] = i2i;
    col->trans[DnaProfile_D2D] = d2d;
    col->trans[DnaProfile_M2M] = (1.0 - m2i - m2d);
    col->trans[DnaProfile_D2M] = (1.0 - d2d);
    col->trans[DnaProfile_I2M] = (1.0 - i2i);
    
    add_DnaProfile(out,col);
  }


  return out;
}


# line 213 "dnaprofile.dy"
DnaProfile * naive_DnaProfile_from_SeqAlign(SeqAlign * sa,double simple_pseudo,double m2i,double m2d,double i2i,double d2d)
{
  int i;
  int j;
  int k;
  DnaProfile * out;
  DnaProfileCol * col;
  int base_n;
  double count[4];
  double total;

  assert(sa != NULL);
  assert(sa->len > 1);
  assert(sa->seq[0]->len > 1);

  out = DnaProfile_alloc_len(sa->seq[0]->len);
  out->sa = hard_link_SeqAlign(sa);

  total = (4*simple_pseudo)+sa->len;

  for(i=0;i<sa->seq[0]->len;i++) {

    for(k=0;k<4;k++) {
      count[k] = simple_pseudo;
    }
    for(j=0;j<sa->len;j++) {
      assert(sa->seq[j]->len > i);
      base_n = base_from_char(sa->seq[j]->seq[i]);
      if( base_n < 4 ) {
	count[base_n]++;
      } else {
	for(k=0;k<4;k++) {
	  count[base_n] += 0.25;
	}
      }
    }

    col = DnaProfileCol_alloc();
    col->seqalign_col = i;

    for(k=0;k<4;k++) {
      col->emit[k] = count[k] / total;
    }

    col->trans[DnaProfile_M2I] = m2i;
    col->trans[DnaProfile_M2D] = m2d;
    col->trans[DnaProfile_I2I] = i2i;
    col->trans[DnaProfile_D2D] = d2d;
    col->trans[DnaProfile_M2M] = (1.0 - m2i - m2d);
    col->trans[DnaProfile_D2M] = (1.0 - d2d);
    col->trans[DnaProfile_I2M] = (1.0 - i2i);


 
    add_DnaProfile(out,col);

    
	
  }

  return out;
}


# line 277 "dnaprofile.dy"
void show_DnaProfile(DnaProfile * dnap,RandomModelDNA * rm,FILE * ofp)
{
  int i;
  int k;
  double max;
  double max_base;
  double ent;

  ent = bits_entropy_DnaProfile(dnap,rm);
  
  fprintf(ofp,"Entropy %.4f (%.4f per column)\n",ent,ent/dnap->len);

  for(i=0;i<dnap->len;i++) {
    max = dnap->col[i]->emit[0];
    max_base =0;

    for(k=1;k<4;k++) {
      if( dnap->col[i]->emit[k] > max ) {
	max = dnap->col[i]->emit[k];
	max_base =k;
      }
    }
    
    fprintf(ofp,"%5d %c ",i,char_from_base(max_base));
    for(k=0;k<4;k++) {
      fprintf(ofp," %c:%.3f ",char_from_base(k),dnap->col[i]->emit[k]);
    }
      
    fprintf(ofp,"\n");
  }

  write_selex_SeqAlign(dnap->sa,15,50,ofp);

  fprintf(ofp,"//\n");

}

# line 283 "dnaprofile.c"
/* Function:  hard_link_DnaProfileCol(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileCol *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileCol *]
 *
 */
DnaProfileCol * hard_link_DnaProfileCol(DnaProfileCol * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaProfileCol object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaProfileCol_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileCol *]
 *
 */
DnaProfileCol * DnaProfileCol_alloc(void) 
{
    DnaProfileCol * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProfileCol *) ckalloc (sizeof(DnaProfileCol))) == NULL)  {  
      warn("DnaProfileCol_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* emit[5] is an array: no default possible */ 
    /* trans[DnaProfile_TRANS_LENGTH] is an array: no default possible */ 
    out->seqalign_col = 0;   


    return out;  
}    


/* Function:  free_DnaProfileCol(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileCol *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileCol *]
 *
 */
DnaProfileCol * free_DnaProfileCol(DnaProfileCol * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProfileCol obj. Should be trappable"); 
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


/* Function:  hard_link_DnaProfileColScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileColScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileColScore *]
 *
 */
DnaProfileColScore * hard_link_DnaProfileColScore(DnaProfileColScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaProfileColScore object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaProfileColScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileColScore *]
 *
 */
DnaProfileColScore * DnaProfileColScore_alloc(void) 
{
    DnaProfileColScore * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProfileColScore *) ckalloc (sizeof(DnaProfileColScore))) == NULL)    {  
      warn("DnaProfileColScore_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* emit[5] is an array: no default possible */ 
    /* trans[DnaProfile_TRANS_LENGTH] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_DnaProfileColScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileColScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileColScore *]
 *
 */
DnaProfileColScore * free_DnaProfileColScore(DnaProfileColScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProfileColScore obj. Should be trappable");    
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


/* Function:  swap_DnaProfileScore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_DnaProfileScore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DnaProfileColScore **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_DnaProfileScore(DnaProfileColScore ** list,int i,int j)  
{
    DnaProfileColScore * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_DnaProfileScore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_DnaProfileScore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DnaProfileColScore **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_DnaProfileScore(DnaProfileColScore ** list,int left,int right,int (*comp)(DnaProfileColScore * ,DnaProfileColScore * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_DnaProfileScore(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_DnaProfileScore (list,++last,i);    
      }  
    swap_DnaProfileScore (list,left,last);   
    qsort_DnaProfileScore(list,left,last-1,comp);    
    qsort_DnaProfileScore(list,last+1,right,comp);   
}    


/* Function:  sort_DnaProfileScore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_DnaProfileScore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DnaProfileScore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_DnaProfileScore(DnaProfileScore * obj,int (*comp)(DnaProfileColScore *, DnaProfileColScore *)) 
{
    qsort_DnaProfileScore(obj->col,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_DnaProfileScore(obj,len)
 *
 * Descrip:    Really an internal function for add_DnaProfileScore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfileScore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DnaProfileScore(DnaProfileScore * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_DnaProfileScore called with no need");    
      return TRUE;   
      }  


    if( (obj->col = (DnaProfileColScore ** ) ckrealloc (obj->col,sizeof(DnaProfileColScore *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_DnaProfileScore, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_DnaProfileScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfileScore *]
 * Arg:        add [OWNER] Object to add to the list [DnaProfileColScore *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_DnaProfileScore(DnaProfileScore * obj,DnaProfileColScore * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_DnaProfileScore(obj,obj->len + DnaProfileScoreLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->col[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_DnaProfileScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaProfileScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DnaProfileScore(DnaProfileScore * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->col[i] != NULL)   {  
        free_DnaProfileColScore(obj->col[i]);    
        obj->col[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  DnaProfileScore_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaProfileScore_alloc_len(DnaProfileScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileScore *]
 *
 */
DnaProfileScore * DnaProfileScore_alloc_std(void) 
{
    return DnaProfileScore_alloc_len(DnaProfileScoreLISTLENGTH); 
}    


/* Function:  DnaProfileScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileScore *]
 *
 */
DnaProfileScore * DnaProfileScore_alloc_len(int len) 
{
    DnaProfileScore * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = DnaProfileScore_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->col = (DnaProfileColScore ** ) ckcalloc (len,sizeof(DnaProfileColScore *))) == NULL)    {  
      warn("Warning, ckcalloc failed in DnaProfileScore_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_DnaProfileScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileScore *]
 *
 */
DnaProfileScore * hard_link_DnaProfileScore(DnaProfileScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaProfileScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaProfileScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileScore *]
 *
 */
DnaProfileScore * DnaProfileScore_alloc(void) 
{
    DnaProfileScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProfileScore *) ckalloc (sizeof(DnaProfileScore))) == NULL)  {  
      warn("DnaProfileScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->col = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_DnaProfileScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileScore *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileScore *]
 *
 */
DnaProfileScore * free_DnaProfileScore(DnaProfileScore * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProfileScore obj. Should be trappable");   
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
    if( obj->col != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->col[i] != NULL) 
          free_DnaProfileColScore(obj->col[i]);  
        }  
      ckfree(obj->col);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_DnaProfile(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_DnaProfile
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DnaProfileCol **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_DnaProfile(DnaProfileCol ** list,int i,int j)  
{
    DnaProfileCol * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_DnaProfile(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_DnaProfile which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DnaProfileCol **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_DnaProfile(DnaProfileCol ** list,int left,int right,int (*comp)(DnaProfileCol * ,DnaProfileCol * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_DnaProfile(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_DnaProfile (list,++last,i); 
      }  
    swap_DnaProfile (list,left,last);    
    qsort_DnaProfile(list,left,last-1,comp); 
    qsort_DnaProfile(list,last+1,right,comp);    
}    


/* Function:  sort_DnaProfile(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_DnaProfile
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DnaProfile *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_DnaProfile(DnaProfile * obj,int (*comp)(DnaProfileCol *, DnaProfileCol *)) 
{
    qsort_DnaProfile(obj->col,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_DnaProfile(obj,len)
 *
 * Descrip:    Really an internal function for add_DnaProfile
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfile *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DnaProfile(DnaProfile * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_DnaProfile called with no need"); 
      return TRUE;   
      }  


    if( (obj->col = (DnaProfileCol ** ) ckrealloc (obj->col,sizeof(DnaProfileCol *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_DnaProfile, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_DnaProfile(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfile *]
 * Arg:        add [OWNER] Object to add to the list [DnaProfileCol *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_DnaProfile(DnaProfile * obj,DnaProfileCol * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_DnaProfile(obj,obj->len + DnaProfileLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->col[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_DnaProfile(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaProfile *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DnaProfile(DnaProfile * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->col[i] != NULL)   {  
        free_DnaProfileCol(obj->col[i]); 
        obj->col[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  DnaProfile_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaProfile_alloc_len(DnaProfileLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfile *]
 *
 */
DnaProfile * DnaProfile_alloc_std(void) 
{
    return DnaProfile_alloc_len(DnaProfileLISTLENGTH);   
}    


/* Function:  DnaProfile_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfile *]
 *
 */
DnaProfile * DnaProfile_alloc_len(int len) 
{
    DnaProfile * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = DnaProfile_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->col = (DnaProfileCol ** ) ckcalloc (len,sizeof(DnaProfileCol *))) == NULL)  {  
      warn("Warning, ckcalloc failed in DnaProfile_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_DnaProfile(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfile *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfile *]
 *
 */
DnaProfile * hard_link_DnaProfile(DnaProfile * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaProfile object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaProfile_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfile *]
 *
 */
DnaProfile * DnaProfile_alloc(void) 
{
    DnaProfile * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProfile *) ckalloc (sizeof(DnaProfile))) == NULL)    {  
      warn("DnaProfile_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->col = NULL; 
    out->len = out->maxlen = 0;  
    out->sa = NULL;  
    out->folded_random = FALSE;  


    return out;  
}    


/* Function:  free_DnaProfile(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfile *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfile *]
 *
 */
DnaProfile * free_DnaProfile(DnaProfile * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProfile obj. Should be trappable");    
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
    if( obj->col != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->col[i] != NULL) 
          free_DnaProfileCol(obj->col[i]);   
        }  
      ckfree(obj->col);  
      }  
    if( obj->sa != NULL) 
      free_SeqAlign(obj->sa);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

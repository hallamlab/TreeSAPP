#ifdef _cplusplus
extern "C" {
#endif
#include "hscore.h"

/* Function:  std_score_Hscore(cut_off,report_stagger)
 *
 * Descrip:    This gives you a standard Hscore
 *             module with a cutoff in score
 *
 *
 * Arg:               cut_off [UNKN ] Undocumented argument [int]
 * Arg:        report_stagger [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Hscore *]
 *
 */
# line 125 "hscore.dy"
Hscore * std_score_Hscore(int cut_off,int report_stagger)
{
  Hscore * out;

  out = Hscore_alloc_std();
  out->his = new_Histogram(-1000,1000,100);
  out->score_level = cut_off;
  out->should_store = raw_should_store_Hscore;
  out->score_to_his = raw_score_to_his;
  out->report_level = report_stagger;

  return out;
}

/* Function:  raw_should_store_Hscore(score,cutoff)
 *
 * Descrip:    This function is for the Hscore std constructor,
 *
 *
 * Arg:         score [UNKN ] Undocumented argument [int]
 * Arg:        cutoff [UNKN ] Undocumented argument [double]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 143 "hscore.dy"
boolean raw_should_store_Hscore(int score,double cutoff)
{
  if( score > cutoff ) {
    return TRUE;
  }
  return FALSE;
}

/* Function:  raw_score_to_his(score)
 *
 * Descrip:    This function is for the Hscore std constructor,
 *
 *
 * Arg:        score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [float]
 *
 */
# line 155 "hscore.dy"
float raw_score_to_his(int score)
{
  return score;
}


/* Function:  should_store_Hscore(hs,score)
 *
 * Descrip:    Tells whether this score should be stored
 *             or not. Also updates Histogram if needed
 *
 *
 * Arg:           hs [UNKN ] Undocumented argument [Hscore *]
 * Arg:        score [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 165 "hscore.dy"
boolean should_store_Hscore(Hscore * hs,int score)
{
  hs->total++;

  if( hs->report_level != -1 && (hs->total % hs->report_level == 0)) {
    if( hs->len > 0) {
      info("Done %d comparisons: last stored comparison was %s to %s",hs->total,hs->ds[hs->len-1]->query->name,hs->ds[hs->len-1]->target->name);
    } else {
      info("Done %d comparisons: No stored comparisons",hs->total);
    }
  }

  if( hs->his != NULL && hs->score_to_his != NULL ) {
    AddToHistogram(hs->his,(*hs->score_to_his)(score));
  }
  if( hs->should_store == NULL ) {
    return TRUE;
  }
  return (*hs->should_store)(score,hs->score_level);
}

 
/* Function:  length_datascore_Hscore(obj)
 *
 * Descrip:    Returns the number of datascores in the hscore
 *             structure
 *
 *
 * Arg:        obj [READ ] Hscore object [Hscore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 194 "hscore.dy"
int length_datascore_Hscore(Hscore * obj)
{
  return obj->len;
}


/* Function:  get_datascore_Hscore(hs,i)
 *
 * Descrip:    Returns the specific datascore held at this
 *             position.
 *
 *             This requires a considerable amount of memory
 *             duplication, so please dont process all your
 *             results by looping through this.
 *
 *
 * Arg:        hs [READ ] Hscore object [Hscore *]
 * Arg:         i [UNKN ] position to be read [int]
 *
 * Return [OWNER]  New datascore object [DataScore *]
 *
 */
# line 213 "hscore.dy"
DataScore * get_datascore_Hscore(Hscore * hs,int i)
{
  DataScore * out;

  out = new_DataScore();
  copy_DataEntry(hs->ds[i]->query,out->query);
  copy_DataEntry(hs->ds[i]->target,out->target);
  out->score = hs->ds[i]->score;
  out->evalue = hs->ds[i]->evalue;
  return out;
}

/* Function:  copy_DataEntry(from,to)
 *
 * Descrip:    Copies the info from one DataEntry to another
 *
 *
 * Arg:        from [UNKN ] Undocumented argument [DataEntry *]
 * Arg:          to [UNKN ] Undocumented argument [DataEntry *]
 *
 */
# line 229 "hscore.dy"
void copy_DataEntry(DataEntry * from,DataEntry * to)
{
  int i;

  to->name = stringalloc(from->name);
  for(i=0;i<DATAENTRYSTDPOINTS;i++) 
    to->data[i] = from->data[i];
  to->is_reversed = from->is_reversed;
  to->byte_position = from->byte_position;
  to->filename = from->filename; /* linked! */
}

/* Function:  get_score_Hscore(hs,i)
 *
 * Descrip: No Description
 *
 * Arg:        hs [READ ] Hscore object [Hscore *]
 * Arg:         i [UNKN ] position to be read [int]
 *
 * Return [UNKN ]  score  [int]
 *
 */
# line 251 "hscore.dy"
int get_score_Hscore(Hscore * hs,int i)
{
  return hs->ds[i]->score;
}


/* Function:  get_evalue_Hscore(hs,i)
 *
 * Descrip:    Returns the evalue of the specific datascore held at this position.
 *
 *
 *
 * Arg:        hs [READ ] Hscore object [Hscore *]
 * Arg:         i [UNKN ] position to be read [int]
 *
 * Return [UNKN ]  evalue  [double]
 *
 */
# line 266 "hscore.dy"
double get_evalue_Hscore(Hscore * hs,int i)
{
  return hs->ds[i]->evalue;
}

/* Function:  fit_Hscore_to_EVD(hs,guess_of_outliers)
 *
 * Descrip:    If a histogram is present, tries to fit the histogram and
 *             then gives evalues to all the scores in the Hscore model
 *
 *
 * Arg:                       hs [UNKN ] Undocumented argument [Hscore *]
 * Arg:        guess_of_outliers [UNKN ] Undocumented argument [float]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 275 "hscore.dy"
boolean fit_Hscore_to_EVD(Hscore * hs,float guess_of_outliers)
{
  int i;


  if( hs->his == NULL ) {
    warn("Your Hscore has no histogram structure, and so no EVD can be fitted");
    return FALSE;
  }
  
  if( ExtremeValueFitHistogram(hs->his,TRUE,guess_of_outliers) == 0 ) {
    warn("Extreme Value distribution is unable to be fitted. Sorry!");
    return FALSE;
  }


  for(i=0;i<hs->len;i++) {
    hs->ds[i]->evalue  = ExtremeValueE((*hs->score_to_his)(hs->ds[i]->score),hs->his->param[EVD_MU],hs->his->param[EVD_LAMBDA],hs->his->total);
  }

  return TRUE;
}

  
/* Function:  minimum_score_Hscore(hs)
 *
 * Descrip:    gets the minimum score from Hscore
 *
 *
 * Arg:        hs [UNKN ] Undocumented argument [Hscore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 302 "hscore.dy"
int minimum_score_Hscore(Hscore * hs)
{
  int i;
  int min;

  if( hs->len == 0) {
    warn("Can't get a minimum score with no entries");
    return 0;
  }

  for(i=1,min=hs->ds[0]->score;i<hs->len;i++) {
    if( min > hs->ds[i]->score ) {
      min = hs->ds[i]->score;
    }
  }

  return min;
}

/* Function:  maximum_score_Hscore(hs)
 *
 * Descrip:    gets the maximum score from Hscore
 *
 *
 * Arg:        hs [UNKN ] Undocumented argument [Hscore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 324 "hscore.dy"
int maximum_score_Hscore(Hscore * hs)
{
  int i;
  int max;

  if( hs->len == 0) {
    warn("Can't get a minimum score with no entries");
    return 0;
  }

  for(i=1,max=hs->ds[0]->score;i<hs->len;i++) {
    if( max < hs->ds[i]->score ) {
      max = hs->ds[i]->score;
    }
  }

  return max;
}
 

/* Function:  basic_show_Hscore(hs,ofp)
 *
 * Descrip:    The most baby-talk showing of Hscore
 *
 *
 * Arg:         hs [UNKN ] Undocumented argument [Hscore *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 349 "hscore.dy"
void basic_show_Hscore(Hscore * hs,FILE * ofp)
{
  int i;

  if( hs == NULL ) {
    warn("parsing in a NULL Hscore object - cannot show!");
    fprintf(ofp,"parsing in a NULL Hscore object - cannot show!");
  }

  for(i=0;i<hs->len;i++) {
    fprintf(ofp,"%3d Query: %12s Target: %12s Score %d\n",i,
	    hs->ds[i]->query->name == NULL ? "NoName" : hs->ds[i]->query->name,
	    hs->ds[i]->target->name == NULL ? "NoName" : hs->ds[i]->target->name,
	    hs->ds[i]->score);
  }
}

/* Function:  sort_Hscore_by_score(hs)
 *
 * Descrip:    As it says, sorts the high score by its score
 *
 *
 * Arg:        hs [UNKN ] Hscore to be sorted [Hscore *]
 *
 */
# line 371 "hscore.dy"
void sort_Hscore_by_score(Hscore * hs)
{
  sort_Hscore(hs,compare_DataScore_by_score);
}


/* Function:  compare_DataScore_by_score(one,two)
 *
 * Descrip:    Used to compare two datascores for
 *             /sort_Hscore_by_score
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [DataScore *]
 * Arg:        two [UNKN ] Undocumented argument [DataScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 383 "hscore.dy"
int compare_DataScore_by_score(DataScore * one,DataScore * two)
{
  if( one->score == two->score ) {
    if( one->query != NULL && one->query->name != NULL && two->query != NULL && two->query->name != NULL )
      return strcmp(one->query->name,two->query->name);
    else return 1;
  }

  return two->score - one->score;
}


/* Function:  new_DataScore(void)
 *
 * Descrip:    The best way to make a new DataScore.
 *             Allocates the query and target DataEntry structures
 *             as well as the DataScore structure.
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DataScore *]
 *
 */
# line 402 "hscore.dy"
DataScore * new_DataScore(void)
{
  DataScore * ds;

  ds = DataScore_alloc();
  ds->query = DataEntry_alloc();
  ds->target = DataEntry_alloc();

  return ds;
}

/* Function:  new_DataScore_from_storage(hs)
 *
 * Descrip:    Gets a new DataScore from Storage
 *
 *
 * Arg:        hs [UNKN ] Undocumented argument [Hscore *]
 *
 * Return [UNKN ]  Undocumented return value [DataScore *]
 *
 */
# line 417 "hscore.dy"
DataScore * new_DataScore_from_storage(Hscore * hs)
{
  DataScoreStorage * new;
  DataScoreStorage * curr;

  if( hs->st_len == 0 ) {
    new = new_DataScoreStorage();
    if( new == NULL ) {
      warn("could not make inital data score storage!");
      return NULL;
    }
    add_st_Hscore(hs,new);
    curr = new;
  } else {
    curr = hs->store[hs->st_len-1];
    if( curr->curr_pos == DATASCORESTORAGE_LENGTH ) {
      new = new_DataScoreStorage();
      if( new == NULL ) {
	warn("could not make data score storage block %d!",hs->st_len-1);
	return NULL;
      }
      add_st_Hscore(hs,new);
      curr = new;
    }
  }

  return &curr->score_array[curr->curr_pos++];
}

      
  
/* Function:  free_DataScore(obj)
 *
 * Descrip:    Correctly handles destruction of a datascore
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [DataScore *]
 *
 * Return [UNKN ]  Undocumented return value [DataScore *]
 *
 */
# line 452 "hscore.dy"
DataScore * free_DataScore(DataScore * obj)
{
  if( obj->is_stored == 1 ) {
    return NULL; /* don't free! */
  }
  if( obj->dynamite_hard_link > 1 ) {
    obj->dynamite_hard_link--;
    return NULL;
  }

  if( obj->query != NULL )
    free_DataEntry(obj->query);

  if( obj->target != NULL )
    free_DataEntry(obj->target);

  return NULL;
}

/* Function:  free_DataScoreStorage(obj)
 *
 * Descrip:    Correctly handles destruction of DataScoreStorage, by
 *             freeing members in data storage
 *
 *
 * Arg:        obj [UNKN ] Undocumented argument [DataScoreStorage *]
 *
 * Return [UNKN ]  Undocumented return value [DataScoreStorage *]
 *
 */
# line 476 "hscore.dy"
DataScoreStorage * free_DataScoreStorage(DataScoreStorage * obj)
{
  int i;

  for(i=0;i<obj->curr_pos;i++) {
    if( obj->query_array[i].name != NULL ) {
      ckfree(obj->query_array[i].name);
    }
    if( obj->target_array[i].name != NULL ) {
      ckfree(obj->target_array[i].name);
    }
  }

  ckfree(obj);

  return NULL;
}
    


/* Function:  new_DataScoreStorage(void)
 *
 * Descrip:    Makes a new DataScoreStorage with all the pointers connected correctly
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DataScoreStorage *]
 *
 */
# line 499 "hscore.dy"
DataScoreStorage * new_DataScoreStorage(void)
{
  DataScoreStorage * out;
  int i;

  out = DataScoreStorage_alloc();
  if( out == NULL ) {
    warn("Unable to make a new DataScoreStorage block with blocksize %d",DATASCORESTORAGE_LENGTH);
    return NULL;
  }

  for(i=0;i<DATASCORESTORAGE_LENGTH;i++) {
    out->score_array[i].query = &out->query_array[i];
    out->score_array[i].target = &out->target_array[i];
    out->score_array[i].is_stored = 1;
  }
    
  return out;
}

# line 490 "hscore.c"
/* Function:  hard_link_DataEntry(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [DataEntry *]
 *
 */
DataEntry * hard_link_DataEntry(DataEntry * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DataEntry object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DataEntry_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DataEntry *]
 *
 */
DataEntry * DataEntry_alloc(void) 
{
    DataEntry * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DataEntry *) ckalloc (sizeof(DataEntry))) == NULL)  {  
      warn("DataEntry_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    /* data[DATAENTRYSTDPOINTS] is an array: no default possible */ 
    out->is_reversed = FALSE;    
    out->byte_position = 0;  


    return out;  
}    


/* Function:  free_DataEntry(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DataEntry *]
 *
 * Return [UNKN ]  Undocumented return value [DataEntry *]
 *
 */
DataEntry * free_DataEntry(DataEntry * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DataEntry obj. Should be trappable"); 
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
    if( obj->name != NULL)   
      ckfree(obj->name);     
    /* obj->filename is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_DataScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DataScore *]
 *
 * Return [UNKN ]  Undocumented return value [DataScore *]
 *
 */
DataScore * hard_link_DataScore(DataScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DataScore object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DataScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DataScore *]
 *
 */
DataScore * DataScore_alloc(void) 
{
    DataScore * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DataScore *) ckalloc (sizeof(DataScore))) == NULL)  {  
      warn("DataScore_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->query = NULL;   
    out->target = NULL;  
    out->score = 0;  
    out->evalue = 0; 
    out->is_stored = 0;  


    return out;  
}    


/* Function:  hard_link_DataScoreStorage(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DataScoreStorage *]
 *
 * Return [UNKN ]  Undocumented return value [DataScoreStorage *]
 *
 */
DataScoreStorage * hard_link_DataScoreStorage(DataScoreStorage * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DataScoreStorage object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DataScoreStorage_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DataScoreStorage *]
 *
 */
DataScoreStorage * DataScoreStorage_alloc(void) 
{
    DataScoreStorage * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DataScoreStorage *) ckalloc (sizeof(DataScoreStorage))) == NULL)    {  
      warn("DataScoreStorage_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* score_array[DATASCORESTORAGE_LENGTH] is an array: no default possible */ 
    /* query_array[DATASCORESTORAGE_LENGTH] is an array: no default possible */ 
    /* target_array[DATASCORESTORAGE_LENGTH] is an array: no default possible */ 
    out->curr_pos = 0;   


    return out;  
}    


/* Function:  swap_Hscore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Hscore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DataScore **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Hscore(DataScore ** list,int i,int j)  
{
    DataScore * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Hscore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Hscore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DataScore **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Hscore(DataScore ** list,int left,int right,int (*comp)(DataScore * ,DataScore * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Hscore(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Hscore (list,++last,i); 
      }  
    swap_Hscore (list,left,last);    
    qsort_Hscore(list,left,last-1,comp); 
    qsort_Hscore(list,last+1,right,comp);    
}    


/* Function:  sort_Hscore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Hscore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Hscore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Hscore(Hscore * obj,int (*comp)(DataScore *, DataScore *)) 
{
    qsort_Hscore(obj->ds,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_Hscore(obj,len)
 *
 * Descrip:    Really an internal function for add_Hscore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Hscore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Hscore(Hscore * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Hscore called with no need"); 
      return TRUE;   
      }  


    if( (obj->ds = (DataScore ** ) ckrealloc (obj->ds,sizeof(DataScore *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Hscore, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Hscore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Hscore *]
 * Arg:        add [OWNER] Object to add to the list [DataScore *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Hscore(Hscore * obj,DataScore * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Hscore(obj,obj->len + HscoreLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->ds[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_Hscore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Hscore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Hscore(Hscore * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->ds[i] != NULL)    {  
        free_DataScore(obj->ds[i]);  
        obj->ds[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_st_Hscore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_st_Hscore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DataScoreStorage **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_st_Hscore(DataScoreStorage ** list,int i,int j)  
{
    DataScoreStorage * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_st_Hscore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_st_Hscore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DataScoreStorage **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_st_Hscore(DataScoreStorage ** list,int left,int right,int (*comp)(DataScoreStorage * ,DataScoreStorage * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_st_Hscore(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_st_Hscore (list,++last,i);  
      }  
    swap_st_Hscore (list,left,last); 
    qsort_st_Hscore(list,left,last-1,comp);  
    qsort_st_Hscore(list,last+1,right,comp); 
}    


/* Function:  sort_st_Hscore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_st_Hscore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Hscore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_st_Hscore(Hscore * obj,int (*comp)(DataScoreStorage *, DataScoreStorage *)) 
{
    qsort_st_Hscore(obj->store,0,obj->st_len-1,comp);    
    return;  
}    


/* Function:  expand_st_Hscore(obj,len)
 *
 * Descrip:    Really an internal function for add_st_Hscore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Hscore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_st_Hscore(Hscore * obj,int len) 
{


    if( obj->st_maxlen > obj->st_len )   {  
      warn("expand_Hscorest_ called with no need");  
      return TRUE;   
      }  


    if( (obj->store = (DataScoreStorage ** ) ckrealloc (obj->store,sizeof(DataScoreStorage *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Hscore, returning FALSE");   
      return FALSE;  
      }  
    obj->st_maxlen = len;    
    return TRUE; 
}    


/* Function:  add_st_Hscore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Hscore *]
 * Arg:        add [OWNER] Object to add to the list [DataScoreStorage *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_st_Hscore(Hscore * obj,DataScoreStorage * add) 
{
    if( obj->st_len >= obj->st_maxlen)   {  
      if( expand_st_Hscore(obj,obj->st_len + HscoreLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->store[obj->st_len++]=add;   
    return TRUE; 
}    


/* Function:  flush_st_Hscore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Hscore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_st_Hscore(Hscore * obj) 
{
    int i;   


    for(i=0;i<obj->st_len;i++)   { /*for i over list length*/ 
      if( obj->store[i] != NULL) {  
        free_DataScoreStorage(obj->store[i]);    
        obj->store[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->st_len = 0; 
    return i;    
}    


/* Function:  Hscore_alloc_std(void)
 *
 * Descrip:    Equivalent to Hscore_alloc_len(HscoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Hscore *]
 *
 */
Hscore * Hscore_alloc_std(void) 
{
    return Hscore_alloc_len(HscoreLISTLENGTH);   
}    


/* Function:  Hscore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Hscore *]
 *
 */
Hscore * Hscore_alloc_len(int len) 
{
    Hscore * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Hscore_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->ds = (DataScore ** ) ckcalloc (len,sizeof(DataScore *))) == NULL)   {  
      warn("Warning, ckcalloc failed in Hscore_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->store = (DataScoreStorage ** ) ckcalloc (len,sizeof(DataScoreStorage *))) == NULL)  {  
      warn("Warning, ckcalloc failed in Hscore_alloc_len");  
      return NULL;   
      }  
    out->st_len = 0; 
    out->st_maxlen = len;    


    return out;  
}    


/* Function:  hard_link_Hscore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Hscore *]
 *
 * Return [UNKN ]  Undocumented return value [Hscore *]
 *
 */
Hscore * hard_link_Hscore(Hscore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Hscore object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Hscore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Hscore *]
 *
 */
Hscore * Hscore_alloc(void) 
{
    Hscore * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Hscore *) ckalloc (sizeof(Hscore))) == NULL)    {  
      warn("Hscore_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->ds = NULL;  
    out->len = out->maxlen = 0;  
    out->store = NULL;   
    out->st_len = out->st_maxlen = 0;    
    out->his = NULL; 
    out->score_level = 0;    
    out->should_store = NULL;    
    out->score_to_his = NULL;    
    out->report_level = 0;   
    out->total = 0;  


    return out;  
}    


/* Function:  free_Hscore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Hscore *]
 *
 * Return [UNKN ]  Undocumented return value [Hscore *]
 *
 */
Hscore * free_Hscore(Hscore * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Hscore obj. Should be trappable");    
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
    if( obj->ds != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->ds[i] != NULL)  
          free_DataScore(obj->ds[i]);    
        }  
      ckfree(obj->ds);   
      }  
    if( obj->store != NULL)  {  
      for(i=0;i<obj->st_len;i++) {  
        if( obj->store[i] != NULL)   
          free_DataScoreStorage(obj->store[i]);  
        }  
      ckfree(obj->store);    
      }  
    if( obj->his != NULL)    
      free_Histogram(obj->his);  
    /* obj->should_store is a function pointer */ 
    /* obj->score_to_his is a function pointer */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  replace_query_DataScore(obj,query)
 *
 * Descrip:    Replace member variable query
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [DataScore *]
 * Arg:        query [OWNER] New value of the variable [DataEntry *]
 *
 * Return [SOFT ]  member variable query [boolean]
 *
 */
boolean replace_query_DataScore(DataScore * obj,DataEntry * query) 
{
    if( obj == NULL)     {  
      warn("In replacement function query for object DataScore, got a NULL object"); 
      return FALSE;  
      }  
    obj->query = query;  
    return TRUE; 
}    


/* Function:  access_query_DataScore(obj)
 *
 * Descrip:    Access member variable query
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DataScore *]
 *
 * Return [SOFT ]  member variable query [DataEntry *]
 *
 */
DataEntry * access_query_DataScore(DataScore * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function query for object DataScore, got a NULL object");    
      return NULL;   
      }  
    return obj->query;   
}    


/* Function:  replace_target_DataScore(obj,target)
 *
 * Descrip:    Replace member variable target
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [DataScore *]
 * Arg:        target [OWNER] New value of the variable [DataEntry *]
 *
 * Return [SOFT ]  member variable target [boolean]
 *
 */
boolean replace_target_DataScore(DataScore * obj,DataEntry * target) 
{
    if( obj == NULL)     {  
      warn("In replacement function target for object DataScore, got a NULL object");    
      return FALSE;  
      }  
    obj->target = target;    
    return TRUE; 
}    


/* Function:  access_target_DataScore(obj)
 *
 * Descrip:    Access member variable target
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DataScore *]
 *
 * Return [SOFT ]  member variable target [DataEntry *]
 *
 */
DataEntry * access_target_DataScore(DataScore * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function target for object DataScore, got a NULL object");   
      return NULL;   
      }  
    return obj->target;  
}    


/* Function:  replace_score_DataScore(obj,score)
 *
 * Descrip:    Replace member variable score
 *             For use principly by API functions
 *
 *
 * Arg:          obj [UNKN ] Object holding the variable [DataScore *]
 * Arg:        score [OWNER] New value of the variable [int]
 *
 * Return [SOFT ]  member variable score [boolean]
 *
 */
boolean replace_score_DataScore(DataScore * obj,int score) 
{
    if( obj == NULL)     {  
      warn("In replacement function score for object DataScore, got a NULL object"); 
      return FALSE;  
      }  
    obj->score = score;  
    return TRUE; 
}    


/* Function:  access_score_DataScore(obj)
 *
 * Descrip:    Access member variable score
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DataScore *]
 *
 * Return [SOFT ]  member variable score [int]
 *
 */
int access_score_DataScore(DataScore * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function score for object DataScore, got a NULL object");    
      return 0;  
      }  
    return obj->score;   
}    


/* Function:  replace_evalue_DataScore(obj,evalue)
 *
 * Descrip:    Replace member variable evalue
 *             For use principly by API functions
 *
 *
 * Arg:           obj [UNKN ] Object holding the variable [DataScore *]
 * Arg:        evalue [OWNER] New value of the variable [double]
 *
 * Return [SOFT ]  member variable evalue [boolean]
 *
 */
boolean replace_evalue_DataScore(DataScore * obj,double evalue) 
{
    if( obj == NULL)     {  
      warn("In replacement function evalue for object DataScore, got a NULL object");    
      return FALSE;  
      }  
    obj->evalue = evalue;    
    return TRUE; 
}    


/* Function:  access_evalue_DataScore(obj)
 *
 * Descrip:    Access member variable evalue
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DataScore *]
 *
 * Return [SOFT ]  member variable evalue [double]
 *
 */
double access_evalue_DataScore(DataScore * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function evalue for object DataScore, got a NULL object");   
      return 0;  
      }  
    return obj->evalue;  
}    


/* Function:  replace_name_DataEntry(obj,name)
 *
 * Descrip:    Replace member variable name
 *             For use principly by API functions
 *
 *
 * Arg:         obj [UNKN ] Object holding the variable [DataEntry *]
 * Arg:        name [OWNER] New value of the variable [char *]
 *
 * Return [SOFT ]  member variable name [boolean]
 *
 */
boolean replace_name_DataEntry(DataEntry * obj,char * name) 
{
    if( obj == NULL)     {  
      warn("In replacement function name for object DataEntry, got a NULL object");  
      return FALSE;  
      }  
    obj->name = name;    
    return TRUE; 
}    


/* Function:  access_name_DataEntry(obj)
 *
 * Descrip:    Access member variable name
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DataEntry *]
 *
 * Return [SOFT ]  member variable name [char *]
 *
 */
char * access_name_DataEntry(DataEntry * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function name for object DataEntry, got a NULL object"); 
      return NULL;   
      }  
    return obj->name;    
}    


/* Function:  replace_is_reversed_DataEntry(obj,is_reversed)
 *
 * Descrip:    Replace member variable is_reversed
 *             For use principly by API functions
 *
 *
 * Arg:                obj [UNKN ] Object holding the variable [DataEntry *]
 * Arg:        is_reversed [OWNER] New value of the variable [boolean]
 *
 * Return [SOFT ]  member variable is_reversed [boolean]
 *
 */
boolean replace_is_reversed_DataEntry(DataEntry * obj,boolean is_reversed) 
{
    if( obj == NULL)     {  
      warn("In replacement function is_reversed for object DataEntry, got a NULL object");   
      return FALSE;  
      }  
    obj->is_reversed = is_reversed;  
    return TRUE; 
}    


/* Function:  access_is_reversed_DataEntry(obj)
 *
 * Descrip:    Access member variable is_reversed
 *             For use principly by API functions
 *
 *
 * Arg:        obj [UNKN ] Object holding the variable [DataEntry *]
 *
 * Return [SOFT ]  member variable is_reversed [boolean]
 *
 */
boolean access_is_reversed_DataEntry(DataEntry * obj) 
{
    if( obj == NULL)     {  
      warn("In accessor function is_reversed for object DataEntry, got a NULL object");  
      return FALSE;  
      }  
    return obj->is_reversed;     
}    



#ifdef _cplusplus
}
#endif

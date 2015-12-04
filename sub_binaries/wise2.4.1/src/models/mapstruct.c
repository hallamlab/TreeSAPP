#ifdef _cplusplus
extern "C" {
#endif
#include "mapstruct.h"


# line 33 "mapstruct.dy"
MappedCloneMatch * new_MappedCloneMatch(MappedCloneSet * iset,MappedCloneSet * jset,int match,int mismatch)
{
  int i;
  int j;
  int ii;
  int jj;
  int k;
  MappedCloneMatch * out;
  int * buffer;

  assert(match >= 0 );
  assert(mismatch <= 0);
  assert(iset);
  assert(jset);


  out = MappedCloneMatch_alloc_matrix(iset->length,jset->length);

  
  buffer = calloc(jset->length,sizeof(int));


  for(i=0;i<iset->len;i++) {

    /* no point even noticing not seen clones */
    if( iset->clone[i]->seen == 0 ) {
      continue;
    }

    for(k=0;k<jset->length;k++) 
      buffer[k] = 0;

    for(j=0;j<jset->len;j++) {
      /* positive scores are easy, as we only have to loop over names */      
      if( strcmp(iset->clone[i]->clone_name,jset->clone[j]->clone_name) == 0 ) {
	for(ii=iset->clone[i]->start;ii<iset->clone[i]->end;ii++) {
	  for(jj=jset->clone[j]->start;jj<jset->clone[j]->end;jj++) {
	    out->matrix[ii][jj] += match;
	    buffer[jj] = 1;
	  }
	}
      }
    }

    /* now handle negative scores. in the j dimension, where buffer == 0, these
       regions do not have this particular i clone. As long as this i clone is actually
       seen (true due to first continue if) then we can substract the mismatch across the i
    */

    for(k=0;k<jset->length;k++) {
      if( buffer[k] == 0 ) {
	for(ii=iset->clone[i]->start;ii<iset->clone[i]->end;ii++) {
	  out->matrix[ii][k] += mismatch;
	}
      }
    }

  }
  

  out->skip_iset = (int*) calloc(iset->length,sizeof(int));
  out->skip_jset = (int*) calloc(jset->length,sizeof(int));

  for(i=0;i<iset->len;i++) {
    for(ii=iset->clone[i]->start;ii<iset->clone[i]->end;ii++) {
      if( iset->clone[i]->seen == 1 ) {
	out->skip_iset += mismatch;
      }
    }
  }

  for(j=0;j<jset->len;j++) {
    for(jj=jset->clone[j]->start;jj<jset->clone[j]->end;jj++) {
      if( jset->clone[j]->seen == 1 ) {
	out->skip_jset += mismatch;
      }
    }
  }


  return out;
}


# line 117 "mapstruct.dy"
int MappedCloneSet_skip(MappedCloneSet * s,int pos,int skip_cost)
{
  int i;
  int score = 0;
  MappedCloneSet * mcs;

  mcs = subsection_MappedCloneSet(s,pos,pos,1);

  free_MappedCloneSet(mcs);

  return skip_cost*mcs->len;
}

  
# line 131 "mapstruct.dy"
int MappedCloneSet_match(MappedCloneSet * weak_query,MappedCloneSet * trusted_target,int qpos,int tpos,int spread,int match,int mismatch)
{
  MappedCloneSet * weak_slice;
  MappedCloneSet * trusted_slice;
  int i;
  int j;
  int score =0;

  weak_slice    = subsection_MappedCloneSet(weak_query,qpos-spread,qpos+spread,1);
  trusted_slice = subsection_MappedCloneSet(trusted_target,tpos-spread,tpos+spread,1);

  if( weak_slice->len == 0 || trusted_slice->len == 0 ) {
    score = mismatch;
  } else {
    for(i=0;i<weak_slice->len;i++) {
      for(j=0;j<trusted_slice->len;j++) {
	if( strcmp(weak_slice->clone[i]->clone_name,trusted_slice->clone[j]->clone_name) == 0 ) {
	  score += match;
	}
      }
    }
  }

  free_MappedCloneSet(weak_slice);
  free_MappedCloneSet(trusted_slice);

  return score;
}



# line 162 "mapstruct.dy"
int old_MappedCloneSet_match(MappedCloneSet * one,MappedCloneSet * two,int qpos,int tpos,int spread,int match,int mismatch)
{
  int i;
  int startj;
  int j;
  int score = 0;
  int has_matched;

  /* sorted by start. If positions are before start - return 0 */

  if( one->clone[0]->start-spread > qpos ) {
    return mismatch;
  }

  if( two->clone[0]->start-spread > tpos ) {
    return mismatch;
  }


  for(i=0;i<one->len;i++) {
    if( one->clone[i]->end+spread >= qpos ) {
      break;
    }
  }

  for(startj=0;startj<two->len;startj++) {
    if( two->clone[startj]->end+spread >= tpos ) {
      break;
    }
  }

  if( i >= one->len ) {
    return mismatch;
  }

  if( startj >= two->len ) {
    return mismatch;
  }


  for(;i<one->len && one->clone[i]->start-spread <= qpos;i++) {
    if( one->clone[i]->seen == 0 ) {
      continue;
    }
    has_matched = 0;
    for(j=startj;j<two->len && j < two->clone[j]->start-spread < tpos;j++) {
      if( two->clone[j]->seen == 0 ) {
	continue;
      }
      if( strcmp(two->clone[j]->clone_name,one->clone[i]->clone_name) == 0 ) {
	has_matched = 1;
	break;
      }
    }
    if( has_matched == 1 ) {
      score += match;    
    } else {
      score -= mismatch;
    }
  }

  return score;
}


/* Function:  synchronise_MappedCloneSets(one,two)
 *
 * Descrip:    updates the internal seen flags for the clone sets in
 *             preparation for the dp
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:        two [UNKN ] Undocumented argument [MappedCloneSet *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 231 "mapstruct.dy"
boolean synchronise_MappedCloneSets(MappedCloneSet * one,MappedCloneSet * two)
{
  int i;
  MappedClone * mc;


  assert(one);
  assert(two);

  for(i=0;i<one->len;i++) {
    mc = find_named_MappedClone(two,one->clone[i]->clone_name);
    if( mc != NULL ) {
      mc->seen = 1;
      one->clone[i]->seen = 1;
    }
  }

  return TRUE;
}
      

/* Function:  subsection_MappedCloneSet(mcs,coord_start,coord_end,only_seen)
 *
 * Descrip:    Returns a sub-section of the MappedClone 
 *
 *
 * Arg:                mcs [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:        coord_start [UNKN ] Undocumented argument [int]
 * Arg:          coord_end [UNKN ] Undocumented argument [int]
 * Arg:          only_seen [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
# line 255 "mapstruct.dy"
MappedCloneSet * subsection_MappedCloneSet(MappedCloneSet * mcs,int coord_start,int coord_end,int only_seen)
{
  MappedCloneSet * out;
  int i;


  out = MappedCloneSet_alloc_std();

  for(i=0;i<mcs->len;i++) {
    if( mcs->clone[i]->end >= coord_start ) {
      break;
    }
  }

  if( i >= mcs->len ) {
    return out;
  }

  for(i=0;i<mcs->len;i++) {
    if( only_seen == 1 && mcs->clone[i]->seen == 0 ) {
      continue;
    }

    if( !(mcs->clone[i]->end < coord_start || mcs->clone[i]->start > coord_end) ) {
      add_MappedCloneSet(out,hard_link_MappedClone(mcs->clone[i]));
    }
    if( mcs->clone[i]->start > coord_end ) {
      break;
    }
  }

  return out;
}

  

/* Function:  find_named_MappedClone(mcs,clone_name)
 *
 * Descrip:    Finds a mapped clone set with this name
 *
 *
 * Arg:               mcs [UNKN ] Undocumented argument [MappedCloneSet *]
 * Arg:        clone_name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [MappedClone *]
 *
 */
# line 294 "mapstruct.dy"
MappedClone * find_named_MappedClone(MappedCloneSet * mcs,char * clone_name)
{
  int i;

  /*we should have hashed. Oooops */

  for(i=0;i<mcs->len;i++) {
    if( strcmp(mcs->clone[i]->clone_name,clone_name) == 0 ) {
      return mcs->clone[i];
    }
  }

  return NULL;
}

/* Function:  start_comp_MappedClone(a,b)
 *
 * Descrip:    sorting for MappedClones
 *
 *
 * Arg:        a [UNKN ] Undocumented argument [MappedClone *]
 * Arg:        b [UNKN ] Undocumented argument [MappedClone *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 312 "mapstruct.dy"
int start_comp_MappedClone(MappedClone * a,MappedClone * b)
{
  if( a->start >= b->start ) {
    return 1;
  } else {
    return -1;
  }
}
    

/* Function:  read_MappedCloneSet(ifp)
 *
 * Descrip:    Reads in a MappedCloneSet file
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
# line 325 "mapstruct.dy"
MappedCloneSet * read_MappedCloneSet(FILE * ifp) 
{
  MappedCloneSet * out;
  MappedClone * temp;
  char buffer[512];

  out = MappedCloneSet_alloc_std();

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strstartcmp(buffer,"#") == 0 ) {
      continue;
    }
    
    temp = read_MappedClone_line(buffer);
    if( temp == NULL ) { 
      continue;
    }

    add_MappedCloneSet(out,temp);
  }

  sort_MappedCloneSet(out,start_comp_MappedClone);

  out->length = out->clone[out->len-1]->end;

  return out;
}


/* Function:  read_MappedClone_line(line)
 *
 * Descrip:    Provides a mapped clone from a name\tstart\tend format
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [MappedClone *]
 *
 */
# line 357 "mapstruct.dy"
MappedClone * read_MappedClone_line(char * line)
{
  MappedClone * out;
  char * a;
  char * b;
  char * c;
  char * d;
  char * e;

  char * copy;

  copy = stringalloc(line);
  out = MappedClone_alloc();

  a = strtok(line,spacestr);
  b = strtok(NULL,spacestr);
  c = strtok(NULL,spacestr);
  d = strtok(NULL,spacestr);
  e = strtok(NULL,spacestr);


  if( a == NULL || b == NULL || c == NULL || d == NULL || e == NULL ) {
    warn("Bad clone line %s",copy);
    ckfree(copy);
    return NULL;
  }



  out->start = atol(a);
  out->end   = atol(b);
  out->clone_name  = stringalloc(c);
  out->accession   = stringalloc(d);
  out->contig      = stringalloc(e);

  ckfree(copy);
  return out;
}


  
# line 417 "mapstruct.c"
/* Function:  hard_link_MappedClone(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MappedClone *]
 *
 * Return [UNKN ]  Undocumented return value [MappedClone *]
 *
 */
MappedClone * hard_link_MappedClone(MappedClone * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MappedClone object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MappedClone_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MappedClone *]
 *
 */
MappedClone * MappedClone_alloc(void) 
{
    MappedClone * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MappedClone *) ckalloc (sizeof(MappedClone))) == NULL)  {  
      warn("MappedClone_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->clone_name = NULL;  
    out->accession = NULL;   
    out->contig = NULL;  
    out->start = 0;  
    out->end = 0;    
    out->seen = 0;   


    return out;  
}    


/* Function:  free_MappedClone(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MappedClone *]
 *
 * Return [UNKN ]  Undocumented return value [MappedClone *]
 *
 */
MappedClone * free_MappedClone(MappedClone * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MappedClone obj. Should be trappable");   
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
    if( obj->clone_name != NULL) 
      ckfree(obj->clone_name);   
    if( obj->accession != NULL)  
      ckfree(obj->accession);    
    if( obj->contig != NULL) 
      ckfree(obj->contig);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_MappedCloneSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_MappedCloneSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [MappedClone **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_MappedCloneSet(MappedClone ** list,int i,int j)  
{
    MappedClone * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_MappedCloneSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_MappedCloneSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [MappedClone **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_MappedCloneSet(MappedClone ** list,int left,int right,int (*comp)(MappedClone * ,MappedClone * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_MappedCloneSet(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_MappedCloneSet (list,++last,i); 
      }  
    swap_MappedCloneSet (list,left,last);    
    qsort_MappedCloneSet(list,left,last-1,comp); 
    qsort_MappedCloneSet(list,last+1,right,comp);    
}    


/* Function:  sort_MappedCloneSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_MappedCloneSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [MappedCloneSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_MappedCloneSet(MappedCloneSet * obj,int (*comp)(MappedClone *, MappedClone *)) 
{
    qsort_MappedCloneSet(obj->clone,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_MappedCloneSet(obj,len)
 *
 * Descrip:    Really an internal function for add_MappedCloneSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MappedCloneSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_MappedCloneSet(MappedCloneSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_MappedCloneSet called with no need"); 
      return TRUE;   
      }  


    if( (obj->clone = (MappedClone ** ) ckrealloc (obj->clone,sizeof(MappedClone *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_MappedCloneSet, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_MappedCloneSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MappedCloneSet *]
 * Arg:        add [OWNER] Object to add to the list [MappedClone *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_MappedCloneSet(MappedCloneSet * obj,MappedClone * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_MappedCloneSet(obj,obj->len + MappedCloneSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->clone[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_MappedCloneSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MappedCloneSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_MappedCloneSet(MappedCloneSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->clone[i] != NULL) {  
        free_MappedClone(obj->clone[i]); 
        obj->clone[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  MappedCloneSet_alloc_std(void)
 *
 * Descrip:    Equivalent to MappedCloneSet_alloc_len(MappedCloneSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * MappedCloneSet_alloc_std(void) 
{
    return MappedCloneSet_alloc_len(MappedCloneSetLISTLENGTH);   
}    


/* Function:  MappedCloneSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * MappedCloneSet_alloc_len(int len) 
{
    MappedCloneSet * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = MappedCloneSet_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->clone = (MappedClone ** ) ckcalloc (len,sizeof(MappedClone *))) == NULL)    {  
      warn("Warning, ckcalloc failed in MappedCloneSet_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_MappedCloneSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MappedCloneSet *]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * hard_link_MappedCloneSet(MappedCloneSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MappedCloneSet object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MappedCloneSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * MappedCloneSet_alloc(void) 
{
    MappedCloneSet * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MappedCloneSet *) ckalloc (sizeof(MappedCloneSet))) == NULL)    {  
      warn("MappedCloneSet_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->clone = NULL;   
    out->len = out->maxlen = 0;  
    out->length = 0; 
    out->cursor = 0; 


    return out;  
}    


/* Function:  free_MappedCloneSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MappedCloneSet *]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneSet *]
 *
 */
MappedCloneSet * free_MappedCloneSet(MappedCloneSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MappedCloneSet obj. Should be trappable");    
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
    if( obj->clone != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->clone[i] != NULL)   
          free_MappedClone(obj->clone[i]);   
        }  
      ckfree(obj->clone);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  MappedCloneMatch_alloc_matrix(leni,lenj)
 *
 * Descrip:    Allocates structure and matrix
 *
 *
 * Arg:        leni [UNKN ] Length of first dimension of matrix [int]
 * Arg:        lenj [UNKN ] Length of second dimension of matrix [int]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneMatch *]
 *
 */
MappedCloneMatch * MappedCloneMatch_alloc_matrix(int leni,int lenj) 
{
    MappedCloneMatch * out; /* out is exported          */ 
    register int i; /* for stepping down matrix */ 
    register int j; /* for stepping across matrix */ 


    /* Call alloc function, return NULL if NULL */ 
    if((out = MappedCloneMatch_alloc()) == NULL) 
      return NULL;   


    /* Allocate memory for matrix  */ 
    if((out->matrix = (int **) ckcalloc (leni,sizeof(int *))) == NULL)   {  
      warn("Memory allocation problem in matrix for MappedCloneMatch matrix, first pointer set");    
      ckfree(out);   
      return NULL;   
      }  


    /* Add NULL to all matrix pointers so free can be called */ 
    for(i=0;i<leni;i++)  
      out->matrix[i] = NULL;     


    /* Allocate each matrix row */ 
    for(i=0;i<leni;i++)  {  
      out->matrix[i] = (int *) ckcalloc (lenj,sizeof(int ));     
      if( out->matrix[i] == NULL)    {  
        warn("Failed alloc on %d, calling free and returning NULL",i);   
        free_MappedCloneMatch(out);  
        return NULL;     
        }  
      }  


    for(i=0;i<leni;i++)  {  
      for(j=0;j<lenj;j++)    
        out->matrix[i][j] = 0;   
      }  


    out->leni=out->maxleni=leni;     
    out->lenj=out->maxlenj=lenj;     


    return out;  
}    


/* Function:  expand_MappedCloneMatch(obj,leni,lenj)
 *
 * Descrip:    Expands matrix. Rarely used
 *
 *
 * Arg:         obj [UNKN ] Undocumented argument [MappedCloneMatch *]
 * Arg:        leni [UNKN ] Undocumented argument [int]
 * Arg:        lenj [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_MappedCloneMatch(MappedCloneMatch * obj,int leni,int lenj) 
{
    int i;   
    int actualj;     


    if( obj == NULL)     {  
      warn("Trying to expand a MappedCloneMatch but is NULL!");  
      return FALSE;  
      }  


    if( leni <= obj->maxleni && lenj <= obj->maxlenj)    
      return TRUE;   


    if( obj->maxleni < leni )    {  
      if( (obj->matrix=(int **) ckrealloc (obj->matrix,sizeof(int *)*leni)) == NULL) 
        return FALSE;    
      obj->maxleni=obj->leni=leni;   
      }  
    if( lenj > obj->maxlenj )    
      actualj = lenj;    
    else actualj = obj->maxlenj; 
    for(i=0;i<obj->leni;i++) {  
      if((obj->matrix[i] = (int *) realloc (obj->matrix[i],sizeof(int ) * actualj)) == NULL) 
        return FALSE;    
      }  
    return TRUE;     
}    


/* Function:  hard_link_MappedCloneMatch(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MappedCloneMatch *]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneMatch *]
 *
 */
MappedCloneMatch * hard_link_MappedCloneMatch(MappedCloneMatch * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MappedCloneMatch object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MappedCloneMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneMatch *]
 *
 */
MappedCloneMatch * MappedCloneMatch_alloc(void) 
{
    MappedCloneMatch * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MappedCloneMatch *) ckalloc (sizeof(MappedCloneMatch))) == NULL)    {  
      warn("MappedCloneMatch_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->matrix = NULL;  
    out->leni=out->maxleni=0;    
    out->lenj=out->maxlenj=0;    
    out->skip_iset = NULL;   
    out->skip_jset = NULL;   


    return out;  
}    


/* Function:  free_MappedCloneMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MappedCloneMatch *]
 *
 * Return [UNKN ]  Undocumented return value [MappedCloneMatch *]
 *
 */
MappedCloneMatch * free_MappedCloneMatch(MappedCloneMatch * obj) 
{
    int return_early = 0;    
    int i;   
    int j;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MappedCloneMatch obj. Should be trappable");  
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
    if( obj->matrix != NULL) {  
      for(i=0;i<obj->leni;i++)   {  
        if( obj->matrix[i] != NULL)  
          ckfree(obj->matrix[i]);    
        }  
      ckfree(obj->matrix);   
      }  
    if( obj->skip_iset != NULL)  
      ckfree(obj->skip_iset);    
    if( obj->skip_jset != NULL)  
      ckfree(obj->skip_jset);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

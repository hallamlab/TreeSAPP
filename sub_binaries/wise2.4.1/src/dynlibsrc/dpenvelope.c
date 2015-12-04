#ifdef _cplusplus
extern "C" {
#endif
#include "dpenvelope.h"

/* Function:  overlap_DPUnit(a,b)
 *
 * Descrip:    Helper function that checks whether things overlap or not
 *
 *
 * Arg:        a [UNKN ] Undocumented argument [DPUnit *]
 * Arg:        b [UNKN ] Undocumented argument [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 37 "dpenvelope.dy"
boolean overlap_DPUnit(DPUnit * a,DPUnit * b)
{
  int diag_a;
  int diag_b;

  if( a->type == DPENV_DIAG && b->type == DPENV_DIAG ) {
    diag_a = a->starti - a->startj;
    diag_b = b->starti - b->startj;

    if( diag_a - a->height > diag_b + b->height || diag_a + a->height < diag_b + b->height ) {
      return FALSE;
    }

    diag_a = a->starti + a->startj;
    diag_b = b->starti + b->startj;

    if( diag_a + a->length < diag_b || diag_a > diag_b + b->length ) {
      return FALSE;
    }

    return TRUE;

  }

  if( a->type == DPENV_RECT && b->type == DPENV_RECT ) {
    fatal("Not implemented rectangle overlap!");
  }


  fatal("Not implemented rectangle - diag overlap!");

}

/* Function:  read_DPEnvelope_file(filename)
 *
 * Descrip:    Helper function that also opens the filename
 *
 *
 * Arg:        filename [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
# line 73 "dpenvelope.dy"
DPEnvelope * read_DPEnvelope_file(char * filename)
{
  FILE * ifp;
  DPEnvelope * dpenv;

  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    error("Could not open file with %s",filename);
    return NULL;
  }

  dpenv = read_DPEnvelope(ifp);

  fclose(ifp);

  return dpenv;
}


/* Function:  read_DPEnvelope(ifp)
 *
 * Descrip:    Reads a DPEnvelope from a file
 *
 *
 * Arg:        ifp [UNKN ] Undocumented argument [FILE *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
# line 95 "dpenvelope.dy"
DPEnvelope * read_DPEnvelope(FILE * ifp)
{
  char buffer[MAXLINE];
  char ** base;
  char ** str;
  DPEnvelope * out;
  DPUnit * unit;

  out = DPEnvelope_alloc_std();



  while( fgets(buffer,MAXLINE,ifp) != NULL ) {

    base = str = breakstring(buffer,spacestr);    
    unit = DPUnit_alloc();
    add_DPEnvelope(out,unit);

    if( strcmp(*str,"rect") == 0 ) {
      unit->type = DPENV_RECT;
    } else if ( strcmp(*str,"diag") == 0 ) {
      unit->type = DPENV_DIAG;
    } else {
      error("Cannot parse DPEnv file");
      continue;
    }

    str++;
    if( *str == NULL ) {
      error("Cannot parse DPEnv file");
      continue;
    }
    unit->starti = atoi(*str);
    str++;
    if( *str == NULL ) {
      error("Cannot parse DPEnv file");
      continue;
    }
    unit->startj = atoi(*str);
    
    str++;
    if( *str == NULL ) {
      error("Cannot parse DPEnv file");
      continue;
    }
    unit->height = atoi(*str);
    
    str++;
    if( *str == NULL ) {
      error("Cannot parse DPEnv file");
      continue;
    }
    unit->length = atoi(*str);
  }

  return out;
}


/* Function:  show_DPEnvelope(dpe,ofp)
 *
 * Descrip:    shows structure. useful for debugging
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 157 "dpenvelope.dy"
void show_DPEnvelope(DPEnvelope * dpe,FILE * ofp)
{
  int i;

  for(i=0;i<dpe->len;i++) {
    fprintf(ofp,"Unit %d [%s] Start %d-%d Height: %d Length: %d\n",i,dpe->dpu[i]->type == DPENV_RECT ? "rect" : "diag",dpe->dpu[i]->starti,dpe->dpu[i]->startj,dpe->dpu[i]->height,dpe->dpu[i]->length);
  }

}

/* Function:  is_in_DPEnvelope(dpe,i,j)
 *
 * Descrip:    Tests whether this i,j position is allowed in the
 *             DPEnvelope
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          i [UNKN ] Undocumented argument [int]
 * Arg:          j [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 171 "dpenvelope.dy"
boolean is_in_DPEnvelope(DPEnvelope * dpe,int i,int j)
{
  int k;

  for(k=0;k<dpe->len;k++) {
    auto DPUnit * u;
    u = dpe->dpu[k];

    switch (u->type) {

    case DPENV_RECT :
      if( i >= u->starti && j >= u->startj && i <= (u->starti+u->height) && j <= (u->startj+u->length) ) 
	return TRUE;
      else 
	break;
    case DPENV_DIAG :
      if( abs( (i-j) - (u->starti-u->startj)) <= u->height && 
	i+j >= u->starti+u->startj && i+j+u->length >= u->starti+u->startj) 
	return TRUE;
      break;

    default :
      warn("Bad DPUnit type put in. Yuk. Bad error... %d",u->type);
      return FALSE;
    }
  }


  return FALSE;
}

/* Function:  prepare_DPEnvelope(dpe)
 *
 * Descrip:    Should run this before using the DPEnvelope
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 205 "dpenvelope.dy"
boolean prepare_DPEnvelope(DPEnvelope * dpe)
{
  int i;
  
  dpe->starti = 1000000;
  dpe->startj = 1000000;
  dpe->endi   = 1;
  dpe->endj   = 1;

  sort_DPEnvelope_by_startj(dpe);

  for(i=0;i<dpe->len;i++) {
    if( dpe->dpu[i]->type == DPENV_RECT ) {
      if( dpe->starti > dpe->dpu[i]->starti-1 ) {
	dpe->starti = dpe->dpu[i]->starti-1;
      }
      if( dpe->startj > dpe->dpu[i]->startj-1 ) {
	dpe->startj = dpe->dpu[i]->startj-1;
      }
      if( dpe->endi < dpe->dpu[i]->starti+ dpe->dpu[i]->height ) {
	dpe->endi = dpe->dpu[i]->starti+ dpe->dpu[i]->height;
      }
      if( dpe->endj < dpe->dpu[i]->startj + dpe->dpu[i]->length ) {
	dpe->endj = dpe->dpu[i]->startj + dpe->dpu[i]->length;
      }
    } else { /* DIAG */
      if( dpe->starti > dpe->dpu[i]->starti-dpe->dpu[i]->height ) {
	dpe->starti = dpe->dpu[i]->starti-dpe->dpu[i]->height;
      }
      if( dpe->startj > dpe->dpu[i]->startj-dpe->dpu[i]->height ) {
	dpe->startj = dpe->dpu[i]->startj-dpe->dpu[i]->height;
      }
      if( dpe->endi < dpe->dpu[i]->starti+ dpe->dpu[i]->length+dpe->dpu[i]->height ) {
	dpe->endi = dpe->dpu[i]->starti+ dpe->dpu[i]->length + dpe->dpu[i]->height;
      }
      if( dpe->endj < dpe->dpu[i]->startj + dpe->dpu[i]->length + dpe->dpu[i]->height ) {
	dpe->endj = dpe->dpu[i]->startj + dpe->dpu[i]->length + dpe->dpu[i]->height;
      }
    }
  }

  return TRUE;
}

/* Function:  sort_DPEnvelope_by_startj(dpe)
 *
 * Descrip:    Sorts by startj
 *
 *
 * Arg:        dpe [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
# line 252 "dpenvelope.dy"
void sort_DPEnvelope_by_startj(DPEnvelope * dpe)
{
  sort_DPEnvelope(dpe,compare_DPUnit_startj);
}

/* Function:  compare_DPUnit_startj(one,two)
 *
 * Descrip:    internal for sort by startj
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [DPUnit *]
 * Arg:        two [UNKN ] Undocumented argument [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 261 "dpenvelope.dy"
int compare_DPUnit_startj(DPUnit * one,DPUnit * two)
{
  return one->startj  - two->startj;
}







# line 291 "dpenvelope.c"
/* Function:  hard_link_DPUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DPUnit *]
 *
 */
DPUnit * hard_link_DPUnit(DPUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DPUnit object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DPUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPUnit *]
 *
 */
DPUnit * DPUnit_alloc(void) 
{
    DPUnit * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DPUnit *) ckalloc (sizeof(DPUnit))) == NULL)    {  
      warn("DPUnit_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = DPENV_RECT;  
    out->starti = 0; 
    out->startj = 0; 
    out->height = 0; 
    out->length = 0; 


    return out;  
}    


/* Function:  free_DPUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [DPUnit *]
 *
 */
DPUnit * free_DPUnit(DPUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DPUnit obj. Should be trappable");    
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


/* Function:  swap_DPEnvelope(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_DPEnvelope
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DPUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_DPEnvelope(DPUnit ** list,int i,int j)  
{
    DPUnit * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_DPEnvelope(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_DPEnvelope which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DPUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_DPEnvelope(DPUnit ** list,int left,int right,int (*comp)(DPUnit * ,DPUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_DPEnvelope(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_DPEnvelope (list,++last,i); 
      }  
    swap_DPEnvelope (list,left,last);    
    qsort_DPEnvelope(list,left,last-1,comp); 
    qsort_DPEnvelope(list,last+1,right,comp);    
}    


/* Function:  sort_DPEnvelope(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_DPEnvelope
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DPEnvelope *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_DPEnvelope(DPEnvelope * obj,int (*comp)(DPUnit *, DPUnit *)) 
{
    qsort_DPEnvelope(obj->dpu,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_DPEnvelope(obj,len)
 *
 * Descrip:    Really an internal function for add_DPEnvelope
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DPEnvelope *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DPEnvelope(DPEnvelope * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_DPEnvelope called with no need"); 
      return TRUE;   
      }  


    if( (obj->dpu = (DPUnit ** ) ckrealloc (obj->dpu,sizeof(DPUnit *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_DPEnvelope, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_DPEnvelope(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DPEnvelope *]
 * Arg:        add [OWNER] Object to add to the list [DPUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_DPEnvelope(DPEnvelope * obj,DPUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_DPEnvelope(obj,obj->len + DPEnvelopeLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->dpu[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_DPEnvelope(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DPEnvelope(DPEnvelope * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->dpu[i] != NULL)   {  
        free_DPUnit(obj->dpu[i]);    
        obj->dpu[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  DPEnvelope_alloc_std(void)
 *
 * Descrip:    Equivalent to DPEnvelope_alloc_len(DPEnvelopeLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * DPEnvelope_alloc_std(void) 
{
    return DPEnvelope_alloc_len(DPEnvelopeLISTLENGTH);   
}    


/* Function:  DPEnvelope_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * DPEnvelope_alloc_len(int len) 
{
    DPEnvelope * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = DPEnvelope_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->dpu = (DPUnit ** ) ckcalloc (len,sizeof(DPUnit *))) == NULL)    {  
      warn("Warning, ckcalloc failed in DPEnvelope_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_DPEnvelope(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * hard_link_DPEnvelope(DPEnvelope * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DPEnvelope object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DPEnvelope_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * DPEnvelope_alloc(void) 
{
    DPEnvelope * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DPEnvelope *) ckalloc (sizeof(DPEnvelope))) == NULL)    {  
      warn("DPEnvelope_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->dpu = NULL; 
    out->len = out->maxlen = 0;  
    out->bbox = NULL;    
    out->starti = 0; 
    out->startj = 0; 
    out->endi = 0;   
    out->endj = 0;   


    return out;  
}    


/* Function:  free_DPEnvelope(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [DPEnvelope *]
 *
 */
DPEnvelope * free_DPEnvelope(DPEnvelope * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DPEnvelope obj. Should be trappable");    
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
    if( obj->dpu != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->dpu[i] != NULL) 
          free_DPUnit(obj->dpu[i]);  
        }  
      ckfree(obj->dpu);  
      }  
    if( obj->bbox != NULL)   
      free_DPUnit(obj->bbox);    


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

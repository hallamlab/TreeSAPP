#ifdef _cplusplus
extern "C" {
#endif
#include "labelmaster.h"

# line 43 "labelmaster.dy"
char * quoted_string_from_Label(Label ** list,int len)
{
  char buf[512];
  char buffer[2048];
  int i;

  sprintf(buffer,"\"%s\"",list[0]->name);

  for(i=1;i<len;i++) {
    sprintf(buf,",\"%s\"",list[i]->name);
    strcat(buffer,buf);

  }

  return stringalloc(buffer);
}

# line 60 "labelmaster.dy"
char * query_quoted_string_from_LabelMaster(LabelMaster * lm)
{
  return quoted_string_from_Label(lm->query,lm->query_len);
}

# line 65 "labelmaster.dy"
char * target_quoted_string_from_LabelMaster(LabelMaster * lm)
{
  return quoted_string_from_Label(lm->target,lm->target_len);
}


# line 71 "labelmaster.dy"
int index_for_label(char * name,Label ** list,int len)
{
  int i;

  for(i=0;i<len;i++) {
    if( strcmp(list[i]->name,name) == 0)
      return i;
  }

  return -1;
}


# line 84 "labelmaster.dy"
int index_for_query_label(char * name,LabelMaster * lm)
{
  return index_for_label(name,lm->query,lm->query_len);
}

# line 89 "labelmaster.dy"
int index_for_target_label(char * name,LabelMaster * lm)
{
  return index_for_label(name,lm->target,lm->target_len);
}


# line 95 "labelmaster.dy"
Label * target_Label_from_name(LabelMaster * lm,char * name)
{
  register int i;

  for(i=0;i<lm->target_len;i++) {
    if( strcmp(name,lm->target[i]->name) == 0 )
      return lm->target[i];
  }

  return NULL;
}

# line 107 "labelmaster.dy"
Label * query_Label_from_name(LabelMaster * lm,char * name)
{
  register int i;

  for(i=0;i<lm->query_len;i++) {
    if( strcmp(name,lm->query[i]->name) == 0 )
      return lm->query[i];
  }

  return NULL;
}


# line 120 "labelmaster.dy"
LabelMaster * LabelMaster_from_GenericMatrix(GenericMatrix * gm)
{
  register int i;
  register int j;
  CellState * state;
  LabelMaster * out;


  out = LabelMaster_alloc_std();


  for(i=0;i<gm->len;i++) {
    state = gm->state[i];

    for(j=0;j<state->len;j++) 
      add_CellSource_to_LabelMaster(out,state->source[j],state->name);
  }

  for(i=0;i<gm->spec_len;i++) {
    state = gm->special[i];

    for(j=0;j<state->len;j++) 
      add_CellSource_to_LabelMaster(out,state->source[j],state->name);
  }

  return out;
}

# line 148 "labelmaster.dy"
boolean add_CellSource_to_LabelMaster(LabelMaster * lm,CellSource * cs,char * state)
{
  Label * temp;
  LabelInstance * li;


  li = LabelInstance_from_CellSource(state,cs);

  if( cs->query_label != NULL ) {
    if( (temp=query_Label_from_name(lm,cs->query_label)) == NULL ) {
      temp = new_query_Label(cs->query_label);
      add_query_LabelMaster(lm,temp);
    }

    add_Label(temp,li);
  }


  li = LabelInstance_from_CellSource(state,cs);

  if( cs->target_label != NULL ) {
    if( (temp=target_Label_from_name(lm,cs->target_label)) == NULL ) {
      temp = new_target_Label(cs->target_label);
      add_target_LabelMaster(lm,temp);
    }

    add_Label(temp,li);
  }

  return TRUE;
}
 

# line 181 "labelmaster.dy"
Label * new_query_Label(char * name)
{
  Label * out;

  out = Label_alloc_std();

  out->name = stringalloc(name);
  out->type = LMQUERYTYPE;

  return out;
}

# line 193 "labelmaster.dy"
Label * new_target_Label(char * name)
{
  Label * out;

  out = Label_alloc_std();

  out->name = stringalloc(name);
  out->type = LMTARGETTYPE;

  return out;
}


# line 206 "labelmaster.dy"
LabelInstance * LabelInstance_from_CellSource(char * state,CellSource * cs)
{
  LabelInstance * out;


  out = LabelInstance_alloc();

  out->state  = stringalloc(state);
  out->source = cs->state_source;
  out->offi = cs->offi;
  out->offj = cs->offj;

  return out;
}







# line 202 "labelmaster.c"
/* Function:  hard_link_LabelInstance(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LabelInstance *]
 *
 * Return [UNKN ]  Undocumented return value [LabelInstance *]
 *
 */
LabelInstance * hard_link_LabelInstance(LabelInstance * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LabelInstance object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LabelInstance_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LabelInstance *]
 *
 */
LabelInstance * LabelInstance_alloc(void) 
{
    LabelInstance * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LabelInstance *) ckalloc (sizeof(LabelInstance))) == NULL)  {  
      warn("LabelInstance_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->state = NULL;   
    out->source = NULL;  
    out->offi = 0;   
    out->offj = 0;   
    out->calc_line = NULL;   


    return out;  
}    


/* Function:  free_LabelInstance(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LabelInstance *]
 *
 * Return [UNKN ]  Undocumented return value [LabelInstance *]
 *
 */
LabelInstance * free_LabelInstance(LabelInstance * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LabelInstance obj. Should be trappable"); 
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
    if( obj->state != NULL)  
      ckfree(obj->state);    
    if( obj->source != NULL) 
      ckfree(obj->source);   
    if( obj->calc_line != NULL)  
      ckfree(obj->calc_line);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_Label(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Label
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [LabelInstance **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Label(LabelInstance ** list,int i,int j)  
{
    LabelInstance * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Label(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Label which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [LabelInstance **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Label(LabelInstance ** list,int left,int right,int (*comp)(LabelInstance * ,LabelInstance * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Label(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Label (list,++last,i);  
      }  
    swap_Label (list,left,last); 
    qsort_Label(list,left,last-1,comp);  
    qsort_Label(list,last+1,right,comp); 
}    


/* Function:  sort_Label(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Label
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Label *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Label(Label * obj,int (*comp)(LabelInstance *, LabelInstance *)) 
{
    qsort_Label(obj->li,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_Label(obj,len)
 *
 * Descrip:    Really an internal function for add_Label
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Label *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Label(Label * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Label called with no need");  
      return TRUE;   
      }  


    if( (obj->li = (LabelInstance ** ) ckrealloc (obj->li,sizeof(LabelInstance *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Label, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Label(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Label *]
 * Arg:        add [OWNER] Object to add to the list [LabelInstance *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Label(Label * obj,LabelInstance * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Label(obj,obj->len + LabelLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->li[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_Label(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Label *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Label(Label * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->li[i] != NULL)    {  
        free_LabelInstance(obj->li[i]);  
        obj->li[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Label_alloc_std(void)
 *
 * Descrip:    Equivalent to Label_alloc_len(LabelLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Label *]
 *
 */
Label * Label_alloc_std(void) 
{
    return Label_alloc_len(LabelLISTLENGTH); 
}    


/* Function:  Label_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Label *]
 *
 */
Label * Label_alloc_len(int len) 
{
    Label * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Label_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->li = (LabelInstance ** ) ckcalloc (len,sizeof(LabelInstance *))) == NULL)   {  
      warn("Warning, ckcalloc failed in Label_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Label(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Label *]
 *
 * Return [UNKN ]  Undocumented return value [Label *]
 *
 */
Label * hard_link_Label(Label * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Label object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Label_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Label *]
 *
 */
Label * Label_alloc(void) 
{
    Label * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Label *) ckalloc (sizeof(Label))) == NULL)  {  
      warn("Label_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->type = 0;   
    out->li = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_Label(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Label *]
 *
 * Return [UNKN ]  Undocumented return value [Label *]
 *
 */
Label * free_Label(Label * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Label obj. Should be trappable"); 
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
    if( obj->li != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->li[i] != NULL)  
          free_LabelInstance(obj->li[i]);    
        }  
      ckfree(obj->li);   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_query_LabelMaster(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_query_LabelMaster
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Label **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_query_LabelMaster(Label ** list,int i,int j)  
{
    Label * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_query_LabelMaster(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_query_LabelMaster which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Label **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_query_LabelMaster(Label ** list,int left,int right,int (*comp)(Label * ,Label * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_query_LabelMaster(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_query_LabelMaster (list,++last,i);  
      }  
    swap_query_LabelMaster (list,left,last); 
    qsort_query_LabelMaster(list,left,last-1,comp);  
    qsort_query_LabelMaster(list,last+1,right,comp); 
}    


/* Function:  sort_query_LabelMaster(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_query_LabelMaster
 *
 *
 * Arg:         obj [UNKN ] Object containing list [LabelMaster *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_query_LabelMaster(LabelMaster * obj,int (*comp)(Label *, Label *)) 
{
    qsort_query_LabelMaster(obj->query,0,obj->query_len-1,comp); 
    return;  
}    


/* Function:  expand_query_LabelMaster(obj,len)
 *
 * Descrip:    Really an internal function for add_query_LabelMaster
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LabelMaster *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_query_LabelMaster(LabelMaster * obj,int len) 
{


    if( obj->query_maxlen > obj->query_len )     {  
      warn("expand_LabelMasterquery_ called with no need");  
      return TRUE;   
      }  


    if( (obj->query = (Label ** ) ckrealloc (obj->query,sizeof(Label *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_LabelMaster, returning FALSE");  
      return FALSE;  
      }  
    obj->query_maxlen = len; 
    return TRUE; 
}    


/* Function:  add_query_LabelMaster(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LabelMaster *]
 * Arg:        add [OWNER] Object to add to the list [Label *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_query_LabelMaster(LabelMaster * obj,Label * add) 
{
    if( obj->query_len >= obj->query_maxlen) {  
      if( expand_query_LabelMaster(obj,obj->query_len + LabelMasterLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->query[obj->query_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_query_LabelMaster(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LabelMaster *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_query_LabelMaster(LabelMaster * obj) 
{
    int i;   


    for(i=0;i<obj->query_len;i++)    { /*for i over list length*/ 
      if( obj->query[i] != NULL) {  
        free_Label(obj->query[i]);   
        obj->query[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->query_len = 0;  
    return i;    
}    


/* Function:  swap_target_LabelMaster(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_target_LabelMaster
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Label **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_target_LabelMaster(Label ** list,int i,int j)  
{
    Label * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_target_LabelMaster(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_target_LabelMaster which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Label **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_target_LabelMaster(Label ** list,int left,int right,int (*comp)(Label * ,Label * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_target_LabelMaster(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_target_LabelMaster (list,++last,i); 
      }  
    swap_target_LabelMaster (list,left,last);    
    qsort_target_LabelMaster(list,left,last-1,comp); 
    qsort_target_LabelMaster(list,last+1,right,comp);    
}    


/* Function:  sort_target_LabelMaster(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_target_LabelMaster
 *
 *
 * Arg:         obj [UNKN ] Object containing list [LabelMaster *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_target_LabelMaster(LabelMaster * obj,int (*comp)(Label *, Label *)) 
{
    qsort_target_LabelMaster(obj->target,0,obj->target_len-1,comp);  
    return;  
}    


/* Function:  expand_target_LabelMaster(obj,len)
 *
 * Descrip:    Really an internal function for add_target_LabelMaster
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LabelMaster *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_target_LabelMaster(LabelMaster * obj,int len) 
{


    if( obj->target_maxlen > obj->target_len )   {  
      warn("expand_LabelMastertarget_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->target = (Label ** ) ckrealloc (obj->target,sizeof(Label *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_LabelMaster, returning FALSE");  
      return FALSE;  
      }  
    obj->target_maxlen = len;    
    return TRUE; 
}    


/* Function:  add_target_LabelMaster(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LabelMaster *]
 * Arg:        add [OWNER] Object to add to the list [Label *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_target_LabelMaster(LabelMaster * obj,Label * add) 
{
    if( obj->target_len >= obj->target_maxlen)   {  
      if( expand_target_LabelMaster(obj,obj->target_len + LabelMasterLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->target[obj->target_len++]=add;  
    return TRUE; 
}    


/* Function:  flush_target_LabelMaster(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LabelMaster *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_target_LabelMaster(LabelMaster * obj) 
{
    int i;   


    for(i=0;i<obj->target_len;i++)   { /*for i over list length*/ 
      if( obj->target[i] != NULL)    {  
        free_Label(obj->target[i]);  
        obj->target[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->target_len = 0; 
    return i;    
}    


/* Function:  LabelMaster_alloc_std(void)
 *
 * Descrip:    Equivalent to LabelMaster_alloc_len(LabelMasterLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LabelMaster *]
 *
 */
LabelMaster * LabelMaster_alloc_std(void) 
{
    return LabelMaster_alloc_len(LabelMasterLISTLENGTH); 
}    


/* Function:  LabelMaster_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LabelMaster *]
 *
 */
LabelMaster * LabelMaster_alloc_len(int len) 
{
    LabelMaster * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = LabelMaster_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->query = (Label ** ) ckcalloc (len,sizeof(Label *))) == NULL)    {  
      warn("Warning, ckcalloc failed in LabelMaster_alloc_len"); 
      return NULL;   
      }  
    out->query_len = 0;  
    out->query_maxlen = len; 


    if((out->target = (Label ** ) ckcalloc (len,sizeof(Label *))) == NULL)   {  
      warn("Warning, ckcalloc failed in LabelMaster_alloc_len"); 
      return NULL;   
      }  
    out->target_len = 0; 
    out->target_maxlen = len;    


    return out;  
}    


/* Function:  hard_link_LabelMaster(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LabelMaster *]
 *
 * Return [UNKN ]  Undocumented return value [LabelMaster *]
 *
 */
LabelMaster * hard_link_LabelMaster(LabelMaster * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LabelMaster object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LabelMaster_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LabelMaster *]
 *
 */
LabelMaster * LabelMaster_alloc(void) 
{
    LabelMaster * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LabelMaster *) ckalloc (sizeof(LabelMaster))) == NULL)  {  
      warn("LabelMaster_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->query = NULL;   
    out->query_len = out->query_maxlen = 0;  
    out->target = NULL;  
    out->target_len = out->target_maxlen = 0;    


    return out;  
}    


/* Function:  free_LabelMaster(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LabelMaster *]
 *
 * Return [UNKN ]  Undocumented return value [LabelMaster *]
 *
 */
LabelMaster * free_LabelMaster(LabelMaster * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LabelMaster obj. Should be trappable");   
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
    if( obj->query != NULL)  {  
      for(i=0;i<obj->query_len;i++)  {  
        if( obj->query[i] != NULL)   
          free_Label(obj->query[i]); 
        }  
      ckfree(obj->query);    
      }  
    if( obj->target != NULL) {  
      for(i=0;i<obj->target_len;i++) {  
        if( obj->target[i] != NULL)  
          free_Label(obj->target[i]);    
        }  
      ckfree(obj->target);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

#ifdef _cplusplus
extern "C" {
#endif
#include "matrixdebug.h"

/* Function:  show_PackAln_Debug(de,pal,ofp)
 *
 * Descrip:    Shows PackAln in debug context
 *
 *
 * Arg:         de [UNKN ] Undocumented argument [DebugMatrix *]
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 63 "matrixdebug.dy"
void show_PackAln_Debug(DebugMatrix * de,PackAln * pal,FILE * ofp)
{
  int i;
  int j,k;
  int cum;
  PackAlnUnit * prev = pal->pau[0];
  PackAlnUnit * c;
  DebugState * state;
  DebugTransition * trans;
 
  for(i=1;i<pal->len;i++) {
    c= pal->pau[i];
    /* have to find the correct state_num */
    if( c->state >= de->len ) {
      /* special state - can't handle */
      fprintf(ofp,"Special state\n");
      prev = c;
      continue;
    }
    if( prev->state >= de->len ) {
      /* special state - can't handle */
      fprintf(ofp,"From Special state\n");
      prev = c;
      continue;
    }
    
    state = de->state[pal->pau[i]->state];
    trans = NULL;
    for(k=0;k<state->len;k++) {
      if( state->trans[k]->from_state_num == prev->state ) {
	trans = state->trans[k];
      }
    }

    fprintf(ofp,"SCORE: %5d %5d %5d %s\n",c->score,c->i,c->j,state->statename);

    if( trans != NULL ) {
      (*trans->show_transition)(de->matrix,c->i,c->j,ofp);
    }
    prev = c;
  }
}

/* Function:  show_DebugMatrix(de,buffer)
 *
 * Descrip:    Show function cell state transition
 *
 *
 * Arg:            de [UNKN ] Undocumented argument [DebugMatrix *]
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 *
 */
# line 109 "matrixdebug.dy"
void show_DebugMatrix(DebugMatrix * de,char * buffer)
{
  char ** base;
  char ** brk;
  DebugState * ds;
  DebugTransition * tr;

  base = brk = breakstring(buffer,spacestr);

  brk++;

  if( *brk == NULL ) {
    (*de->show_cell)(de->matrix,de->currenti,de->currentj,de->out);
  } else {
    ds = find_DebugState(de,*brk);
    if( ds == NULL ) {
      fprintf(de->out,"No state called %s\n",*brk);
    } else {
      brk++;
      if( *brk != NULL ) {
	tr = find_DebugTransition(ds,*brk);
	if( tr == NULL ) {
	  fprintf(de->out,"In state %s, no transition from %s\n",ds->statename,*brk);
	} else {
	  (*tr->show_transition)(de->matrix,de->currenti,de->currentj,de->out);
	}
      } else {
	/* show state */
	(*ds->show_state)(de->matrix,de->currenti,de->currentj,de->out);
      }
    }
  }

  ckfree(base);

}


/* Function:  find_DebugState(de,name)
 *
 * Descrip:    Finds a DebugState 
 *
 *
 * Arg:          de [UNKN ] Undocumented argument [DebugMatrix *]
 * Arg:        name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
# line 150 "matrixdebug.dy"
DebugState * find_DebugState(DebugMatrix * de,char * name)
{
  int i;
  
  for(i=0;i<de->len;i++) {
    if( strcmp(de->state[i]->statename,name) == 0 ) {
      return de->state[i];
    } 
  }

  return NULL;
}
 
/* Function:  find_DebugTransition(state,name)
 *
 * Descrip:    Finds a DebugTransition
 *
 *
 * Arg:        state [UNKN ] Undocumented argument [DebugState *]
 * Arg:         name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [DebugTransition *]
 *
 */
# line 166 "matrixdebug.dy"
DebugTransition * find_DebugTransition(DebugState * state,char * name)
{
  int i;

  for(i=0;i<state->len;i++) {
    if( strcmp(state->trans[i]->fromstate,name) == 0 ) {
      return state->trans[i];
    }
  }

  return NULL;
}


/* Function:  user_DebugMatrix(de)
 *
 * Descrip:    Main function to talk to user
 *
 *
 * Arg:        de [UNKN ] Undocumented argument [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 183 "matrixdebug.dy"
boolean user_DebugMatrix(DebugMatrix * de)
{
  char buffer[MAXLINE];
  char ** base;
  char ** brk;
  FILE * in;
  FILE * out;

  DebugBreakPoint * bp;

  assert(de);
  assert(de->in);
  assert(de->out);

  in  = de->in;
  out = de->out;

  /* set reset to FALSE */
  de->reset = 0;

  fprintf(out,"Entering dynamite debugger. Type help for help\n");

  while(1) {
    fprintf(out,"Dy %5d:%5d max: %d >",de->currenti,de->currentj,de->max_score);
    
    fgets(buffer,MAXLINE,in);
    
    if( strstartcmp(buffer,"quit") == 0 ) {
      exit(1);
    }

    if( strstartcmp(buffer,"show") == 0 ) {
      show_DebugMatrix(de,buffer);
      continue;
    }
    
    if( strstartcmp(buffer,"run") == 0 ) {
      /* return to calling function */
      return 0;
    }

    if( strstartcmp(buffer,"break") == 0 ) {
      base = brk = breakstring(buffer,spacestr);
      brk++;
      
      if( *brk == NULL || *(brk+1) == NULL ) {
	fprintf(out,">>>> break must have i j positions\n");
	continue;
      } else {
	bp = DebugBreakPoint_alloc();
	bp->type = MDBP_Cursor;
	/* reset cursor positions */
	if( is_integer_string(*brk,&bp->posi) == FALSE ) {
	  fprintf(out,">>>> i position not an integer.\n");
	} 
	brk++;
	if( is_integer_string(*brk,&bp->posj) == FALSE ) {
	  fprintf(out,">>>> j position not an integer.\n");
	} 	
	fprintf(out,"Adding cursor break point %d,%d\n",bp->posi,bp->posj);
	add_bp_DebugMatrix(de,bp);
      }
      ckfree(base);
      continue;
    }

    if( strstartcmp(buffer,"set") == 0 ) {
      base = brk = breakstring(buffer,spacestr);
      brk++;
      if( *brk == NULL || *(brk+1) == NULL ) {
	fprintf(out,">>>> set must have i j positions\n");
      } else {
	/* reset cursor positions */
	if( is_integer_string(*brk,&de->currenti) == FALSE ) {
	  fprintf(out,">>>> i position not an integer.\n");
	} 
	if( is_integer_string(*brk,&de->currenti) == FALSE ) {
	  fprintf(out,">>>> j position not an integer.\n");
	} 	
      }
      de->reset = 1;
      ckfree(base);
      continue;
    }

  }

}


/* Function:  should_break_DebugMatrix(de)
 *
 * Descrip:    Indicates whether we should break at this point. Assummes
 *             the de datastructure has been updated correctly
 *
 *
 * Arg:        de [UNKN ] Undocumented argument [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [MatrixDebugBreakPoint_type]
 *
 */
# line 277 "matrixdebug.dy"
MatrixDebugBreakPoint_type should_break_DebugMatrix(DebugMatrix * de)
{
  int i;

  for(i=0;i<de->bp_len;i++) {
    auto DebugBreakPoint * p;
    p = de->point[i];

    switch(p->type) {
    case MDBP_Cursor :
      if( de->currenti == p->posi && de->currentj == p->posj ) {
	return MDBP_Cursor;
      } 
      break;
    case MDBP_Overflow :
      if( de->max_score > p->overflow ) {
	return MDBP_Overflow;
      }
      break;
    case MDBP_Underflow :
      if( de->min_score < p->overflow ) {
	return MDBP_Overflow;
      }
      break;
    default :
      warn("Weird break point type %d",p->type);
    }

  }

  return MDBP_NoBreakPoint;
      
}


/* Function:  std_DebugMatrix(void)
 *
 * Descrip:    Builds a "standard" Debug matrix
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
# line 315 "matrixdebug.dy"
DebugMatrix * std_DebugMatrix(void)
{
  DebugMatrix * out;
  DebugBreakPoint * bp;

  out = DebugMatrix_alloc_std();

  bp = DebugBreakPoint_alloc();

  bp->type = MDBP_Cursor;
  bp->posi = -1;
  bp->posj = -1;

  add_bp_DebugMatrix(out,bp);

  bp = DebugBreakPoint_alloc();

  bp->type = MDBP_Overflow;
  bp->overflow = 1000000;

  add_bp_DebugMatrix(out,bp);

  out->in = stdin;
  out->out = stderr;
  return out;
}

# line 329 "matrixdebug.c"
/* Function:  hard_link_DebugTransition(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DebugTransition *]
 *
 * Return [UNKN ]  Undocumented return value [DebugTransition *]
 *
 */
DebugTransition * hard_link_DebugTransition(DebugTransition * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DebugTransition object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DebugTransition_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugTransition *]
 *
 */
DebugTransition * DebugTransition_alloc(void) 
{
    DebugTransition * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DebugTransition *) ckalloc (sizeof(DebugTransition))) == NULL)  {  
      warn("DebugTransition_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->from_state_num = 0; 
    out->show_transition = NULL; 


    return out;  
}    


/* Function:  free_DebugTransition(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DebugTransition *]
 *
 * Return [UNKN ]  Undocumented return value [DebugTransition *]
 *
 */
DebugTransition * free_DebugTransition(DebugTransition * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DebugTransition obj. Should be trappable");   
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
    /* obj->fromstate is linked in */ 
    /* obj->show_transition is a function pointer */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_DebugState(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_DebugState
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DebugTransition **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_DebugState(DebugTransition ** list,int i,int j)  
{
    DebugTransition * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_DebugState(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_DebugState which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DebugTransition **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_DebugState(DebugTransition ** list,int left,int right,int (*comp)(DebugTransition * ,DebugTransition * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_DebugState(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_DebugState (list,++last,i); 
      }  
    swap_DebugState (list,left,last);    
    qsort_DebugState(list,left,last-1,comp); 
    qsort_DebugState(list,last+1,right,comp);    
}    


/* Function:  sort_DebugState(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_DebugState
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DebugState *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_DebugState(DebugState * obj,int (*comp)(DebugTransition *, DebugTransition *)) 
{
    qsort_DebugState(obj->trans,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_DebugState(obj,len)
 *
 * Descrip:    Really an internal function for add_DebugState
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DebugState *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DebugState(DebugState * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_DebugState called with no need"); 
      return TRUE;   
      }  


    if( (obj->trans = (DebugTransition ** ) ckrealloc (obj->trans,sizeof(DebugTransition *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_DebugState, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_DebugState(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DebugState *]
 * Arg:        add [OWNER] Object to add to the list [DebugTransition *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_DebugState(DebugState * obj,DebugTransition * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_DebugState(obj,obj->len + DebugStateLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->trans[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_DebugState(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DebugState *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DebugState(DebugState * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->trans[i] != NULL) {  
        free_DebugTransition(obj->trans[i]); 
        obj->trans[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  DebugState_alloc_std(void)
 *
 * Descrip:    Equivalent to DebugState_alloc_len(DebugStateLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
DebugState * DebugState_alloc_std(void) 
{
    return DebugState_alloc_len(DebugStateLISTLENGTH);   
}    


/* Function:  DebugState_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
DebugState * DebugState_alloc_len(int len) 
{
    DebugState * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = DebugState_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->trans = (DebugTransition ** ) ckcalloc (len,sizeof(DebugTransition *))) == NULL)    {  
      warn("Warning, ckcalloc failed in DebugState_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_DebugState(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DebugState *]
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
DebugState * hard_link_DebugState(DebugState * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DebugState object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DebugState_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
DebugState * DebugState_alloc(void) 
{
    DebugState * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DebugState *) ckalloc (sizeof(DebugState))) == NULL)    {  
      warn("DebugState_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->state_num = 0;  
    out->trans = NULL;   
    out->len = out->maxlen = 0;  
    out->show_state = NULL;  


    return out;  
}    


/* Function:  free_DebugState(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DebugState *]
 *
 * Return [UNKN ]  Undocumented return value [DebugState *]
 *
 */
DebugState * free_DebugState(DebugState * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DebugState obj. Should be trappable");    
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
    /* obj->statename is linked in */ 
    if( obj->trans != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->trans[i] != NULL)   
          free_DebugTransition(obj->trans[i]);   
        }  
      ckfree(obj->trans);    
      }  
    /* obj->show_state is a function pointer */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_DebugBreakPoint(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DebugBreakPoint *]
 *
 * Return [UNKN ]  Undocumented return value [DebugBreakPoint *]
 *
 */
DebugBreakPoint * hard_link_DebugBreakPoint(DebugBreakPoint * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DebugBreakPoint object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DebugBreakPoint_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugBreakPoint *]
 *
 */
DebugBreakPoint * DebugBreakPoint_alloc(void) 
{
    DebugBreakPoint * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DebugBreakPoint *) ckalloc (sizeof(DebugBreakPoint))) == NULL)  {  
      warn("DebugBreakPoint_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = MDBP_Cursor; 
    out->posi = 0;   
    out->posj = 0;   
    out->overflow = 0;   
    out->underflow = 0;  


    return out;  
}    


/* Function:  free_DebugBreakPoint(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DebugBreakPoint *]
 *
 * Return [UNKN ]  Undocumented return value [DebugBreakPoint *]
 *
 */
DebugBreakPoint * free_DebugBreakPoint(DebugBreakPoint * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DebugBreakPoint obj. Should be trappable");   
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


/* Function:  swap_DebugMatrix(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_DebugMatrix
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DebugState **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_DebugMatrix(DebugState ** list,int i,int j)  
{
    DebugState * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_DebugMatrix(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_DebugMatrix which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DebugState **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_DebugMatrix(DebugState ** list,int left,int right,int (*comp)(DebugState * ,DebugState * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_DebugMatrix(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_DebugMatrix (list,++last,i);    
      }  
    swap_DebugMatrix (list,left,last);   
    qsort_DebugMatrix(list,left,last-1,comp);    
    qsort_DebugMatrix(list,last+1,right,comp);   
}    


/* Function:  sort_DebugMatrix(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_DebugMatrix
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DebugMatrix *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_DebugMatrix(DebugMatrix * obj,int (*comp)(DebugState *, DebugState *)) 
{
    qsort_DebugMatrix(obj->state,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_DebugMatrix(obj,len)
 *
 * Descrip:    Really an internal function for add_DebugMatrix
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DebugMatrix *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DebugMatrix(DebugMatrix * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_DebugMatrix called with no need");    
      return TRUE;   
      }  


    if( (obj->state = (DebugState ** ) ckrealloc (obj->state,sizeof(DebugState *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_DebugMatrix, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_DebugMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DebugMatrix *]
 * Arg:        add [OWNER] Object to add to the list [DebugState *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_DebugMatrix(DebugMatrix * obj,DebugState * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_DebugMatrix(obj,obj->len + DebugMatrixLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->state[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_DebugMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DebugMatrix(DebugMatrix * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->state[i] != NULL) {  
        free_DebugState(obj->state[i]);  
        obj->state[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_bp_DebugMatrix(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_bp_DebugMatrix
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DebugBreakPoint **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_bp_DebugMatrix(DebugBreakPoint ** list,int i,int j)  
{
    DebugBreakPoint * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_bp_DebugMatrix(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_bp_DebugMatrix which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DebugBreakPoint **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_bp_DebugMatrix(DebugBreakPoint ** list,int left,int right,int (*comp)(DebugBreakPoint * ,DebugBreakPoint * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_bp_DebugMatrix(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_bp_DebugMatrix (list,++last,i); 
      }  
    swap_bp_DebugMatrix (list,left,last);    
    qsort_bp_DebugMatrix(list,left,last-1,comp); 
    qsort_bp_DebugMatrix(list,last+1,right,comp);    
}    


/* Function:  sort_bp_DebugMatrix(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_bp_DebugMatrix
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DebugMatrix *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_bp_DebugMatrix(DebugMatrix * obj,int (*comp)(DebugBreakPoint *, DebugBreakPoint *)) 
{
    qsort_bp_DebugMatrix(obj->point,0,obj->bp_len-1,comp);   
    return;  
}    


/* Function:  expand_bp_DebugMatrix(obj,len)
 *
 * Descrip:    Really an internal function for add_bp_DebugMatrix
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DebugMatrix *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_bp_DebugMatrix(DebugMatrix * obj,int len) 
{


    if( obj->bp_maxlen > obj->bp_len )   {  
      warn("expand_DebugMatrixbp_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->point = (DebugBreakPoint ** ) ckrealloc (obj->point,sizeof(DebugBreakPoint *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_DebugMatrix, returning FALSE");  
      return FALSE;  
      }  
    obj->bp_maxlen = len;    
    return TRUE; 
}    


/* Function:  add_bp_DebugMatrix(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DebugMatrix *]
 * Arg:        add [OWNER] Object to add to the list [DebugBreakPoint *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_bp_DebugMatrix(DebugMatrix * obj,DebugBreakPoint * add) 
{
    if( obj->bp_len >= obj->bp_maxlen)   {  
      if( expand_bp_DebugMatrix(obj,obj->bp_len + DebugMatrixLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->point[obj->bp_len++]=add;   
    return TRUE; 
}    


/* Function:  flush_bp_DebugMatrix(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_bp_DebugMatrix(DebugMatrix * obj) 
{
    int i;   


    for(i=0;i<obj->bp_len;i++)   { /*for i over list length*/ 
      if( obj->point[i] != NULL) {  
        free_DebugBreakPoint(obj->point[i]); 
        obj->point[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->bp_len = 0; 
    return i;    
}    


/* Function:  DebugMatrix_alloc_std(void)
 *
 * Descrip:    Equivalent to DebugMatrix_alloc_len(DebugMatrixLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
DebugMatrix * DebugMatrix_alloc_std(void) 
{
    return DebugMatrix_alloc_len(DebugMatrixLISTLENGTH); 
}    


/* Function:  DebugMatrix_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
DebugMatrix * DebugMatrix_alloc_len(int len) 
{
    DebugMatrix * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = DebugMatrix_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->state = (DebugState ** ) ckcalloc (len,sizeof(DebugState *))) == NULL)  {  
      warn("Warning, ckcalloc failed in DebugMatrix_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->point = (DebugBreakPoint ** ) ckcalloc (len,sizeof(DebugBreakPoint *))) == NULL)    {  
      warn("Warning, ckcalloc failed in DebugMatrix_alloc_len"); 
      return NULL;   
      }  
    out->bp_len = 0; 
    out->bp_maxlen = len;    


    return out;  
}    


/* Function:  hard_link_DebugMatrix(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
DebugMatrix * hard_link_DebugMatrix(DebugMatrix * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DebugMatrix object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DebugMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
DebugMatrix * DebugMatrix_alloc(void) 
{
    DebugMatrix * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DebugMatrix *) ckalloc (sizeof(DebugMatrix))) == NULL)  {  
      warn("DebugMatrix_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->currenti = -1;  
    out->currentj = -1;  
    out->max_score = 0;  
    out->min_score = 100000; 
    out->max_score_i = -1;   
    out->max_score_j = -1;   
    out->max_score_cell = 0; 
    out->state = NULL;   
    out->len = out->maxlen = 0;  
    out->point = NULL;   
    out->bp_len = out->bp_maxlen = 0;    
    out->reset = FALSE;  
    out->show_cell = NULL;   


    return out;  
}    


/* Function:  free_DebugMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DebugMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [DebugMatrix *]
 *
 */
DebugMatrix * free_DebugMatrix(DebugMatrix * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DebugMatrix obj. Should be trappable");   
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
    /* obj->matrix is linked in */ 
    if( obj->state != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->state[i] != NULL)   
          free_DebugState(obj->state[i]);    
        }  
      ckfree(obj->state);    
      }  
    if( obj->point != NULL)  {  
      for(i=0;i<obj->bp_len;i++) {  
        if( obj->point[i] != NULL)   
          free_DebugBreakPoint(obj->point[i]);   
        }  
      ckfree(obj->point);    
      }  
    /* obj->in is linked in */ 
    /* obj->out is linked in */ 
    /* obj->show_cell is a function pointer */ 


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

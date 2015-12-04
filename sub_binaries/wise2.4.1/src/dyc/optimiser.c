#ifdef _cplusplus
extern "C" {
#endif
#include "optimiser.h"

/* Function:  expr_ij_dependence(expr,ijdep)
 *
 * Descrip:    figures out whether this cse is i or j dependent
 *
 *
 * Arg:         expr [UNKN ] Undocumented argument [ExprTree *]
 * Arg:        ijdep [UNKN ] Undocumented argument [int *]
 *
 */
# line 33 "optimiser.dy"
void expr_ij_dependence(ExprTree * expr,int * ijdep)
{
  int i;

  if( expr->type == ETR_NAME ) {
    if( strcmp(expr->word,"i") == 0 ) {
      *ijdep = (*ijdep | CSE_I_DEP);
    } 
    if( strcmp(expr->word,"j") == 0 ) {
      *ijdep = (*ijdep | CSE_J_DEP);
    } 
  }

  for(i=0;i<expr->nochild;i++) {
    expr_ij_dependence(expr->child[i],ijdep);
  }
}

/* Function:  strcat_cses_ExprTree(epr,buffer,sc,mts,dpi)
 *
 * Descrip:    Writes code with sub expressions in the correct places
 *
 *
 * Arg:           epr [UNKN ] Undocumented argument [ExprTree *]
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 * Arg:            sc [UNKN ] Undocumented argument [Scope *]
 * Arg:           mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:           dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
# line 54 "optimiser.dy"
void strcat_cses_ExprTree(ExprTree * epr,char * buffer,Scope * sc,MethodTypeSet * mts,DPImplementation * dpi)
{
  strcat_ExprTree_Scoped(epr,buffer,sc,mts,dpi->dycw,cses_expr_placer,NULL);
}

/* Function:  cses_expr_placer(etr,buffer,data)
 *
 * Descrip:    pointer to function that does the magic on the cses placer system
 *
 *
 * Arg:           etr [UNKN ] Undocumented argument [ExprTree *]
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 * Arg:          data [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 63 "optimiser.dy"
boolean cses_expr_placer(ExprTree * etr,char * buffer,void * data)
{
  char b[20];
  /*  fprintf(stdout,"Getting into scoping, but with %d and %d\n",etr->type,etr->cse);*/
  if( etr->cse == NULL ) 
    return FALSE;

  strcat(buffer,"subexpr");
  sprintf(b,"%d",etr->cse->id);
  strcat(buffer,b);

  return TRUE;
}

/* Function:  id_CommonSubExpressionSet(cses)
 *
 * Descrip:    Places ids into the cses
 *
 *
 * Arg:        cses [UNKN ] Undocumented argument [CommonSubExpressionSet *]
 *
 */
# line 81 "optimiser.dy"
void id_CommonSubExpressionSet(CommonSubExpressionSet * cses)
{
  int i;

  for(i=0;i<cses->len;i++) 
    cses->cse[i]->id = i;

}




/* Function:  show_CommonSubExpressionSet(cses,ofp)
 *
 * Descrip:    Shows common sub expression set
 *
 *
 * Arg:        cses [UNKN ] Undocumented argument [CommonSubExpressionSet *]
 * Arg:         ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 96 "optimiser.dy"
void show_CommonSubExpressionSet(CommonSubExpressionSet * cses,FILE * ofp)
{
  int i;

  for(i=0;i<cses->len;i++) {
    fprintf(ofp,"Sub Expression %d - Number of times seen %d IJdep %d\n   ",i,cses->cse[i]->number,cses->cse[i]->type);
    print_ExprTree(cses->cse[i]->expr,ofp);
    fprintf(ofp,"\n");
  }
}


/* Function:  find_CommonSubExpressions(gm,source_ind_promote)
 *
 * Descrip:    Makes a CommonSubExpressionSet from the
 *             GenericMatrix structure. 
 *
 *             if source_ind_promote is true, then source_independent type systems
 *             are automatically promoted as a common sub expression, regardless of 
 *             number of subexpressions found
 *
 *
 * Arg:                        gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:        source_ind_promote [UNKN ] Undocumented argument [boolean]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
# line 116 "optimiser.dy"
CommonSubExpressionSet * find_CommonSubExpressions(GenericMatrix * gm,boolean source_ind_promote)
{
  int i,j;
  ExprTree * exlist[512];
  int notags =0;

  CommonSubExpressionSet * out;
  CommonSubExpression * temp;

  
  out = CommonSubExpressionSet_alloc_std();

  for(i=0;i<gm->len;i++) {
    for(j=0;j<gm->state[i]->len;j++) {
      attach_expressions_to_list(exlist,&notags,gm->state[i]->source[j]->etr);
    }
    if( source_ind_promote == TRUE && gm->state[i]->etr != NULL) {
      temp = CommonSubExpression_alloc();
      temp->expr = gm->state[i]->etr;
      gm->state[i]->etr->cse = temp;
      temp->number = 1;
      add_CommonSubExpressionSet(out,temp);
    } else {
      attach_expressions_to_list(exlist,&notags,gm->state[i]->etr);
    }
  }

  for(i=0;i<gm->spec_len;i++) {
    for(j=0;j<gm->special[i]->len;j++) {
      attach_expressions_to_list(exlist,&notags,gm->special[i]->source[j]->etr);
    }
    attach_expressions_to_list(exlist,&notags,gm->special[i]->etr);
  }

  for(i=0;i<notags;i++) {
    /*   fprintf(stderr,"Looking at %d  ",i);
    print_ExprTree(exlist[i],stderr);
    fprintf(stderr,"\n");
    */

    if( exlist[i]->cse != NULL )
      continue;
    temp = NULL;
    for(j=i+1;j<notags;j++) {
      if( exlist[j]->cse != NULL ) 
	continue;
      
      if( identical_ExprTree(exlist[i],exlist[j]) == TRUE ) {
	if( temp == NULL ) {
	  temp = CommonSubExpression_alloc();
	  temp->expr = exlist[i];
	  exlist[i]->cse = temp;
	  temp->number = 1;
	  add_CommonSubExpressionSet(out,temp);
	}
	
	/* we do this for every subexpression that we find */
	exlist[j]->cse = temp;
	temp->number++;
      }
    }
  }

  /* find i,j dependence */
  for(i=0;i<out->len;i++) {
    expr_ij_dependence(out->cse[i]->expr,&out->cse[i]->type);
  }

  id_CommonSubExpressionSet(out);
  return out;

}

/* Function:  should_store_ExprTree_cse(start)
 *
 * Descrip:    whether we should consider this a possible common sub expression or not
 *
 *
 * Arg:        start [UNKN ] Undocumented argument [ExprTree *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 193 "optimiser.dy"
boolean should_store_ExprTree_cse(ExprTree * start)
{
  if( start->type == ETR_EXPRESSION ||  start->type == ETR_METHOD )
    return TRUE;

  if( start->type == ETR_TAG || start->type == ETR_ARRAY || start->type == ETR_STRUCTREF || start->type == ETR_REFERENCE ) {
    if( start->parent->type == ETR_EXPRESSION || start->parent->type == ETR_STATEMENT) {
      return TRUE;
    }
  }

  return FALSE;
}


/* Function:  attach_expressions_to_list(list,no,start)
 *
 * Descrip:    Attaches only ETR_EXPRESSIONS into the list.
 *             Really a subroutine for find_CommonSubExpressions
 *
 *
 * Arg:         list [UNKN ] Undocumented argument [ExprTree **]
 * Arg:           no [UNKN ] Undocumented argument [int *]
 * Arg:        start [UNKN ] Undocumented argument [ExprTree *]
 *
 */
# line 212 "optimiser.dy"
void attach_expressions_to_list(ExprTree ** list,int * no,ExprTree * start)
{
  int i;
  if( start == NULL )
    return;

  /* fprintf(stderr,"Looking at expression %d... to %d\n",start->type,ETR_EXPRESSION); */
  if( should_store_ExprTree_cse(start) == TRUE ) {
    list[*no] = start;
    *no = *no +1;
  }

  if( start->type == ETR_METHOD && start->nochild == 2) {
    for(i=0;i<start->child[1]->nochild;i++) 
      attach_expressions_to_list(list,no,start->child[1]->child[i]);
  } else {
    for(i=0;i<start->nochild;i++) 
      attach_expressions_to_list(list,no,start->child[i]);
  }
}
  

/* Function:  identical_ExprTree(one,two)
 *
 * Descrip:    Says whether two ExprTrees are completely identical - ie,
 *             same trees with same tags in the same order. 
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [ExprTree *]
 * Arg:        two [UNKN ] Undocumented argument [ExprTree *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 239 "optimiser.dy"
boolean identical_ExprTree(ExprTree * one,ExprTree * two)
{
  int i;
  if( one->type != two->type ) 
    return FALSE;

  if( one->word != NULL || two->word != NULL ) {
    if( strcmp(one->word,two->word) != 0 )
      return FALSE;
  }

  if( one->nochild != two->nochild )
    return FALSE;

  for(i=0;i<one->nochild;i++ ) {
    if( identical_ExprTree(one->child[i],two->child[i]) == FALSE ) {
      return FALSE;
    }
  }

  return TRUE;
}
  

# line 290 "optimiser.c"
/* Function:  hard_link_CommonSubExpression(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CommonSubExpression *]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpression *]
 *
 */
CommonSubExpression * hard_link_CommonSubExpression(CommonSubExpression * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CommonSubExpression object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CommonSubExpression_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpression *]
 *
 */
CommonSubExpression * CommonSubExpression_alloc(void) 
{
    CommonSubExpression * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CommonSubExpression *) ckalloc (sizeof(CommonSubExpression))) == NULL)  {  
      warn("CommonSubExpression_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->id = 0; 
    out->type = 0;   
    out->number = 0; 


    return out;  
}    


/* Function:  free_CommonSubExpression(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CommonSubExpression *]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpression *]
 *
 */
CommonSubExpression * free_CommonSubExpression(CommonSubExpression * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CommonSubExpression obj. Should be trappable");   
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
    /* obj->expr is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_CommonSubExpressionSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_CommonSubExpressionSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [CommonSubExpression **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_CommonSubExpressionSet(CommonSubExpression ** list,int i,int j)  
{
    CommonSubExpression * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_CommonSubExpressionSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_CommonSubExpressionSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [CommonSubExpression **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_CommonSubExpressionSet(CommonSubExpression ** list,int left,int right,int (*comp)(CommonSubExpression * ,CommonSubExpression * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_CommonSubExpressionSet(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_CommonSubExpressionSet (list,++last,i); 
      }  
    swap_CommonSubExpressionSet (list,left,last);    
    qsort_CommonSubExpressionSet(list,left,last-1,comp); 
    qsort_CommonSubExpressionSet(list,last+1,right,comp);    
}    


/* Function:  sort_CommonSubExpressionSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_CommonSubExpressionSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [CommonSubExpressionSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_CommonSubExpressionSet(CommonSubExpressionSet * obj,int (*comp)(CommonSubExpression *, CommonSubExpression *)) 
{
    qsort_CommonSubExpressionSet(obj->cse,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_CommonSubExpressionSet(obj,len)
 *
 * Descrip:    Really an internal function for add_CommonSubExpressionSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [CommonSubExpressionSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_CommonSubExpressionSet(CommonSubExpressionSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_CommonSubExpressionSet called with no need"); 
      return TRUE;   
      }  


    if( (obj->cse = (CommonSubExpression ** ) ckrealloc (obj->cse,sizeof(CommonSubExpression *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_CommonSubExpressionSet, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_CommonSubExpressionSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [CommonSubExpressionSet *]
 * Arg:        add [OWNER] Object to add to the list [CommonSubExpression *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_CommonSubExpressionSet(CommonSubExpressionSet * obj,CommonSubExpression * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_CommonSubExpressionSet(obj,obj->len + CommonSubExpressionSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->cse[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_CommonSubExpressionSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [CommonSubExpressionSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_CommonSubExpressionSet(CommonSubExpressionSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->cse[i] != NULL)   {  
        free_CommonSubExpression(obj->cse[i]);   
        obj->cse[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  CommonSubExpressionSet_alloc_std(void)
 *
 * Descrip:    Equivalent to CommonSubExpressionSet_alloc_len(CommonSubExpressionSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
CommonSubExpressionSet * CommonSubExpressionSet_alloc_std(void) 
{
    return CommonSubExpressionSet_alloc_len(CommonSubExpressionSetLISTLENGTH);   
}    


/* Function:  CommonSubExpressionSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
CommonSubExpressionSet * CommonSubExpressionSet_alloc_len(int len) 
{
    CommonSubExpressionSet * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = CommonSubExpressionSet_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->cse = (CommonSubExpression ** ) ckcalloc (len,sizeof(CommonSubExpression *))) == NULL)  {  
      warn("Warning, ckcalloc failed in CommonSubExpressionSet_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_CommonSubExpressionSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [CommonSubExpressionSet *]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
CommonSubExpressionSet * hard_link_CommonSubExpressionSet(CommonSubExpressionSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a CommonSubExpressionSet object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  CommonSubExpressionSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
CommonSubExpressionSet * CommonSubExpressionSet_alloc(void) 
{
    CommonSubExpressionSet * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(CommonSubExpressionSet *) ckalloc (sizeof(CommonSubExpressionSet))) == NULL)    {  
      warn("CommonSubExpressionSet_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->cse = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_CommonSubExpressionSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [CommonSubExpressionSet *]
 *
 * Return [UNKN ]  Undocumented return value [CommonSubExpressionSet *]
 *
 */
CommonSubExpressionSet * free_CommonSubExpressionSet(CommonSubExpressionSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a CommonSubExpressionSet obj. Should be trappable");    
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
    if( obj->cse != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->cse[i] != NULL) 
          free_CommonSubExpression(obj->cse[i]); 
        }  
      ckfree(obj->cse);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

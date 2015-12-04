#ifdef _cplusplus
extern "C" {
#endif
#include "dynadebug.h"

# line 20 "dynadebug.dy"
void write_debug_funcs(DYNFILE * dfp,GenericMatrix * gm)
{
  
  explicit_debug_func(dfp,gm);

  transition_debug_func(dfp,gm);

  make_debug_struct_func(dfp,gm);

  state_debug_func(dfp,gm);
}


# line 33 "dynadebug.dy"
ExprSet * build_ExprSet(ExprTree * root)
{
  ExprSet * out;


  assert(root);
  assert(root->type == ETR_STATEMENT);
  assert(root->child[0]);

  

  out = ExprSet_alloc_std();

  /* deal with simple cases: number,tag, method */
  
  if( root->child[0]->type != ETR_EXPRESSION ) {
    /* only one child */
    add_ExprSet(out,hard_link_ExprTree(root->child[0]));
  } else {
    if( descend_ExprTree_for_ExprSet(out,root->child[0]) == FALSE ) {
      free_ExprSet(out);
      return NULL;
    }
  }

  return out;
}

# line 61 "dynadebug.dy"
boolean descend_ExprTree_for_ExprSet(ExprSet * set,ExprTree * t)
{
  boolean ret = TRUE;

  assert(set);
  assert(t);

  if( t->type != ETR_EXPRESSION ) {
    warn("In descending ExprTree for making an ExprSet, failed!");
    return FALSE;
  }

  if( t->child[0]->type == ETR_EXPRESSION ) {
    if( descend_ExprTree_for_ExprSet(set,t->child[0]) == FALSE ) {
      ret = FALSE;
    }
  } else {
    add_ExprSet(set,hard_link_ExprTree(t->child[0]));
  }

  if( t->child[2]->type == ETR_EXPRESSION ) {
    if( descend_ExprTree_for_ExprSet(set,t->child[2]) == FALSE ) {
      ret = FALSE;
    }
  } else {
    add_ExprSet(set,hard_link_ExprTree(t->child[2]));
  }

  return ret;
}

# line 92 "dynadebug.dy"
void make_debug_struct_func(DYNFILE * dfp,GenericMatrix * gm)
{
  FuncInfo * fi;
  int k,l,m;

  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"create_debug_%s",gm->name);
  start_function_FuncInfo(fi,dfp,"DebugMatrix * create_debug_%s(void)",gm->name);
  
  expr(dfp,"DebugMatrix * out;");
  expr(dfp,"DebugState * st;");
  expr(dfp,"DebugTransition * tr;");

  add_break(dfp);

  expr(dfp,"out = std_DebugMatrix();");
  expr(dfp,"out->show_cell = cell_debug_%s;",gm->name);
  add_break(dfp);

  for(k=0;k<gm->len;k++) {
    expr(dfp,"st = DebugState_alloc_std();");
    expr(dfp,"st->statename = \"%s\";",gm->state[k]->name);
    expr(dfp,"st->state_num = %d",k);
    expr(dfp,"st->show_state = state_debug_%s_%s;",gm->name,gm->state[k]->name);
    expr(dfp,"add_DebugMatrix(out,st);");
    
    for(l=0;l<gm->state[k]->len;l++) {
      expr(dfp,"tr = DebugTransition_alloc();");
      expr(dfp,"tr->fromstate = \"%s\";",gm->state[k]->source[l]->state_source);
      expr(dfp,"tr->from_state_num = %d;",gm->state[k]->source[l]->from_state_no); 
      expr(dfp,"tr->show_transition = transition_debug_%s_%s_%s_%d_%d;",gm->name,gm->state[k]->name,gm->state[k]->source[l]->state_source,gm->state[k]->source[l]->offi,gm->state[k]->source[l]->offj);
      expr(dfp,"add_DebugState(st,tr);");
    }
  }

  expr(dfp,"return out;");
  close_function(dfp);
  add_break(dfp);

}


# line 133 "dynadebug.dy"
void state_debug_func(DYNFILE * dfp,GenericMatrix * gm)
{
  FuncInfo * fi;
  int k,l;

  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"cell_debug_%s",gm->name);
  start_function_FuncInfo(fi,dfp,"void cell_debug_%s(void * matrixpoi,int i,int j,FILE * ofp)",gm->name);
  expr(dfp,"%s * mat",gm->name);
  add_break(dfp);
  
  expr(dfp,"mat = (%s *) matrixpoi;",gm->name);
  add_break(dfp);
  
  for(k=0;k<gm->len;k++) {
    expr(dfp,"fprintf(ofp,\"State %9s %%d\\n\",%s_EXPL_MATRIX(mat,i,j,%s));",gm->state[k]->name,gm->name,gm->state[k]->name);
  }
  
  close_function(dfp);
  add_break(dfp);

  for(k=0;k<gm->len;k++) {
    fi = FuncInfo_named_from_varstr(FI_INTERNAL,"state_debug_%s_%s",gm->name,gm->state[k]->name);
    start_function_FuncInfo(fi,dfp,"void state_debug_%s_%s(void * matrixpoi,int i,int j,FILE * ofp)",gm->name,gm->state[k]->name);
    expr(dfp,"%s * mat",gm->name);
    add_break(dfp);

    expr(dfp,"mat = (%s *) matrixpoi;",gm->name);
    add_break(dfp);

    expr(dfp,"fprintf(ofp,\"State %s Score %%d\\n\",%s_EXPL_MATRIX(mat,i,j,%s));",gm->state[k]->name,gm->name,gm->state[k]->name);
    for(l=0;l<gm->state[k]->len;l++) {
      auto CellSource * source;
      source = gm->state[k]->source[l];
      expr(dfp,"fprintf(ofp,\"     From %8s Matrix %%5d Score %%d\\n\",%s_EXPL_%s(mat,i-%d,j-%d,%s),%s);",source->state_source,gm->name,source->isspecial == TRUE ? "SPECIAL" : "MATRIX",source->offi,source->offj,source->state_source,source->calc_expr);
    }
    close_function(dfp);
    add_break;
  }
}


# line 174 "dynadebug.dy"
void transition_debug_func(DYNFILE * dfp,GenericMatrix * gm)
{
  int k,l,m;
  FuncInfo * fi;
  ExprSet * set;
  char buffer[512];

  assert(dfp);
  assert(gm);

  for(k=0;k<gm->len;k++) {
    for(l=0;l<gm->state[k]->len;l++) {
      auto CellSource * source;
      source = gm->state[k]->source[l];
      
      fi = FuncInfo_named_from_varstr(FI_INTERNAL,"transition_debug_%s_%s_%s_%d_%d",gm->name,gm->state[k]->name,source->state_source,source->offi,source->offj);
      start_function_FuncInfo(fi,dfp,"void transition_debug_%s_%s_%s_%d_%d(void * matrixpoi,int i,int j,FILE * ofp)",gm->name,gm->state[k]->name,source->state_source,source->offi,source->offj);
      expr(dfp,"%s * mat;",gm->name);
      add_break(dfp);
      expr(dfp,"mat = (%s *) matrixpoi;",gm->name);
      add_break(dfp);

      set = build_ExprSet(source->etr);
      for(m=0;m<set->len;m++) {
	buffer[0] = '\0';
	strcat_ExprTree_Scoped(set->expr[m],buffer,gm->sc,gm->mts,NULL,NULL,NULL);
	expr(dfp,"fprintf(ofp,\"%s : %%d\\n\",%s);",buffer,buffer);
      }
      free_ExprSet(set);

      close_function(dfp);
      add_break(dfp);
    }
  }
      
}

# line 211 "dynadebug.dy"
void explicit_debug_func(DYNFILE * dfp,GenericMatrix * gm)
{
  FuncInfo * fi;
  ArgInfo  * ai;
  char * arg_str;
  char * chainstr;

  int k,l;

  /*** prepare function information ***/

  
  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"explicit_debug_%s",gm->name);
  add_line_to_Ftext(fi->ft,"This function provides debugging support",gm->name);
  add_line_to_Ftext(fi->ft,"for the dynamic programming models");
  add_break_to_Ftext(fi->ft);


  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"mat");
  ai->desc=stringallocf("%s which contains explicit basematrix memory",gm->name);


  start_function_FuncInfo(fi,dfp,"boolean explicit_debug_%s(DebugMatrix * de,boolean do_one_cell)",gm->name);
  expr(dfp,"int i;");
  expr(dfp,"int j;");
  expr(dfp,"int leni;");
  expr(dfp,"int lenj;");
  expr(dfp,"int num;");
  expr(dfp,"%s * mat;",gm->name);

  add_break(dfp);

  expr(dfp,"mat = (%s *) de->matrix;",gm->name);

  expr(dfp,"if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )");
  startbrace(dfp);
  warn_expr(dfp,"in explicit_debug_%s, passed a non Explicit matrix type, cannot calculate!",gm->name);
  expr(dfp,"return FALSE;");
  closebrace(dfp);
  add_break(dfp);
  expr(dfp,"leni = mat->leni;");
  expr(dfp,"lenj = mat->lenj;");

  

  /*** see if there any specials to specials to do ***/
  
  add_break(dfp);

  add_block_comment(dfp,"Enter debug matrix before matrix, reset");
  expr(dfp,"user_DebugMatrix(de);");
  add_break(dfp);
  expr(dfp,"if( de->reset == FALSE)");
  startbrace(dfp);
  expr(dfp,"de->currenti = de->currentj = 0;");
  closebrace(dfp);

  
  expr(dfp,"for(j=de->currentj;j<lenj;j++)");
  startbrace(dfp);
  expr(dfp,"auto int score");
  expr(dfp,"auto int temp");
  expr(dfp,"for(i=de->currenti;i<leni;i++)");
  startbrace(dfp);


  /* test for break points */

  expr(dfp,"if( do_one_cell == FALSE && should_break_DebugMatrix(de) != MDBP_NoBreakPoint ) ");
  startbrace(dfp);
  expr(dfp,"user_DebugMatrix(de)");
  add_block_comment(dfp,"May well have reset i and j here");
  expr(dfp,"if( de->reset == TRUE )");
  startbrace(dfp);
  expr(dfp,"i = de->currenti;");
  expr(dfp,"j = de->currentj;");
  expr(dfp,"continue;");
  closebrace(dfp);
  closebrace(dfp);

  
  write_score_block(dfp,gm,"EXPL_MATRIX","mat","EXPL_SPECIAL",TRUE);


  add_break(dfp);
  add_block_comment(dfp,"Now update DebugMatrix datastructure.");
  add_break(dfp);  

  for(k=0;k<gm->len;k++) {
    expr(dfp,"if( %s_EXPL_MATRIX(mat,i,j,%s) > de->max_score )",gm->name,gm->state[k]->name);
    startbrace(dfp);
    expr(dfp,"de->max_score = %s_EXPL_MATRIX(mat,i,j,%s);",gm->name,gm->state[k]->name);
    expr(dfp,"de->max_score_i = i;");
    expr(dfp,"de->max_score_j = j;");
    expr(dfp,"de->max_score_cell = %s;",gm->state[k]->name);
    closebrace(dfp);
    add_break(dfp);
  }
  expr(dfp,"de->currenti = i;");
  expr(dfp,"de->currentj = j;");

  add_break(dfp);
	
  closebrace(dfp);

  write_special_block(dfp,gm,"EXPL_MATRIX","EXPL_SPECIAL",NULL);

  add_block_comment(dfp,"Update currenti");
  expr(dfp,"de->currenti = 0;");
  add_break(dfp);


  closebrace(dfp);

  expr(dfp,"return TRUE");
  
  close_function(dfp);
  
  add_break(dfp);

}


# line 326 "dynadebug.c"
/* Function:  swap_ExprSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ExprSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ExprTree **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ExprSet(ExprTree ** list,int i,int j)  
{
    ExprTree * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ExprSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ExprSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ExprTree **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ExprSet(ExprTree ** list,int left,int right,int (*comp)(ExprTree * ,ExprTree * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ExprSet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ExprSet (list,++last,i);    
      }  
    swap_ExprSet (list,left,last);   
    qsort_ExprSet(list,left,last-1,comp);    
    qsort_ExprSet(list,last+1,right,comp);   
}    


/* Function:  sort_ExprSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ExprSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ExprSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ExprSet(ExprSet * obj,int (*comp)(ExprTree *, ExprTree *)) 
{
    qsort_ExprSet(obj->expr,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_ExprSet(obj,len)
 *
 * Descrip:    Really an internal function for add_ExprSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ExprSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ExprSet(ExprSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ExprSet called with no need");    
      return TRUE;   
      }  


    if( (obj->expr = (ExprTree ** ) ckrealloc (obj->expr,sizeof(ExprTree *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_ExprSet, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ExprSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ExprSet *]
 * Arg:        add [OWNER] Object to add to the list [ExprTree *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ExprSet(ExprSet * obj,ExprTree * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ExprSet(obj,obj->len + ExprSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->expr[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_ExprSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ExprSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ExprSet(ExprSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->expr[i] != NULL)  {  
        free_ExprTree(obj->expr[i]); 
        obj->expr[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ExprSet_alloc_std(void)
 *
 * Descrip:    Equivalent to ExprSet_alloc_len(ExprSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ExprSet *]
 *
 */
ExprSet * ExprSet_alloc_std(void) 
{
    return ExprSet_alloc_len(ExprSetLISTLENGTH); 
}    


/* Function:  ExprSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ExprSet *]
 *
 */
ExprSet * ExprSet_alloc_len(int len) 
{
    ExprSet * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ExprSet_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->expr = (ExprTree ** ) ckcalloc (len,sizeof(ExprTree *))) == NULL)   {  
      warn("Warning, ckcalloc failed in ExprSet_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ExprSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ExprSet *]
 *
 * Return [UNKN ]  Undocumented return value [ExprSet *]
 *
 */
ExprSet * hard_link_ExprSet(ExprSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ExprSet object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ExprSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ExprSet *]
 *
 */
ExprSet * ExprSet_alloc(void) 
{
    ExprSet * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ExprSet *) ckalloc (sizeof(ExprSet))) == NULL)  {  
      warn("ExprSet_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->expr = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_ExprSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ExprSet *]
 *
 * Return [UNKN ]  Undocumented return value [ExprSet *]
 *
 */
ExprSet * free_ExprSet(ExprSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ExprSet obj. Should be trappable");   
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
    if( obj->expr != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->expr[i] != NULL)    
          free_ExprTree(obj->expr[i]);   
        }  
      ckfree(obj->expr); 
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

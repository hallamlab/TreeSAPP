#ifdef _cplusplus
extern "C" {
#endif
#include "kbestsearch.h"



      
/* Function:  write_kbest_score_GenericMatrix(dfp,gm,sc,mts,dpi)
 *
 * Descrip:    Produces a kbest search type single score function
 *
 *             kbest algorithm used here is that each cell is 
 *             reduced to a single score + state number that it
 *             came from. The kbest heuristic is to provide only k paths
 *             onto the next position in the sequence. In our case we
 *             have k = length of model / number of states. This sort
 *             of kbest heurisitc is good because it cuts down on excessive
 *             book keeping of the alignments, by being able to store
 *             information of the state only
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:         sc [UNKN ] Undocumented argument [Scope *]
 * Arg:        mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
# line 63 "kbestsearch.dy"
void write_kbest_score_GenericMatrix(DYNFILE * dfp,GenericMatrix * gm,Scope * sc,MethodTypeSet * mts,DPImplementation * dpi)
{
  char subexpr_buffer[MAXLINE];
  int i,j;
  FuncInfo * fi;
  char * arg_str;
  char * chain_str;
  CommonSubExpressionSet * cses;

  cses = find_CommonSubExpressions(gm,TRUE);
  show_CommonSubExpressionSet(cses,stdout);
  

  for(j=0;j<gm->spec_len;j++)
    if( gm->special[j]->is_start == TRUE ) 
      break;

  /*** prepare function information ***/

  
  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"kbest_score_only_%s",gm->name);
  add_line_to_Ftext(fi->ft,"This function just calculates the score for the matrix",gm->name);
  add_line_to_Ftext(fi->ft,"It uses a kbest style algorithm to compute the score");
  add_line_to_Ftext(fi->ft,"It calls /allocate_%s_only",gm->name);

  arg_str = get_argstr_GenericMatrix(gm);
  add_args_GenericMatrix_FuncInfo(fi,gm);


  start_function_FuncInfo(fi,dfp,"int kbest_score_only_%s(%s)",gm->name,arg_str);

  /*** clean up ***/
  ckfree(arg_str);



  /*** into function body ***/


  expr(dfp,"int bestscore = NEGI;");
  expr(dfp,"int i;");
  expr(dfp,"int j;");
  expr(dfp,"%s * mat",gm->name);
  expr(dfp,"int * score_mat[%d];",gm->window_j+1);
  expr(dfp,"char * state_mat;",gm->window_j+1);
  expr(dfp,"int score_special[%d][%d];",gm->window_j+1,gm->spec_len);
  
  /** for the moment, ignore the possibility of static
      allocation
  */

  if( 0 && gm->qtype != NULL && gm->qtype->maxlen != NULL) {
    expr(dfp,"int internal_matrix[%d][(%s+%d) * %d];",gm->window_j+1,gm->qtype->maxlen,gm->window_i,gm->len);
    expr(dfp,"int internal_specials[%d][%d];",gm->window_j+1,gm->spec_len);
  }
  

  if(0 && dpi->largemem == TRUE ) {
    expr(dfp,"int * internal_pointer_db;");
    expr(dfp,"int * internal_special_db;");
  }

  /** kbest optimisation stuff **/

  if( dpi->dokbestcse == TRUE ) {
    for(i=0;i<cses->len;i++) {
      expr(dfp,"int subexpr%d;",i);
    }
  }

  
  add_break(dfp);
  chain_str = get_chainstr_GenericMatrix(gm);
  expr(dfp,"mat = allocate_%s_only(%s);",gm->name,chain_str);
  ckfree(chain_str);

  expr(dfp,"if( mat == NULL )");
  startbrace(dfp);
  warn_expr(dfp,"Memory allocation error in the db search - unable to communicate to calling function. this spells DISASTER!");
  expr(dfp,"return NEGI");
  closebrace(dfp);
  /*
  if(0 && gm->qtype != NULL && gm->qtype->maxlen != NULL) {
    add_block_comment(dfp,"Ok,don't need to allocate matrix as it is internal, because we have a max length");
  } else  {
    expr(dfp,"if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(%d,(mat->leni + %d)*2,%d,%d)) == NULL)",gm->window_j+1,gm->window_i,gm->window_j+1,gm->spec_len);
    startbrace(dfp);
    expr(dfp,"warn(\"Score only matrix for %s cannot be allocated, (asking for %d  by %%d  cells)\",mat->leni);",gm->name,gm->window_j);
    expr(dfp,"mat = free_%s(mat)",gm->name);
    expr(dfp,"return NEGI;");
    closebrace(dfp);
    expr(dfp,"mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;");
  }
  */

  add_block_comment(dfp,"Allocate memory for matrix");
  expr(dfp,"for(j=0;j<%d;j++)",gm->window_j+1);
  startbrace(dfp);
  expr(dfp,"score_mat[j] = (int *) ckalloc (sizeof(int) * (mat->leni+%d));",gm->window_i);
  closebrace(dfp);
  expr(dfp,"state_mat = (char *) ckalloc (sizeof(char) * (mat->leni+%d)*%d *%d);",gm->window_i,gm->len,gm->window_j+1);
  add_break(dfp);

  add_block_comment(dfp,"Now, initiate matrix");

  expr(dfp,"for(j=0;j<%d;j++)",gm->window_j+2);
  startbrace(dfp);
  expr(dfp,"for(i=(-%d);i<mat->leni;i++)",gm->window_i);
  startbrace(dfp);
  expr(dfp,"%s_KBEST_MATRIX_SCORE(mat,i,j,KBEST_SCORE) = NEGI;",gm->name);
  for(i=0;i<gm->spec_len;i++) {
    expr(dfp,"%s_KBEST_SPECIAL(mat,0,j,%s) = %s;",gm->name,gm->special[i]->name,gm->special[i]->def_score);
  }
  closebrace(dfp);
  closebrace(dfp);

  add_break(dfp);

  add_block_comment(dfp,"Ok, lets do-o-o-o-o it");

  add_break(dfp);

  if( dpi->dokbestcse == TRUE ) {
    for(i=0;i<cses->len;i++) {
      if( IS_NON_IJ_DEP_CSE(cses->cse[i]) ) {  
	subexpr_buffer[0]='\0';
	strcat_ExprTree_Scoped(cses->cse[i]->expr,subexpr_buffer,sc,mts,dpi->dycw,NULL,NULL);
	expr(dfp,"subexpr%d = %s;",i,subexpr_buffer);
      }
    }
  }

  expr(dfp,"for(j=0;j<mat->lenj;j++)");
  startbrace_tag(dfp,"for all target positions");
  expr(dfp,"auto int score");
  expr(dfp,"auto int temp");
  expr(dfp,"auto int state");
  expr(dfp,"auto int temp_state");

  if( dpi->dokbestcse == TRUE ) {
    for(i=0;i<cses->len;i++) {
      if( IS_J_DEP_CSE(cses->cse[i]) == TRUE && IS_I_DEP_CSE(cses->cse[i]) == FALSE ) {  
	subexpr_buffer[0]='\0';
	strcat_ExprTree_Scoped(cses->cse[i]->expr,subexpr_buffer,sc,mts,dpi->dycw,NULL,NULL);
	expr(dfp,"subexpr%d = %s;",i,subexpr_buffer);
      }
    }
  }

  add_break(dfp);
  add_block_comment(dfp,"Initialise these specials");
  for(i=0;i<gm->spec_len;i++) {
    expr(dfp,"%s_KBEST_SPECIAL(mat,0,j,%s) = %s;",gm->name,gm->special[i]->name,gm->special[i]->def_score);
  }
  expr(dfp,"for(i=0;i<mat->leni;i++)");
  startbrace_tag(dfp,"for all query positions");


  /** kbest cse optimisations **/

  if( dpi->dokbestcse == TRUE ) {
    for(i=0;i<cses->len;i++) {
      if( IS_J_DEP_CSE(cses->cse[i]) == TRUE && IS_I_DEP_CSE(cses->cse[i]) == TRUE ) {  
	subexpr_buffer[0]='\0';
	strcat_ExprTree_Scoped(cses->cse[i]->expr,subexpr_buffer,sc,mts,dpi->dycw,NULL,NULL);
	expr(dfp,"subexpr%d = %s;",i,subexpr_buffer);
      }
    }
  }

  write_kbest_block(dfp,gm,"KBEST_MATRIX","mat","KBEST_SPECIAL",TRUE,cses,mts,dpi);
  
  closebrace(dfp);

  add_break(dfp);

  write_special_block(dfp,gm,"KBEST_MATRIX","KBEST_SPECIAL","bestscore");

  closebrace(dfp);

  add_break(dfp);

  expr(dfp,"mat = free_%s(mat)",gm->name);

  expr(dfp,"return bestscore;");

  close_function(dfp);

  add_break(dfp);
}

      
/* Function:  write_kbest_block(dfp,gm,matrixtag,pointertag,specialtag,use_special,cses,mts,dpi)
 *
 * Descrip:    Produces the actual kbest scoring inner loop
 *
 *
 * Arg:                dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:                 gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:          matrixtag [UNKN ] Undocumented argument [char *]
 * Arg:         pointertag [UNKN ] Undocumented argument [char *]
 * Arg:         specialtag [UNKN ] Undocumented argument [char *]
 * Arg:        use_special [UNKN ] Undocumented argument [boolean]
 * Arg:               cses [UNKN ] Undocumented argument [CommonSubExpressionSet *]
 * Arg:                mts [UNKN ] Undocumented argument [MethodTypeSet *]
 * Arg:                dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
# line 258 "kbestsearch.dy"
void write_kbest_block(DYNFILE * dfp,GenericMatrix * gm,char * matrixtag,char * pointertag,char * specialtag,boolean use_special,CommonSubExpressionSet * cses,MethodTypeSet * mts,DPImplementation * dpi)
{
  char buffer[MAXLINE];
  register int i;
  register int j;
  register int k;
  TransitionSet * ts;

  int current_i_off;
  int current_j_off;

  ts = TransitionSet_from_GenericMatrix(gm);
  sort_TransitionSet_offset(ts);
  /*  show_TransitionSet(ts,stderr); */
  add_break(dfp);
  expr(dfp,"score = NEGI;");
  expr(dfp,"state = 300;");
  add_break(dfp);

  /* do switch statements from specials first */

  for(k=0;k<ts->len && Transition_from_special(ts->trans[k]) == 1;k++) {
    if( ts->trans[k]->trans_type == TRANSITION_FROM_START ) {
      add_block_comment(dfp,"Source is the start state, and hence no need to look, must be 0");
      expr(dfp,"temp = %s",
	   ts->trans[k]->calc
	   );
    } else {
      expr(dfp,"temp = %s_KBEST_SPECIAL(mat,i-%d,j-%d,%s) + %s",gm->name,
	   ts->trans[k]->offi,
	   ts->trans[k]->offj,
	   ts->trans[k]->from,
	   ts->trans[k]->calc
	   );
    }

    expr(dfp,"if(temp > score)");
    startbrace(dfp);
    expr(dfp,"score = temp;");
    expr(dfp,"state = %s;",ts->trans[k]->to);
    closebrace(dfp);
  }

  /* do switch statements over the possibilities for each i,j offset */

  for(;k<ts->len;) {
    current_i_off = ts->trans[k]->offi;
    current_j_off = ts->trans[k]->offj;
    expr(dfp,"switch ( %s_KBEST_MATRIX_STATE(mat,i-%d,j-%d,DUMMY_KBEST_STATE) )",gm->name,current_i_off,current_j_off);
    startbrace(dfp);
    for(;k < ts->len && ts->trans[k]->offi == current_i_off && ts->trans[k]->offj == current_j_off;k++) {
      expr(dfp,"case %s :",ts->trans[k]->from);
      startcase(dfp);
      if( dpi->dokbestcse == TRUE ) {
	buffer[0] = '\0';
	strcat(buffer,"(");
	strcat_cses_ExprTree(ts->trans[k]->expr,buffer,gm->sc,mts,dpi);
	strcat(buffer,")");
	if( ts->trans[k]->expr_state != NULL ) {
	  strcat(buffer,"+ (");
	  strcat_cses_ExprTree(ts->trans[k]->expr_state,buffer,gm->sc,mts,dpi);
	  strcat(buffer,")");
	}
	expr(dfp,"temp = %s;",buffer);
      } else {
	expr(dfp,"temp = %s;",ts->trans[k]->calc);
      }
      expr(dfp,"temp_state = %s;",ts->trans[k]->to);
      expr(dfp,"break;");
      closecase(dfp);
    }
    expr(dfp,"default :");
    startcase(dfp);
    expr(dfp,"temp = NEGI;");
    expr(dfp,"temp_state = (200);");
    expr(dfp,"break;");
    closecase(dfp);
    closebrace(dfp); /* close switch statement */
    add_break(dfp);

    expr(dfp,"if( temp > NEGI )");
    startbrace(dfp);
    expr(dfp,"temp += %s_KBEST_MATRIX_SCORE(mat,i-%d,j-%d,KBEST_SCORE);",gm->name,current_i_off,current_j_off);
    expr(dfp,"if( temp > score )");
    startbrace(dfp);
    expr(dfp,"score = temp;");
    expr(dfp,"state = temp_state;");
    closebrace(dfp);
    closebrace(dfp);
    add_break(dfp);
  }
      
  expr(dfp,"%s_KBEST_MATRIX_SCORE(mat,i,j,KBEST_SCORE) = score;",gm->name);
  expr(dfp,"%s_KBEST_MATRIX_STATE(mat,i,j,KBEST_STATE) = state;",gm->name);

  add_break(dfp);
  add_block_comment(dfp,"Now do any potential main to special movements");
  
  for(i=0;i<gm->spec_len;i++) {
    for(j=0;j<gm->special[i]->len;j++) {
      if( gm->special[i]->source[j]->isspecial == FALSE ) {
	expr(dfp,"if( state == %s)",gm->special[i]->source[j]->state_source);
	startbrace(dfp);
	expr(dfp,"temp = score + (%s) + (%s);",
	     gm->special[i]->source[j]->calc_expr,
	     gm->special[i]->calc_expr == NULL ? "0" : gm->special[i]->calc_expr);
	expr(dfp,"if( temp > %s_KBEST_SPECIAL(mat,i,j,%s) )",gm->name,gm->special[i]->name);
	hang_expr(dfp,"%s_KBEST_SPECIAL(mat,i,j,%s) = temp;",gm->name,gm->special[i]->name);
	closebrace(dfp);
      }
    }
  }

}


/* Function:  TransitionSet_from_GenericMatrix(gm)
 *
 * Descrip:    Makes a transition set from a generic matrix
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
# line 377 "kbestsearch.dy"
TransitionSet * TransitionSet_from_GenericMatrix(GenericMatrix * gm)
{
  TransitionSet *out;
  int i,j,k;
  Transition * trans;
  char buffer[1024];
  
  out = TransitionSet_alloc_std();

  for(i=0;i<gm->len;i++)
    for(j=0;j<gm->state[i]->len;j++) {
      trans = Transition_alloc();
      trans->offi = gm->state[i]->source[j]->offi;
      trans->offj = gm->state[i]->source[j]->offj;
      trans->to   = stringalloc(gm->state[i]->name);
      trans->from = stringalloc(gm->state[i]->source[j]->state_source);
      trans->expr = gm->state[i]->source[j]->etr;
      trans->expr_state = gm->state[i]->etr;
      /* if it is from a special, see if it is from start */
      if( gm->state[i]->source[j]->isspecial ) {

	/* this is painful, but because I hadn't got the previous
	   data structure right! */

	for(k=0;k<gm->spec_len;k++) {
	  if( gm->special[k]->is_start == TRUE && strcmp(gm->state[i]->source[j]->state_source,gm->special[k]->name) == 0) {
	    trans->trans_type = TRANSITION_FROM_START;
	    break;
	  }
	}
	if( k == gm->spec_len ) {
	  trans->trans_type = TRANSITION_FROM_SPECIAL;
	}
      }

      if( gm->state[i]->calc_expr != NULL ) {
	sprintf(buffer,"(%s) + (%s)",gm->state[i]->source[j]->calc_expr,gm->state[i]->calc_expr);
      } else {
	sprintf(buffer,"%s",gm->state[i]->source[j]->calc_expr);
      }
      trans->calc = stringalloc(buffer);
      add_TransitionSet(out,trans);
    }

  return out;
}

/* Function:  can_kbest_GenericMatrix(gm)
 *
 * Descrip:    sees whether the generic matrix is suitable for kbest optimisations
 *
 *
 * Arg:        gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 427 "kbestsearch.dy"
boolean can_kbest_GenericMatrix(GenericMatrix * gm)
{
  int i,j;
  int k,l;

  /* have to find cases in which the source,i,j is 
     the same twice */

  for(i=0;i<gm->len;i++) {
    auto CellState * state;
    state = gm->state[i];
    
    for(j=0;j<state->len;j++) {
      /* now loop over all other transitions in this state */
      for(l=j+1;l<state->len;l++)
	if( state->source[j]->offi == state->source[l]->offi &&
	    state->source[j]->offj == state->source[l]->offj &&
	    strcmp(state->source[j]->state_source,state->source[l]->state_source) == 0 )
	  return FALSE;
      

      /* now loop over all other states */

      for(k=i+1;k<gm->len;k++) {
	for(l=0;l<gm->state[k]->len;l++)
	  if( state->source[j]->offi == gm->state[k]->source[l]->offi &&
	      state->source[j]->offj == gm->state[k]->source[l]->offj &&
	      strcmp(state->source[j]->state_source,gm->state[k]->source[l]->state_source) == 0 )
	    return FALSE;
      }/* all other states */
    } /* all transitions in this particular i'th state */
  }/* all states */

  return TRUE;
}


/* Function:  sort_TransitionSet_offset(ts)
 *
 * Descrip:    Sorts by offi then offj
 *
 *
 * Arg:        ts [UNKN ] Undocumented argument [TransitionSet *]
 *
 */
# line 467 "kbestsearch.dy"
void sort_TransitionSet_offset(TransitionSet * ts)
{
  sort_TransitionSet(ts,comp_Transition);
}

/* Function:  comp_Transition(two,one)
 *
 * Descrip:    comparison by offi/offj
 *
 *
 * Arg:        two [UNKN ] Undocumented argument [Transition *]
 * Arg:        one [UNKN ] Undocumented argument [Transition *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 475 "kbestsearch.dy"
int comp_Transition(Transition * two,Transition * one)
{
  if( Transition_from_special(one) != Transition_from_special(two) ) {
    if( Transition_from_special(two) ) {
      return -1;
    } else {
      return 1;
    }
  }

  if( one->offi == two->offi ) {
    return one->offj - two->offj;
  }

  return one->offi - two->offi;
}

	
/* Function:  show_TransitionSet(tset,ofp)
 *
 * Descrip:    Shows a transition set
 *
 *
 * Arg:        tset [UNKN ] Undocumented argument [TransitionSet *]
 * Arg:         ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 496 "kbestsearch.dy"
void show_TransitionSet(TransitionSet * tset,FILE * ofp)
{
  int i;

  for(i=0;i<tset->len;i++)
    show_Transition(tset->trans[i],ofp);
}

/* Function:  show_Transition(trans,ofp)
 *
 * Descrip:    Shows a transition
 *
 *
 * Arg:        trans [UNKN ] Undocumented argument [Transition *]
 * Arg:          ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 508 "kbestsearch.dy"
void show_Transition(Transition * trans,FILE * ofp)
{
  fprintf(ofp,"Transition [%s %d %d %s]\n",trans->from,trans->offi,trans->offj,trans->to);
  fprintf(ofp,"   calc %s\n",trans->calc);
}




# line 526 "kbestsearch.c"
/* Function:  hard_link_Transition(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Transition *]
 *
 * Return [UNKN ]  Undocumented return value [Transition *]
 *
 */
Transition * hard_link_Transition(Transition * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Transition object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Transition_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Transition *]
 *
 */
Transition * Transition_alloc(void) 
{
    Transition * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Transition *) ckalloc (sizeof(Transition))) == NULL)    {  
      warn("Transition_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->to = NULL;  
    out->from = NULL;    
    out->offi = 0;   
    out->offj = 0;   
    out->calc = NULL;    
    out->trans_type = TRANSITION_NORMAL; 
    out->expr = NULL;    
    out->expr_state = NULL;  


    return out;  
}    


/* Function:  free_Transition(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Transition *]
 *
 * Return [UNKN ]  Undocumented return value [Transition *]
 *
 */
Transition * free_Transition(Transition * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Transition obj. Should be trappable");    
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
    if( obj->to != NULL) 
      ckfree(obj->to);   
    if( obj->from != NULL)   
      ckfree(obj->from);     
    if( obj->calc != NULL)   
      ckfree(obj->calc);     
    if( obj->expr != NULL)   
      free_ExprTree(obj->expr);  
    if( obj->expr_state != NULL) 
      free_ExprTree(obj->expr_state);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_TransitionSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_TransitionSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Transition **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_TransitionSet(Transition ** list,int i,int j)  
{
    Transition * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_TransitionSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_TransitionSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Transition **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_TransitionSet(Transition ** list,int left,int right,int (*comp)(Transition * ,Transition * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_TransitionSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_TransitionSet (list,++last,i);  
      }  
    swap_TransitionSet (list,left,last); 
    qsort_TransitionSet(list,left,last-1,comp);  
    qsort_TransitionSet(list,last+1,right,comp); 
}    


/* Function:  sort_TransitionSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_TransitionSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [TransitionSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_TransitionSet(TransitionSet * obj,int (*comp)(Transition *, Transition *)) 
{
    qsort_TransitionSet(obj->trans,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_TransitionSet(obj,len)
 *
 * Descrip:    Really an internal function for add_TransitionSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransitionSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_TransitionSet(TransitionSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_TransitionSet called with no need");  
      return TRUE;   
      }  


    if( (obj->trans = (Transition ** ) ckrealloc (obj->trans,sizeof(Transition *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_TransitionSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_TransitionSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransitionSet *]
 * Arg:        add [OWNER] Object to add to the list [Transition *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_TransitionSet(TransitionSet * obj,Transition * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_TransitionSet(obj,obj->len + TransitionSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->trans[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_TransitionSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransitionSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_TransitionSet(TransitionSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->trans[i] != NULL) {  
        free_Transition(obj->trans[i]);  
        obj->trans[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  TransitionSet_alloc_std(void)
 *
 * Descrip:    Equivalent to TransitionSet_alloc_len(TransitionSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
TransitionSet * TransitionSet_alloc_std(void) 
{
    return TransitionSet_alloc_len(TransitionSetLISTLENGTH); 
}    


/* Function:  TransitionSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
TransitionSet * TransitionSet_alloc_len(int len) 
{
    TransitionSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = TransitionSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->trans = (Transition ** ) ckcalloc (len,sizeof(Transition *))) == NULL)  {  
      warn("Warning, ckcalloc failed in TransitionSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_TransitionSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransitionSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
TransitionSet * hard_link_TransitionSet(TransitionSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransitionSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransitionSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
TransitionSet * TransitionSet_alloc(void) 
{
    TransitionSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransitionSet *) ckalloc (sizeof(TransitionSet))) == NULL)  {  
      warn("TransitionSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->trans = NULL;   
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_TransitionSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransitionSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransitionSet *]
 *
 */
TransitionSet * free_TransitionSet(TransitionSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransitionSet obj. Should be trappable"); 
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
    if( obj->trans != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->trans[i] != NULL)   
          free_Transition(obj->trans[i]);    
        }  
      ckfree(obj->trans);    
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

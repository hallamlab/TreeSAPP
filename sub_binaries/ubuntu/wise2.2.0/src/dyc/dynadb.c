#ifdef _cplusplus
extern "C" {
#endif
#include "dynadb.h"

/* Function:  make_search_loop_function(dfp,gm)
 *
 * Descrip:    Makes the serial search function, which loops
 *             over databases 
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
# line 20 "dynadb.dy"
void make_search_loop_function(DYNFILE * dfp,GenericMatrix * gm)
{
  int i;
  char buffer[MAXLINE];
  FuncInfo * fi;
  boolean qdb = FALSE;
  boolean tdb = FALSE;
  char * cstr;

  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"serial_search_%s",gm->name);
  add_line_to_Ftext(fi->ft,"This function makes a database search of %s",gm->name);
  add_line_to_Ftext(fi->ft,"It is a single processor implementation");

  if( gm->qtype != NULL && gm->qtype->is_database == TRUE) 
    qdb = TRUE;
  if( gm->ttype != NULL && gm->ttype->is_database == TRUE) 
    tdb = TRUE;

  if( qdb == TRUE)
    sprintf(buffer,"Search_Return_Type serial_search_%s(Hscore * out,%s querydb,",gm->name,gm->qtype->database_type);
  else sprintf(buffer,"Search_Return_Type serial_search_%s(Hscore * out,%s %s,",gm->name,gm->query->element_type,gm->query->name);
  
  if( tdb == TRUE) {
    strcat(buffer,gm->ttype->database_type);
    strcat(buffer," targetdb ");
  } else {
    strcat(buffer,gm->target->element_type);
    strcat(buffer," ");
    strcat(buffer,gm->target->name);
    strcat(buffer," ");
  }
  
  for(i=0;i<gm->res_len;i++) {
    strcat(buffer,",");
    strcat(buffer,gm->resource[i]->element_type);
    strcat(buffer," ");
    strcat(buffer,gm->resource[i]->name);
  }
  strcat(buffer,")");


  start_function_FuncInfo(fi,dfp,buffer);

  if( qdb == TRUE )
    expr(dfp,"%s %s",gm->query->element_type,gm->query->name);

  if( tdb == TRUE )
    expr(dfp,"%s %s",gm->target->element_type,gm->target->name);
  expr(dfp,"int db_status");
  expr(dfp,"int score");
  expr(dfp,"int query_pos = 0;");
  expr(dfp,"int target_pos = 0;");
  expr(dfp,"DataScore * ds;");

  add_break(dfp);
  expr(dfp,"push_errormsg_stack(\"Before any actual search in db searching\");");

  if( qdb == TRUE) {
    expr(dfp,"%s = %s(querydb,&db_status);",gm->query->name,gm->qtype->init_func);
    expr(dfp,"if( db_status == DB_RETURN_ERROR ) ");
    startbrace(dfp);
    warn_expr(dfp,"In searching %s, got a database reload error on the query [%s] database",gm->name,gm->query->name);
    expr(dfp,"return SEARCH_ERROR;");
    closebrace(dfp);
    expr(dfp,"for(;;)");
    startbrace_tag(dfp,"For all query entries");
  }
  add_break(dfp);

  expr(dfp,"target_pos = 0");


  add_break(dfp);

  if( tdb == TRUE) {
    expr(dfp,"%s = %s(targetdb,&db_status);",gm->target->name,gm->ttype->init_func);
    
    expr(dfp,"if( db_status == DB_RETURN_ERROR ) ");
    startbrace(dfp);
    warn_expr(dfp,"In searching %s, got a database init error on the target [%s] database",gm->name,gm->target->name);
    expr(dfp,"return SEARCH_ERROR;");
    closebrace(dfp);
    expr(dfp,"for(;;)");
    startbrace_tag(dfp,"For all target entries");
  }

  add_break(dfp);

  cstr = get_chainstr_GenericMatrix(gm);

  if(gm->qtype != NULL && gm->qtype->maxlen != NULL) {
    add_block_comment(dfp,"Maximum length to search - should check");
    expr(dfp,"if( %s->%s > %s ) ",gm->query->name,gm->query_len,gm->qtype->maxlen);
    startbrace_tag(dfp,"if over length");
    warn_expr(dfp,"A query over the length when maxlen provided. Problem!");
    expr(dfp,"score = -10000;\n");
    closebrace(dfp);
    expr(dfp,"else");
    hang_expr(dfp,"score = score_only_%s(%s)",gm->name,cstr);
  } else {
    add_block_comment(dfp,"No maximum length - allocated on-the-fly");
    expr(dfp,"score = score_only_%s(%s)",gm->name,cstr);
  }

  ckfree(cstr);
  expr(dfp,"if( should_store_Hscore(out,score) == TRUE ) ");
  startbrace_tag(dfp,"if storing datascore");
  expr(dfp,"ds = new_DataScore_from_storage(out)");
  expr(dfp,"if( ds == NULL ) ");
  startbrace(dfp);
  warn_expr(dfp,"%s search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure",gm->name);
  expr(dfp,"return SEARCH_ERROR;");
  closebrace(dfp);
  add_block_comment(dfp,"Now: add query/target information to the entry");

  if( qdb == TRUE)
    expr(dfp,"%s(ds->query,%s,querydb)",gm->qtype->dataentry_add,gm->query->name);
  if( tdb == TRUE)
    expr(dfp,"%s(ds->target,%s,targetdb)",gm->ttype->dataentry_add,gm->target->name);

  expr(dfp,"ds->score = score");
  expr(dfp,"add_Hscore(out,ds);");
  closebrace(dfp); /* end of if stores */

  expr(dfp,"pop_errormsg_stack()");
  expr(dfp,"push_errormsg_stack(\"DB searching: just finished [Query Pos: %%d] [Target Pos: %%d]\",query_pos,target_pos);");
  add_break(dfp);

  if( tdb == TRUE ) {
    expr(dfp," %s = %s(%s,targetdb,&db_status)",gm->target->name,gm->ttype->reload_func,gm->target->name);
    expr(dfp,"if( db_status == DB_RETURN_ERROR )");
    startbrace(dfp);
    expr(dfp,"warn(\"In searching %s, Reload error on database %s, position %%d,%%d\",query_pos,target_pos);",gm->name,gm->target->name); 
    expr(dfp,"return SEARCH_ERROR;");
    closebrace(dfp);
    expr(dfp,"if( db_status == DB_RETURN_END )");
    hang_expr(dfp,"break;");
    add_end_comment(dfp,"Out of target loop");
    expr(dfp,"target_pos++;");
    closebrace(dfp);
    expr(dfp,"%s(%s,targetdb)",gm->ttype->close_func,gm->target->name);
  }

  if( qdb == TRUE ) {
    expr(dfp," %s = %s(%s,querydb,&db_status)",gm->query->name,gm->qtype->reload_func,gm->query->name);
    expr(dfp,"if( db_status == DB_RETURN_ERROR)");
    startbrace(dfp);
    expr(dfp,"warn(\"In searching %s, Reload error on database %s, position %%d,%%d\",query_pos,target_pos);",gm->name,gm->query->name); 
    expr(dfp,"return SEARCH_ERROR;");
    closebrace(dfp);
    expr(dfp,"if( db_status == DB_RETURN_END)");
    hang_expr(dfp,"break;");
    add_end_comment(dfp,"Out of query loop");
    expr(dfp,"query_pos++;");

    closebrace(dfp);
    expr(dfp,"%s(%s,querydb)",gm->qtype->close_func,gm->query->name);
  }
  expr(dfp,"pop_errormsg_stack()");
  expr(dfp,"return SEARCH_OK;");

  close_function(dfp);
  add_break(dfp);
}


/* Function:  write_one_score_GenericMatrix(dfp,gm,dpi)
 *
 * Descrip:    Makes the score only function, which gives the score
 *             for two objects. Used in the serial and the pthreads
 *             ports
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
# line 191 "dynadb.dy"
void write_one_score_GenericMatrix(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi)
{
  int i,j;
  FuncInfo * fi;
  char * arg_str;
  char * chain_str;

  for(j=0;j<gm->spec_len;j++)
    if( gm->special[j]->is_start == TRUE ) 
      break;

  /*** prepare function information ***/

  
  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"score_only_%s",gm->name);
  add_line_to_Ftext(fi->ft,"This function just calculates the score for the matrix",gm->name);
  add_line_to_Ftext(fi->ft,"I am pretty sure we can do this better, but hey, for the moment...");
  add_line_to_Ftext(fi->ft,"It calls /allocate_%s_only",gm->name);

  arg_str = get_argstr_GenericMatrix(gm);
  add_args_GenericMatrix_FuncInfo(fi,gm);


  start_function_FuncInfo(fi,dfp,"int score_only_%s(%s)",gm->name,arg_str);

  /*** clean up ***/
  ckfree(arg_str);



  /*** into function body ***/


  expr(dfp,"int bestscore = NEGI;");
  expr(dfp,"int i;");
  expr(dfp,"int j;");
  expr(dfp,"int k;");
  expr(dfp,"%s * mat",gm->name);

  if(gm->qtype != NULL && gm->qtype->maxlen != NULL) {
    expr(dfp,"int internal_matrix[%d][(%s+%d) * %d];",gm->window_j+1,gm->qtype->maxlen,gm->window_i,gm->len);
    expr(dfp,"int internal_specials[%d][%d];",gm->window_j+1,gm->spec_len);
  }

  if( dpi->largemem == TRUE ) {
    expr(dfp,"int * internal_pointer_db;");
    expr(dfp,"int * internal_special_db;");
  }
  
  add_break(dfp);
  chain_str = get_chainstr_GenericMatrix(gm);
  expr(dfp,"mat = allocate_%s_only(%s);",gm->name,chain_str);
  ckfree(chain_str);

  expr(dfp,"if( mat == NULL )");
  startbrace(dfp);
  warn_expr(dfp,"Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");
  expr(dfp,"return NEGI");
  closebrace(dfp);

  if(gm->qtype != NULL && gm->qtype->maxlen != NULL) {
    add_block_comment(dfp,"Ok,don't need to allocate matrix as it is internal, because we have a max length");
  } else if ( dpi->largemem == TRUE ) {
    expr(dfp,"if( (internal_pointer_db = (int *)ckcalloc(((mat->leni+%d) * %d),sizeof(int))) == NULL)",
	 gm->window_i,((gm->window_j+1)*gm->len));
    hang_expr(dfp,"fatal(\"could not allocate internal matrix in long mode\");");
    expr(dfp,"if( (internal_special_db = (int *)ckcalloc(%d,sizeof(int))) == NULL)",
	 ((gm->window_j+1)*gm->spec_len));
    hang_expr(dfp,"fatal(\"could not allocate internal matrix in long mode\");"); 
  } else {
    expr(dfp,"if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(%d,(mat->leni + %d) * %d,%d,%d)) == NULL)",gm->window_j+1,gm->window_i,gm->len,gm->window_j+1,gm->spec_len);
    startbrace(dfp);
    expr(dfp,"warn(\"Score only matrix for %s cannot be allocated, (asking for %d  by %%d  cells)\",mat->leni*%d);",gm->name,gm->window_j,gm->len);
    expr(dfp,"mat = free_%s(mat)",gm->name);
    expr(dfp,"return 0;");
    closebrace(dfp);
    expr(dfp,"mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;");
  }

  add_break(dfp);

  add_block_comment(dfp,"Now, initiate matrix");

  expr(dfp,"for(j=0;j<%d;j++)",gm->window_j+2);
  startbrace(dfp);
  expr(dfp,"for(i=(-%d);i<mat->leni;i++)",gm->window_i);
  startbrace(dfp);
  expr(dfp,"for(k=0;k<%d;k++)",gm->len);
  hang_expr(dfp,"%s_VSMALL_MATRIX(mat,i,j,k) = NEGI;",gm->name);
  closebrace(dfp);
  for(i=0;i<gm->spec_len;i++) {
    expr(dfp,"%s_VSMALL_SPECIAL(mat,i,j,%s) = %s;",gm->name,gm->special[i]->name,gm->special[i]->def_score);
  }
  closebrace(dfp);

  add_break(dfp);

  add_block_comment(dfp,"Ok, lets do-o-o-o-o it");

  add_break(dfp);

  expr(dfp,"for(j=0;j<mat->lenj;j++)");
  startbrace_tag(dfp,"for all target positions");
  expr(dfp,"auto int score");
  expr(dfp,"auto int temp");
  if( dpi->largemem == TRUE ) {
    add_block_comment(dfp,"Need to reset START to 0");
    expr(dfp,"%s_VSMALL_SPECIAL(mat,0,j,%s) = 0",gm->name,gm->special[j]->name);
  }
  expr(dfp,"for(i=0;i<mat->leni;i++)");
  startbrace_tag(dfp,"for all query positions");
  

  write_score_block(dfp,gm,"VSMALL_MATRIX","mat","VSMALL_SPECIAL",TRUE);
  
  closebrace(dfp);

  add_break(dfp);

  write_special_block(dfp,gm,"VSMALL_MATRIX","VSMALL_SPECIAL","bestscore");

  closebrace(dfp);

  add_break(dfp);

  expr(dfp,"mat = free_%s(mat)",gm->name);

  expr(dfp,"return bestscore;");

  close_function(dfp);

  add_break(dfp);
}
  

  

 

  

  

  




  

# line 341 "dynadb.c"

#ifdef _cplusplus
}
#endif

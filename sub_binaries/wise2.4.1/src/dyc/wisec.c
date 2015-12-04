#ifdef _cplusplus
extern "C" {
#endif
#include "wisec.h"



# line 45 "wisec.dy"
char * name_from_func_type(char * ty)
{
  char buffer[MAXLINE];
  char * runner;
  char * type = ty;



  for(;*type && *type != '(';type++)
    ;
  for(;*type && *type != '*';type++) 
    ;
  type++;
  for(;*type && isspace(*type);type++)
    ;

  if( *type == '\0') {
    warn("Cannot make %s into a function name. not good",ty);
    return NULL;
  }

  for(runner=buffer;!isspace(*type) && *type && *type != ')' ;type++,runner++)
    *runner = (*type);

  *runner = '\0';

  return stringalloc(buffer);
}


# line 75 "wisec.dy"
void write_StructHolder_header  (DYNFILE * dfp,StructHolder * sh)
{
  if( sh->oi != NULL ) {
    fprintf(dfp->head,"/* Object %s\n *\n",sh->name);
    write_C_ObjectInfo(sh->oi,dfp->head);
    fprintf(dfp->head," *\n */\n");
  }

  write_StructHolder_typedef(dfp,sh,NULL);
}

# line 86 "wisec.dy"
void write_StructHolder_function(DYNFILE * dfp,StructHolder * sh,ModuleFunctionList * mfl)
{
  register int j;
  ModuleFunction * mf;
  boolean  is_hard_link = FALSE;	

	
  is_hard_link = dfp->should_hard_link;

  if( sh == NULL || sh->name == NULL ) {
    warn("Very bad problem. You are trying to build a struct holder function without a valid structholder object or holder");
  }


  mf = get_ModuleFunction_from_name(mfl,sh->name);


  if( is_listfunc(sh) == TRUE)
    {
      for(j=0;j<sh->len;j++)
	if( sh->el[j]->islist == TRUE)
	  {
	    write_list_sort_function(dfp,sh,sh->el[j]);
	    write_list_expand_function(dfp,sh,sh->el[j]);
	    write_list_add_function(dfp,sh,sh->el[j]);
	    write_list_flush_function(dfp,sh,sh->el[j]);
	  }
      write_alloc_std_function(dfp,sh);
      write_alloc_len_function(dfp,sh);
    }
  else if( is_matrixfunc(sh) == TRUE)
    {
      write_alloc_matrix_function(dfp,sh);
      write_expand_matrix_function(dfp,sh);
    }

  /*** currently disabled
    write_binary_dump_function(sh);
    write_binary_read_function(sh);
    ****/

  if( is_hard_link == TRUE)
    write_hard_link_function(dfp,sh);

  if( mf == NULL || mf->has_cons == FALSE)
    write_simplealloc_function(dfp,sh);
  
  if( mf == NULL || mf->has_decons == FALSE)
    write_free_function(dfp,sh);


}


# line 140 "wisec.dy"
void write_StructHolder_typedef(DYNFILE * dfp,StructHolder * sh,FILE * ofp)
{
  register int i;
  boolean  is_hard_link = FALSE;

  is_hard_link = dfp->should_hard_link;

  start_struct(dfp,"struct %s%s",dfp->package_name == NULL ? "" : dfp->package_name,sh->name);

  if( is_hard_link == TRUE ) {
    struct_expr(dfp,"int dynamite_hard_link;");
    struct_macro(dfp,"#ifdef PTHREAD");
    struct_expr(dfp,"pthread_mutex_t dynamite_mutex");
    struct_macro(dfp,"#endif");
  }

  for(i=0;i<sh->len;i++)
    write_StructElement_typedef(dfp,sh->el[i]);
	
  close_struct(dfp,"; ",sh->name);
  add_end_comment_header(dfp,"%s defined",sh->name);

  fprintf(dfp->head,"#ifndef DYNAMITE_DEFINED_%s\n",sh->name);
  if( dfp->package_name != NULL ) {
    fprintf(dfp->head,"typedef struct %s%s %s%s;\n",dfp->package_name,sh->name,dfp->package_name,sh->name);
    fprintf(dfp->head,"#define %s %s%s\n",sh->name,dfp->package_name,sh->name);
  } else {
    fprintf(dfp->head,"typedef struct %s %s;\n",sh->name,sh->name);
  }

  fprintf(dfp->head,"#define DYNAMITE_DEFINED_%s\n",sh->name);
  fprintf(dfp->head,"#endif\n");


  fprintf(dfp->head,"\n\n");

}

# line 178 "wisec.dy"
void write_StructElement_typedef(DYNFILE * dfp,StructElement * se)
{
  if( se->isfunc == TRUE )
    {
      struct_expr(dfp,"%s;",se->element_type);
      return;
    }

  struct_expr(dfp,"%s %s; ",se->element_type,se->name);
  striptoprint(se->comment);
  if( se->comment != NULL)
    add_end_comment_header(dfp,se->comment);

  if( se->islist == TRUE)
    {
      struct_expr(dfp,"int %slen;",CKN(se->len_append));
      add_end_comment_header(dfp,"len for above %s ",se->name);
      struct_expr(dfp,"int %smaxlen;",CKN(se->len_append));
      add_end_comment_header(dfp,"maxlen for above %s",se->name);
    }
  else if( se->ismatrix == TRUE)
    {
      struct_expr(dfp,"int leni");
      add_end_comment_header(dfp,"leni for above %s ",se->name);
      struct_expr(dfp,"int maxleni;");
      add_end_comment_header(dfp,"max length for above pointer set");
		
      struct_expr(dfp,"int lenj");
      add_end_comment_header(dfp,"lenj for above %s ",se->name);
      struct_expr(dfp,"int maxlenj;");
      add_end_comment_header(dfp,"max length for above pointer set");
    }
	
}

# line 213 "wisec.dy"
char * free_function(char * type)
{
  char * runner;
  char * temp;
  char buffer[MAXLINE];

  temp = depointer_element(type);

  if( temp == NULL ) {
    warn("Unable to depointer %s",type);
    return NULL;
  }

  if( strchr(temp,'*') != NULL)
    return stringalloc("ckfree");

  if( strstr(temp,"struct") != NULL)
    return NULL;
	
  if( (runner=strtok(temp," *\t\n")) == NULL)
    return NULL;

  if( is_simple_type_string(temp) == TRUE)
    return stringalloc("ckfree");

  strcpy(buffer,"free_");
  strcat(buffer,runner);

  free(temp);

  return stringalloc(buffer);
}


# line 247 "wisec.dy"
void write_hard_link_function(DYNFILE * dfp,StructHolder * sh)
{
  FuncInfo * fi;
  ArgInfo * ai;

  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"hard_link_%s",sh->name);
  add_line_to_Ftext(fi->ft,"Bumps up the reference count of the object");
  add_line_to_Ftext(fi->ft,"Meaning that multiple pointers can 'own' it");
  
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"obj");
  ai->desc=stringalloc("Object to be hard linked");


  start_function_FuncInfo(fi,dfp,"%s * hard_link_%s(%s * obj)",sh->name,sh->name,sh->name);

  expr(dfp,"if( obj == NULL ) ");
  startbrace(dfp);
  warn_expr(dfp,"Trying to hard link to a %s object: passed a NULL object",sh->name);
  expr(dfp,"return NULL;");
  closebrace(dfp);
  expr(dfp,"obj->dynamite_hard_link++;");
  expr(dfp,"return obj;");
  close_function(dfp);
  add_break(dfp);
}

# line 273 "wisec.dy"
void write_free_function(DYNFILE * dfp,StructHolder * sh)
{
  FuncInfo * fi;
  ArgInfo * ai;
  StructElement * temp;
  register int i;
  register int islist=FALSE;
  register int ismat =FALSE;
  char * run;
  char * run2;
  boolean  is_hard_link = FALSE;

  is_hard_link = dfp->should_hard_link;

  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"free_%s",sh->name);
  add_line_to_Ftext(fi->ft,"Free Function: removes the memory held by obj");
  add_line_to_Ftext(fi->ft,"Will chain up to owned members and clear all lists");
  
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"obj");
  ai->desc=stringalloc("Object that is free'd");


  start_function_FuncInfo(fi,dfp,"%s * free_%s(%s * obj)",sh->name,sh->name,sh->name);
  expr(dfp,"int return_early = 0;");
  if( (islist=is_listfunc(sh)) == TRUE)
    {
      expr(dfp,"int i;");
    }
  else if ( (ismat=is_matrixfunc(sh)) == TRUE)
    {
      expr(dfp,"int i;");
      expr(dfp,"int j;");
    }

  add_break(dfp);

  expr(dfp,"if( obj == NULL)");
  startbrace(dfp);
  warn_expr(dfp,"Attempting to free a NULL pointer to a %s obj. Should be trappable",sh->name);	
  expr(dfp,"return NULL;");
  closebrace(dfp);
  add_break(dfp);


  if( is_hard_link == TRUE ) {
    macro(dfp,"#ifdef PTHREAD");
    expr(dfp,"assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0);");
    macro(dfp,"#endif");
    expr(dfp,"if( obj->dynamite_hard_link > 1) ");
    startbrace(dfp);
    expr(dfp,"return_early = 1;");
    expr(dfp,"obj->dynamite_hard_link--;");
    closebrace(dfp);

    macro(dfp,"#ifdef PTHREAD");
    expr(dfp,"assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);");
    macro(dfp,"#endif");
    expr(dfp,"if( return_early == 1)");
    hang_expr(dfp,"return NULL;");

  }

  for(i=0;i<sh->len;i++)
    {
      temp = sh->el[i];
      if( temp->islist == TRUE)
	{
	  run = depointer_element(temp->element_type);
	  if((run2 = free_function(run)) == NULL)
	    {
	      log_full_error(WARNING,0,"Unable to make function for %s (%s)",temp->name,temp->element_type);
	      add_block_comment(dfp,"Unable to make function for %s",temp->name);
	      continue;
	    }

	  expr(dfp,"if( obj->%s != NULL)",temp->name);
	  startbrace(dfp);
	  expr(dfp,"for(i=0;i<obj->%slen;i++)",CKN(temp->len_append));
	  startbrace(dfp);
	  expr(dfp,"if( obj->%s[i] != NULL)",temp->name);
	  hang_expr(dfp,"%s(obj->%s[i]);",run2,temp->name);
	  ckfree(run);
	  ckfree(run2);
	  closebrace(dfp);
	  expr(dfp,"ckfree(obj->%s);",temp->name);
	  closebrace(dfp);
	}
      else if( temp->ismatrix == TRUE)
	{
	  expr(dfp,"if( obj->%s != NULL)",temp->name);
	  startbrace(dfp);
	  expr(dfp,"for(i=0;i<obj->leni;i++)");
	  startbrace(dfp);
	  expr(dfp,"if( obj->%s[i] != NULL)",temp->name);
	  hang_expr(dfp,"ckfree(obj->%s[i]);",temp->name);
	  closebrace(dfp);
	  expr(dfp,"ckfree(obj->%s)",temp->name);
	  closebrace(dfp);
	}
      else if( temp->islinked == TRUE)
	{
	  add_block_comment(dfp,"obj->%s is linked in",temp->name);
	}
      else if( temp->isfunc == TRUE)
	{
	  add_block_comment(dfp,"obj->%s is a function pointer",temp->name);
	}
      else	{
		  if( strchr(temp->element_type,'*') == NULL)
		    continue; 
		  if( (run=free_function(temp->element_type)) == NULL)
		    {
		      add_block_comment(dfp,"Unable to make free function for obj->%s",temp->name);
		      continue;
		    }
		  expr(dfp,"if( obj->%s != NULL)",temp->name);
		  hang_expr(dfp,"%s(obj->%s)",run,temp->name);
		  ckfree(run);
		}
    }
  add_break(dfp);

  expr(dfp,"ckfree(obj);");


  expr(dfp,"return NULL;");
  close_function(dfp);

  add_break(dfp);
}


# line 405 "wisec.dy"
boolean is_simple_type(StructElement * se)
{
  if( se->isfunc == TRUE )
    return FALSE;
  if( strchr(se->element_type,'*') != NULL )
    return FALSE;

  return is_simple_type_string(se->element_type);
}

# line 415 "wisec.dy"
boolean is_simple_type_string(char * str)
{
  if( strcmp(str,"int") == 0)
    return TRUE;
  if( strcmp(str,"char") == 0)
    return TRUE;
  if( strcmp(str,"float") == 0)
    return TRUE;
  if( strcmp(str,"double") == 0)
    return TRUE;
  if( strcmp(str,"boolean") == 0)
    return TRUE;
	
  return FALSE;
}


# line 432 "wisec.dy"
boolean is_string_type(StructElement * se)
{
  if( se->isfunc == TRUE )
    return FALSE;

  if( strcmp(se->element_type,"char *") == 0)
    return TRUE;

  return FALSE;
}

# line 443 "wisec.dy"
boolean is_binary_dump_type(StructElement * se)
{
  if( se->isfunc == TRUE)
    return FALSE;

  if( strchr(se->element_type,'*') == NULL )
    return FALSE;
  if( strstr(se->element_type,"**") != NULL)
    return FALSE;

  return TRUE;
}

# line 456 "wisec.dy"
void write_binary_dump_function(DYNFILE * dfp,StructHolder * sh)
{
  register int i;
	
  start_function(dfp,"boolean binary_dump_%s(%s * obj,FILE * ofp)",sh->name,sh->name);
  if( is_listfunc(sh) == TRUE )
    {
      expr(dfp,"register int i;");
    }
  else if ( is_matrixfunc(sh) == TRUE )
    {
      expr(dfp,"register int i;");
      expr(dfp,"register int j;");
    }

  expr(dfp,"char buffer[MAXBINARYDUMP]");
  expr(dfp,"int strl;");
  add_break(dfp);
	
  /*** first dump tag line ***/

  expr(dfp,"memset(buffer,0,MAXBINARYDUMP)");
  expr(dfp,"strcpy(buffer,\"DYN %s\")",sh->name);
  expr(dfp,"buffer[3] = '\\0';");
  expr(dfp,"buffer[63] = '\\0';");
  expr(dfp,"buffer[64] = '\\0';");
	
  /*** dump ***/

  expr(dfp,"fwrite(buffer,sizeof(char),64,ofp)");

  add_break(dfp);
  add_block_comment(dfp,"written block tag, write out elements");
  add_break(dfp);
	
  /*** ok - now chain to elements ***/

  for(i=0;i<sh->len;i++)
    {
      auto char * run;
      auto char * run2;
      auto StructElement * temp;

      temp = sh->el[i];

      if( temp->isfunc == TRUE || temp->islinked == TRUE || temp->islocal == TRUE)
	continue;

		
      if( temp->islist == TRUE )
	{
	  run = depointer_element(temp->element_type);
	  run2 = depointer_element(run);

	  expr(dfp,"fwrite(&obj->%slen,sizeof(int),1,ofp)",CKN(temp->len_append));
	  expr(dfp,"for(i=0;i<obj->%slen;i++)",CKN(temp->len_append));
	  startbrace(dfp);
	  expr(dfp,"binary_dump_%s(obj->%s[i],ofp)",run2,temp->name);
	  closebrace(dfp);
	  ckfree(run);
	  ckfree(run2);
	  continue;
	}	

      if( is_simple_type(temp) == TRUE )
	{
	  expr(dfp,"fwrite(&obj->%s,sizeof(%s),1,ofp)",temp->name,temp->element_type);
	}
      else if( is_string_type(temp) == TRUE )
	{
	  expr(dfp,"strl = strlen(obj->%s)",temp->name);
	  expr(dfp,"fwrite(&strl,sizeof(int),1,ofp);");
	  expr(dfp,"fwrite(obj->%s,sizeof(char),strlen(obj->%s),ofp)",temp->name,temp->name);
	}
      else if( is_binary_dump_type(temp) == TRUE )
	{
	  run = depointer_element(temp->element_type);
	  /** call binary dump func **/
	  expr(dfp,"if( obj->%s == NULL )",temp->name);
	  hang_expr(dfp,"strl = 0;");
	  expr(dfp,"else strl = 1;");
	  expr(dfp,"fwrite(&strl,sizeof(int),1,ofp)");
	  expr(dfp,"if( obj->%s != NULL )",temp->name);
	  hang_expr(dfp,"binary_dump_%s(obj->%s,ofp)",run,temp->name);
	  ckfree(run);
	}
      else	{
		  warn("Unable to handle binary dump for %s element of %s",temp->name,sh->name);
		}
    }

  expr(dfp,"return TRUE");
  close_function(dfp);
  add_break(dfp);
}

# line 552 "wisec.dy"
void write_binary_read_function(DYNFILE * dfp,StructHolder * sh)
{
  register int i;
  char * lenapp;

	
	
  start_function(dfp,"%s * binary_read_%s(FILE * ifp)",sh->name,sh->name);
  if( is_listfunc(sh) == TRUE )
    {
      expr(dfp,"register int i;");
      expr(dfp,"int count");
    }
  else if ( is_matrixfunc(sh) == TRUE )
    {
      expr(dfp,"register int i;");
      expr(dfp,"register int j;");
    }

  expr(dfp,"char buffer[MAXBINARYDUMP]");
  expr(dfp,"int strl;");
  expr(dfp,"%s * obj",sh->name);
  add_break(dfp);
	
  /*** first read dump tag line ***/

  expr(dfp,"fread(buffer,sizeof(char),64,ifp)");
  expr(dfp,"if( strcmp(buffer,\"DYN\") != 0 )");
  startbrace(dfp);
  warn_expr(dfp,"Unable to read DYN tag in binary read for %s",sh->name);
  expr(dfp,"return NULL;");
  closebrace(dfp);

  expr(dfp,"if( strcmp(buffer+4,\"%s\") != 0 )",sh->name);
  startbrace(dfp);
  expr(dfp,"warn(\"Unable to pick up %s tag in reading binary dump, got %%s\",buffer+4)",sh->name);
  expr(dfp,"return NULL");
  closebrace(dfp);

	


  add_break(dfp);
  add_block_comment(dfp,"read block tag,time to read in elements");
  add_break(dfp);

  if( is_listfunc(sh) == TRUE)
    expr(dfp,"obj = %s_alloc_std()",sh->name);
  else	expr(dfp,"obj = %s_alloc()",sh->name);
  expr(dfp,"if( obj == NULL )");
  startbrace(dfp);
  warn_expr(dfp,"Memory allocation problem in %s binary read",sh->name);
  expr(dfp,"return NULL;");
  closebrace(dfp);
	
  /*** ok - now chain to elements ***/

  for(i=0;i<sh->len;i++)
    {
      auto char * run;
      auto char * run2;
      auto StructElement * temp;

      temp = sh->el[i];

      if( temp->isfunc == TRUE || temp->islinked == TRUE || temp->islocal == TRUE)
	continue;

		
      if( temp->islist == TRUE )
	{
	  if( temp->len_append == NULL)
	    lenapp = "";
	  else	lenapp = temp->len_append;

	  run = depointer_element(temp->element_type);
	  run2 = depointer_element(run);

	  expr(dfp,"fread(&count,sizeof(int),1,ifp)");
	  expr(dfp,"for(;count > 0 ;count--)",temp->len_append);
	  startbrace(dfp);
	  expr(dfp,"add_%s%s(obj,binary_read_%s(ifp))",lenapp,temp->name,run2);
	  closebrace(dfp);
	  ckfree(run);
	  ckfree(run2);
	  continue;
	}	

      if( is_simple_type(temp) == TRUE )
	{
	  expr(dfp,"fread(&obj->%s,sizeof(%s),1,ifp)",temp->name,temp->element_type);
	}
      else if( is_string_type(temp) == TRUE )
	{
	  expr(dfp,"memset(buffer,0,MAXBINARYDUMP)");
	  expr(dfp,"fread(&strl,sizeof(int),1,ifp)");
	  expr(dfp,"fread(buffer,sizeof(char),strl,ifp)");
	  expr(dfp,"obj->%s = stringalloc(buffer)",temp->name,temp->name);
	}
      else if( is_binary_dump_type(temp) == TRUE )
	{
	  run = depointer_element(temp->element_type);
	  /** call binary dump func **/
	  expr(dfp,"fread(&strl,sizeof(int),1,ifp)");
	  expr(dfp,"if( strl == 1 )");
	  hang_expr(dfp,"obj->%s = binary_read_%s(ifp)",temp->name,run);
	  ckfree(run);
	}
      else	{
		  warn("Unable to handle binary dump for %s element of %s",temp->name,sh->name);
		}
    }

  expr(dfp,"return obj");
  close_function(dfp);
  add_break(dfp);
}

# line 670 "wisec.dy"
void write_expand_matrix_function(DYNFILE * dfp,StructHolder * sh)
{
  FuncInfo * fi;

  register int i;
  char * runner;
  char * run2;
  StructElement * temp;


  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"expand_%s",sh->name);
  add_line_to_Ftext(fi->ft,"Expands matrix. Rarely used");
  

  start_function_FuncInfo(fi,dfp,"boolean expand_%s(%s * obj,int leni,int lenj)",sh->name,sh->name);

  for(i=0;i<sh->len;i++)
    if(sh->el[i]->ismatrix == TRUE)
      break;

  if( i >= sh->len)
    return;

  temp=sh->el[i];
	


  expr(dfp,"int i");
  expr(dfp,"int actualj");

  add_break(dfp);

  expr(dfp,"if( obj == NULL) ");
  startbrace(dfp);
  expr(dfp,"warn(\"Trying to expand a %s but is NULL!\");",sh->name);
  expr(dfp,"return FALSE;");
  closebrace(dfp);
	
  add_break(dfp);

  expr(dfp,"if( leni <= obj->maxleni && lenj <= obj->maxlenj)");
  hang_expr(dfp,"return TRUE");

  add_break(dfp);

  expr(dfp,"if( obj->maxleni < leni )");
  startbrace(dfp);
  if( (runner = depointer_element(temp->element_type)) == NULL)
    log_full_error(WARNING,0,"Oooops");
  if( (run2 = depointer_element(runner)) == NULL)
    log_full_error(WARNING,0,"Oooops");
	
  expr(dfp,"if( (obj->%s=(%s) ckrealloc (obj->%s,sizeof(%s)*leni)) == NULL)",temp->name,temp->element_type,temp->name,runner);
  hang_expr(dfp,"return FALSE;");
  expr(dfp,"obj->maxleni=obj->leni=leni");
  closebrace(dfp);	

  expr(dfp,"if( lenj > obj->maxlenj )");
  hang_expr(dfp,"actualj = lenj");
  expr(dfp,"else actualj = obj->maxlenj;");

  expr(dfp,"for(i=0;i<obj->leni;i++)");
  startbrace(dfp);
  expr(dfp,"if((obj->%s[i] = (%s) realloc (obj->%s[i],sizeof(%s) * actualj)) == NULL)",temp->name,runner,temp->name,run2) ;
  hang_expr(dfp,"return FALSE;");
	
  closebrace(dfp);
	
  expr(dfp,"return TRUE");

  close_function(dfp);

  ckfree(runner);
  ckfree(run2);

  add_break(dfp);

  return; 
}

# line 750 "wisec.dy"
void write_alloc_matrix_function(DYNFILE * dfp,StructHolder * sh)
{
  FuncInfo * fi;
  ArgInfo * ai;

  register int i;
  char * run2;
  char * runner;


  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"%s_alloc_matrix",sh->name);
  add_line_to_Ftext(fi->ft,"Allocates structure and matrix");
  
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"leni");
  ai->desc=stringalloc("Length of first dimension of matrix");
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"lenj");
  ai->desc=stringalloc("Length of second dimension of matrix");


  start_function_FuncInfo(fi,dfp,"%s * %s_alloc_matrix(int leni,int lenj)",sh->name,sh->name);
	

  expr(dfp,"%s * out",sh->name);
  add_end_comment(dfp,"out is exported         ");
  expr(dfp,"register int i");
  add_end_comment(dfp,"for stepping down matrix");	
  expr(dfp,"register int j");
  add_end_comment(dfp,"for stepping across matrix");

  add_break(dfp);

  add_block_comment(dfp,"Call alloc function, return NULL if NULL");
  expr(dfp,"if((out = %s_alloc()) == NULL)",sh->name);
  hang_expr(dfp,"return NULL");
	
  add_break(dfp);


  for(i=0;i<sh->len;i++)
    if( sh->el[i]->ismatrix == TRUE)
      {
	if( (runner=depointer_element(sh->el[i]->element_type)) == NULL)
	  continue;
	if( (run2=depointer_element(runner)) == NULL)
	  continue;

	add_block_comment(dfp,"Allocate memory for %s ",sh->el[i]->name);

	expr(dfp,"if((out->%s = (%s) ckcalloc (leni,sizeof(%s))) == NULL)",sh->el[i]->name,sh->el[i]->element_type,runner);
	startbrace(dfp);
	expr(dfp,"warn(\"Memory allocation problem in matrix for %s %s, first pointer set\");",sh->name,sh->el[i]->name);
	expr(dfp,"ckfree(out)");
	expr(dfp,"return NULL");

	closebrace(dfp);
	add_break(dfp);
	add_block_comment(dfp,"Add NULL to all matrix pointers so free can be called");
	expr(dfp,"for(i=0;i<leni;i++)");
	hang_expr(dfp,"out->%s[i] = NULL",sh->el[i]->name);
	add_break(dfp);

	add_block_comment(dfp,"Allocate each matrix row");

	expr(dfp,"for(i=0;i<leni;i++)");
	startbrace(dfp);
	expr(dfp,"out->%s[i] = (%s) ckcalloc (lenj,sizeof(%s))",sh->el[i]->name,runner,run2);
	expr(dfp,"if( out->%s[i] == NULL)",sh->el[i]->name);
	startbrace(dfp);
	expr(dfp,"warn(\"Failed alloc on %%d, calling free and returning NULL\",i);");
	expr(dfp,"free_%s(out)",sh->name);
	expr(dfp,"return NULL");
	closebrace(dfp);
	closebrace(dfp);

	add_break(dfp);
	if( sh->el[i]->def != NULL )
	  {
	    expr(dfp,"for(i=0;i<leni;i++)");
	    startbrace(dfp);
	    expr(dfp,"for(j=0;j<lenj;j++)");
	    hang_expr(dfp,"out->%s[i][j] = %s",sh->el[i]->name,sh->el[i]->def);
	    closebrace(dfp);
	    add_break(dfp);
	  }
	expr(dfp,"out->leni=out->maxleni=leni");
	expr(dfp,"out->lenj=out->maxlenj=lenj");

	add_break(dfp);
      }

  expr(dfp,"return out");

  close_function(dfp);
  add_break(dfp);
}



# line 848 "wisec.dy"
void write_list_sort_function(DYNFILE * dfp,StructHolder * sh,StructElement * temp)
{
  FuncInfo * fi;
  ArgInfo * ai;
  register char * runner;
  char * listappend;

  listappend=CKN(temp->len_append);

  runner = depointer_element(temp->element_type);

  if( runner == NULL)
    {
      log_full_error(WARNING,0,"Unable to depointer %s",temp->element_type);
      return;
    }


  add_block_comment(dfp,"swap function for qsort function");


  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"swap_%s%s",listappend,sh->name);
  add_line_to_Ftext(fi->ft,"swap function: an internal for qsort_%s%s",listappend,sh->name);
  add_line_to_Ftext(fi->ft,"swaps two positions in the array");
  
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"list");
  ai->desc=stringalloc("List of structures to swap in");
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"i");
  ai->desc=stringalloc("swap position");
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"j");
  ai->desc=stringalloc("swap position");


  start_function_FuncInfo(fi,dfp,"void swap_%s%s(%s list,int i,int j) ",listappend,sh->name,temp->element_type);
  expr(dfp,"%s temp;",runner);

  expr(dfp,"temp=list[i];");
  expr(dfp,"list[i]=list[j];");
  expr(dfp,"list[j]=temp;");

  close_function(dfp);
  add_break(dfp);


  /*** sort function actual type ***/

  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"qsort_%s%s",listappend,sh->name);
  add_line_to_Ftext(fi->ft,"qsort - lifted from K&R ",listappend,sh->name);
  add_line_to_Ftext(fi->ft,"sorts the array using quicksort");
  add_line_to_Ftext(fi->ft,"Probably much better to call sort_%s%s which sorts from start to end",listappend,sh->name);
  
  
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"list");
  ai->desc=stringalloc("List of structures to swap in");
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"left");
  ai->desc=stringalloc("left position");
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"right");
  ai->desc=stringalloc("right position");
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"comp");
  ai->desc=stringalloc("Function which returns -1 or 1 to sort on");

  start_function_FuncInfo(fi,dfp,"void qsort_%s%s(%s list,int left,int right,int (*comp)(%s ,%s ))",listappend,sh->name,temp->element_type,runner,runner);
  expr(dfp,"int i,last;");
  expr(dfp,"if( left >= right )");
  hang_expr(dfp,"return;");
  add_break(dfp);

  expr(dfp,"swap_%s%s(list,left,(left+right)/2);",listappend,sh->name);
  expr(dfp,"last = left;");
  expr(dfp,"for ( i=left+1; i <= right;i++)");
  startbrace(dfp);
  expr(dfp,"if( (*comp)(list[i],list[left]) < 0)");
  hang_expr(dfp,"swap_%s%s (list,++last,i);",listappend,sh->name);
  closebrace(dfp);

  expr(dfp,"swap_%s%s (list,left,last);",listappend,sh->name);
  expr(dfp,"qsort_%s%s(list,left,last-1,comp);",listappend,sh->name);
  expr(dfp,"qsort_%s%s(list,last+1,right,comp);",listappend,sh->name);
  close_function(dfp);
  add_break(dfp);


  /*** actual function to be called ***/



  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"sort_%s%s",listappend,sh->name);
  add_line_to_Ftext(fi->ft,"sorts from start to end using comp ",listappend,sh->name);
  add_line_to_Ftext(fi->ft,"sorts the array using quicksort by calling qsort_%s%s",listappend,sh->name);
  
  
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"obj");
  ai->desc=stringalloc("Object containing list");
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"comp");
  ai->desc=stringalloc("Function which returns -1 or 1 to sort on");

  start_function_FuncInfo(fi,dfp,"void sort_%s%s(%s * obj,int (*comp)(%s, %s))",listappend,sh->name,sh->name,runner,runner);
  expr(dfp,"qsort_%s%s(obj->%s,0,obj->%slen-1,comp);",listappend,sh->name,temp->name,listappend);
  expr(dfp,"return;");
  close_function(dfp);

  add_break(dfp);
}	

# line 952 "wisec.dy"
void write_list_flush_function(DYNFILE * dfp,StructHolder * sh,StructElement * temp)
{
  FuncInfo * fi;
  ArgInfo * ai;
  register char * runner;
  char * listapp;
  char * run2;

  listapp = CKN(temp->len_append);

  runner = depointer_element(temp->element_type);
  run2   = free_function(runner);


  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"flush_%s%s",listapp,sh->name);
  add_line_to_Ftext(fi->ft,"Frees the list elements, sets length to 0");
  add_line_to_Ftext(fi->ft,"If you want to save some elements, use hard_link_xxx");
  add_line_to_Ftext(fi->ft,"to protect them from being actually destroyed in the free");
  fi->simple= stringallocf("flush_%s",temp->name);

  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"obj");
  ai->desc=stringalloc("Object which contains the list ");

  start_function_FuncInfo(fi,dfp,"int flush_%s%s(%s * obj)",listapp,sh->name,sh->name);
  expr(dfp,"int i;");
  
  add_break(dfp);
  expr(dfp,"for(i=0;i<obj->%slen;i++)",listapp);
  startbrace_tag(dfp,"for i over list length");
  expr(dfp,"if( obj->%s[i] != NULL)",temp->name);
  startbrace(dfp);
  expr(dfp,"%s(obj->%s[i]);",run2,temp->name);
  expr(dfp,"obj->%s[i] = NULL;",temp->name);
  closebrace(dfp);
  closebrace(dfp);
  add_break(dfp);
  expr(dfp,"obj->%slen = 0;",listapp);

  expr(dfp,"return i;");
  close_function(dfp);
  add_break(dfp);
}
  

# line 996 "wisec.dy"
void write_list_add_function(DYNFILE * dfp,StructHolder * sh,StructElement * temp)
{
  FuncInfo * fi;
  ArgInfo * ai;
  register char * runner;
  char * listappend;


  listappend = CKN(temp->len_append);


  runner = depointer_element(temp->element_type);

  add_block_comment(dfp,"will expand function if necessary");


  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"add_%s%s",listappend,sh->name);
  add_line_to_Ftext(fi->ft,"Adds another object to the list. It will expand the list if necessary");
  fi->simple = stringallocf("add_%s",temp->name);

  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"obj");
  ai->desc=stringalloc("Object which contains the list");
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"add");
  ai->desc=stringalloc("Object to add to the list");
  ai->argtype = ARGTYPE_OWNER;
	
  start_function_FuncInfo(fi,dfp,"boolean add_%s%s(%s * obj,%s add)",listappend,sh->name,sh->name,runner);
  expr(dfp,"if( obj->%slen >= obj->%smaxlen)",listappend,listappend);
  startbrace(dfp);
  expr(dfp,"if( expand_%s%s(obj,obj->%slen + %sLISTLENGTH) == FALSE)",listappend,sh->name,listappend,sh->name);
  hang_expr(dfp,"return FALSE;");
  closebrace(dfp);
  add_break(dfp);
  expr(dfp,"obj->%s[obj->%slen++]=add;",temp->name,listappend);
  expr(dfp,"return TRUE;");

  close_function(dfp);
  add_break(dfp);
}


# line 1037 "wisec.dy"
void write_list_expand_function(DYNFILE * dfp,StructHolder * sh,StructElement * temp)
{
  FuncInfo * fi;
  ArgInfo * ai;
  register char * runner;

  char * listappend;

  listappend = CKN(temp->len_append);

  runner = depointer_element(temp->element_type);

  if( runner == NULL)
    {
      log_full_error(WARNING,0,"Unable to depointer %s",temp->element_type);
      return;
    }

  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"expand_%s%s",listappend,sh->name);
  add_line_to_Ftext(fi->ft,"Really an internal function for add_%s%s",listappend,sh->name);
  
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"obj");
  ai->desc=stringalloc("Object which contains the list");
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"len");
  ai->desc=stringalloc("Length to add one");

  start_function_FuncInfo(fi,dfp,"boolean expand_%s%s(%s * obj,int len)",CKN(listappend),sh->name,sh->name);
	

  add_break(dfp);

  expr(dfp,"if( obj->%smaxlen > obj->%slen ) ",listappend,listappend);
  startbrace(dfp);
  expr(dfp,"warn(\"expand_%s%s called with no need\");",sh->name,CKN(listappend));
  expr(dfp,"return TRUE;");
  closebrace(dfp);


  add_break(dfp);

  expr(dfp,"if( (obj->%s = (%s ) ckrealloc (obj->%s,sizeof(%s)*len)) == NULL) ",
       temp->name,temp->element_type,temp->name,runner);

  startbrace(dfp);
  expr(dfp,"warn(\"ckrealloc failed for expand_%s, returning FALSE\");",sh->name);
  expr(dfp,"return FALSE;");
  closebrace(dfp);

  expr(dfp,"obj->%smaxlen = len;",listappend,listappend);
  expr(dfp,"return TRUE;");

  close_function(dfp);

  add_break(dfp);
}

# line 1093 "wisec.dy"
void write_alloc_std_function(DYNFILE * dfp,StructHolder * sh)
{
  FuncInfo * fi;

  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"%s_alloc_std",sh->name);
  add_line_to_Ftext(fi->ft,"Equivalent to %s_alloc_len(%sLISTLENGTH)",sh->name,sh->name);
  
  start_function_FuncInfo(fi,dfp,"%s * %s_alloc_std(void)",sh->name,sh->name);
  expr(dfp,"return %s_alloc_len(%sLISTLENGTH);",sh->name,sh->name);
  close_function(dfp);
  add_break(dfp);
}

# line 1106 "wisec.dy"
void write_alloc_len_function(DYNFILE * dfp,StructHolder * sh)
{
  FuncInfo * fi;
  ArgInfo * ai;
  register int i;

  char * runner;

  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"%s_alloc_len",sh->name);
  add_line_to_Ftext(fi->ft,"Allocates len length to all lists",sh->name,sh->name);
  ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"len");
  ai->desc=stringalloc("Length of lists to allocate");



  start_function_FuncInfo(fi,dfp,"%s * %s_alloc_len(int len)",sh->name,sh->name);
  expr(dfp,"%s * out;",sh->name);
  add_end_comment(dfp,"out is exported at the end of function");

  add_break(dfp);

  add_block_comment(dfp,"Call alloc function: return NULL if NULL");
  add_block_comment(dfp,"Warning message alread in alloc function");
  expr(dfp,"if((out = %s_alloc()) == NULL)",sh->name);
  hang_expr(dfp,"return NULL;");

  add_break(dfp);

  add_block_comment(dfp,"Calling ckcalloc for list elements");

  for(i=0;i<sh->len;i++)
    {
      if( sh->el[i]->islist == FALSE)
	continue;
      if( (runner=depointer_element(sh->el[i]->element_type)) ==
	 NULL)
	continue;
      expr(dfp,"if((out->%s = (%s ) ckcalloc (len,sizeof(%s))) == NULL)",sh->el[i]->name,sh->el[i]->element_type,runner);
      free(runner);
      startbrace(dfp);
      expr(dfp,"warn(\"Warning, ckcalloc failed in %s_alloc_len\");",sh->name);
      expr(dfp,"return NULL;");
      closebrace(dfp);
      expr(dfp,"out->%slen = 0;",CKN(sh->el[i]->len_append));
      expr(dfp,"out->%smaxlen = len;",CKN(sh->el[i]->len_append));
      add_break(dfp);
    }

  expr(dfp,"return out;");

  close_function(dfp);

  add_break(dfp);

} 


# line 1163 "wisec.dy"
char * depointer_element(char * element)
{
  char buffer[MAXLINE];
  char * runner;


  for(runner=buffer;*element && *element != '*';)
    *(runner++)=(*(element++));

  if( *element == '\0')
    return NULL;

  element++;			/* depointers it */

  while( *element )
    *(runner++)=(*(element++));
  *(runner)='\0';

  return stringalloc(buffer);
}

# line 1184 "wisec.dy"
boolean isarray(StructElement * se)
{
  if( strchr(se->name,'[') != NULL )
    return TRUE;
  return FALSE;
}

# line 1191 "wisec.dy"
void write_simplealloc_function(DYNFILE * dfp,StructHolder * sh)
{
  FuncInfo * fi;
  register int i;
  boolean  is_hard_link = FALSE;

  is_hard_link = dfp->should_hard_link;

  fi = FuncInfo_named_from_varstr(FI_CALLABLE,"%s_alloc",sh->name);
  add_line_to_Ftext(fi->ft,"Allocates structure: assigns defaults if given ",sh->name,sh->name);
  fi->simple = stringalloc("alloc");

  start_function_FuncInfo(fi,dfp,"%s * %s_alloc(void)",sh->name,sh->name);

  expr(dfp,"%s * out;",sh->name);
  add_end_comment(dfp,"out is exported at end of function");

  add_break(dfp);

  add_block_comment(dfp,"call ckalloc and see if NULL");

  expr(dfp,"if((out=(%s *) ckalloc (sizeof(%s))) == NULL)",
       sh->name,sh->name);

  startbrace(dfp);

  expr(dfp,"warn(\"%s_alloc failed \");", sh->name);

  expr(dfp,"return NULL;");
  add_end_comment(dfp,"calling function should respond!");


  closebrace(dfp);

  if( is_hard_link == TRUE ) {
    expr(dfp,"out->dynamite_hard_link = 1;");
    macro(dfp,"#ifdef PTHREAD");
    expr(dfp,"pthread_mutex_init(&(out->dynamite_mutex),NULL)");
    macro(dfp,"#endif");
  }

  for(i=0;i<sh->len;i++)
    {
      if( sh->el[i]->islinked == TRUE)
	continue;
      else if( isarray(sh->el[i]) == TRUE ) {
	add_block_comment(dfp,"%s is an array: no default possible",sh->el[i]->name);
      }
      else if( sh->el[i]->isfunc == TRUE)
	{
	  expr(dfp,"out->%s = NULL;",sh->el[i]->name);
	}
      else if( sh->el[i]->islist == TRUE)
	{
	  expr(dfp,"out->%s = NULL;",sh->el[i]->name);
	  expr(dfp,"out->%slen = out->%smaxlen = 0;",CKN(sh->el[i]->len_append),CKN(sh->el[i]->len_append));
	}
      else if( sh->el[i]->ismatrix == TRUE)
	{
	  expr(dfp,"out->%s = NULL;",sh->el[i]->name);
	  expr(dfp,"out->leni=out->maxleni=0;");
	  expr(dfp,"out->lenj=out->maxlenj=0;");
	}
      else if( sh->el[i]->def != NULL)
	expr(dfp,"out->%s = %s;",sh->el[i]->name,
	     sh->el[i]->def);
      else if( strcmp(sh->el[i]->element_type,"PSPosition") == 0)
	expr(dfp,"out->%s.x = out->%s.y = 0;",sh->el[i]->name,sh->el[i]->name);
      else	add_block_comment(dfp,"Unable to find def for %s",
				  sh->el[i]->name);
    }

  add_break(dfp);


  expr(dfp,"return out;");

  close_function(dfp);

  add_break(dfp);
}


# line 1274 "wisec.dy"
boolean is_listfunc(StructHolder * sh)
{
  register int i;

  for(i=0;i<sh->len;i++)
    if( sh->el[i]->islist == TRUE)
      return TRUE;

  return FALSE;
}


# line 1286 "wisec.dy"
boolean is_matrixfunc(StructHolder * sh)
{
  register int i;

  for(i=0;i<sh->len;i++)
    if( sh->el[i]->ismatrix == TRUE)
      return TRUE;
  return FALSE;
}

/* IO */

# line 1298 "wisec.dy"
char * def_from_element(char * str)
{
  char * runner;

  if( strchr(str,'*') != NULL)
    return stringalloc("NULL");

  if( (runner=strtok(str," \t*")) == NULL)
    return NULL;

  if( strcmp(str,"int") == 0 || strcmp(str,"float") == 0 || strcmp(str,"double") == 0)
    return stringalloc("0");
  if( strcmp(str,"char") == 0)
    return stringalloc("'u'");
  if( strcmp(str,"boolean") == 0 ) 
    return stringalloc("FALSE");
  if( strcmp(str,"Score") == 0 ) 
    return stringalloc("0");
  if( strcmp(str,"Probability") == 0 ) 
    return stringalloc("0.0");
  

  warn("For element %s - got no good default. Returning 0",str); 
  return stringalloc("0");
}

# line 1324 "wisec.dy"
char * get_word_from_bang(char * str,char ** word)
{

  if( str == NULL || *str == '\0')
    return NULL;
	



  while( *str != '\0'  && *str != '!')
    str++;

  if( *str == '\0')
    {
      *word=NULL;
      return NULL;
    }

  *word=str+1;

  while( *str && !isspace(*str) )
    str++;
  *str='\0';

  str++;


  return str;
}

# line 1354 "wisec.dy"
boolean read_StructHolder_elements(StructHolder * sh,FILE * ifp)
{
  char buffer[MAXLINE];
  char * runner;
  char * comment;
  boolean islink;
  boolean islist;
  boolean isfunc=FALSE;
  boolean ishidden=FALSE;
  boolean ismatrix;
  boolean islocal;
  char * def;
  char * lenapp;
  char * name;
  char ** base;
  char ** brk;
  StructElement * temp;
  ObjectInfo * oi;

  while( get_watched_line(buffer,MAXLINE,ifp) != NULL)
    {
      temp=NULL;
      comment=NULL;
      runner=NULL;
      def=NULL;
      lenapp=NULL;
      isfunc=islink=islist=ismatrix=islocal=FALSE;
      if( strstartcmp(buffer,"%info") == 0 ) {
	oi = read_ObjectInfo_line_func(buffer,MAXLINE,ifp,get_watched_line);
	if( oi == NULL ) {
	  warn("In Object %s, bad info!",sh->name);
	  continue;
	}
	sh->oi = oi;
	continue;
      }

      if( (runner=strstr(buffer,"//")) != NULL)
	{
	  *(runner++)='\0';
	  while( *runner == '/')
	    runner++;
	  comment = stringalloc (runner);
	  /*	striptoprint(comment); */
	}

      if((runner=strchr(buffer,'!')) != NULL)
	{
	  *(runner-1)='\0';
	  base=brk=breakstring(runner,spacestr);
	  for(;*brk != NULL;brk++)
	    {
	      if( strstartcmp(*brk,"!link") == 0)
		islink=TRUE;
	      else if( strstartcmp(*brk,"!matrix") == 0)
		ismatrix=TRUE;
	      else if( strstartcmp(*brk,"!func") == 0)
		isfunc=TRUE;
	      else if( strstartcmp(*brk,"!list") == 0)
		islist=TRUE;
	      else if( strstartcmp(*brk,"!local") == 0)
		islocal=TRUE;
	      else if( strstartcmp(*brk,"!hidden") == 0)
		ishidden=TRUE;
	      else if( strstartcmp(*brk,"!def") == 0)
		{
		  def=string_from_quoted_equality(*brk);
		}
	      else if( strstartcmp(*brk,"!len") == 0)
		{
		  lenapp=string_from_quoted_equality(*brk);
		}	
	    }
	  ckfree(base);
	}


      runner=buffer+strlen(buffer);

      while( runner != NULL && !isalnum(*runner) && *runner != ']' && *runner != ')' )
	{
	  if( runner == buffer)
	    runner=NULL;
	  else	*(runner--) = '\0';
	}

      if( runner == NULL)
	break;
		

      while( isalnum(*runner) || strchr("*_[])(",*runner) != NULL )  {
	if( *runner == ')') {
	  while( *runner != '(' )
	    runner--;
	}
	else runner--;
      }

      name=runner+1;
      *runner='\0';		/* better not be flush! */

      /*	fprintf(stderr,"Adding %s and %s\n",name,buffer); */


      if((temp=basic_add_StructElement( sh,name,buffer)) == NULL)
	fatal("Memeory problem... in adding struct element. Ugh!");

      if( comment != NULL)
	temp->comment=comment;
      if( def != NULL)
	temp->def=def;
      else	temp->def=def_from_element(temp->element_type);

      if( lenapp != NULL)
	temp->len_append=lenapp;

      if( islocal == TRUE)
	temp->islocal = TRUE;

      if( islist == TRUE) {
	temp->islist=TRUE;
	
	/*** make list components ***/

	/*** disabled: handled in write 
	sprintf(buf2,"%slen",temp->len_append == NULL ? "" : temp->len_append );
	basic_add_StructElement(sh,buf2,"int");
	
	sprintf(buf2,"%smaxlen",temp->len_append == NULL ? "" : temp->len_append );
	basic_add_StructElement(sh,buf2,"int");
	***/

      }

      if( islink == TRUE)
	temp->islinked = TRUE;
      if( ismatrix == TRUE) {
	temp->ismatrix = TRUE;

	/*** make matrix components ***/

	/*** disabled: handled in write 
	sprintf(buf2,"%sleni",temp->len_append == NULL ? "" : temp->len_append);
	basic_add_StructElement(sh,buf2,"int");

	sprintf(buf2,"%smaxleni",temp->len_append == NULL ? "" : temp->len_append);
	basic_add_StructElement(sh,buf2,"int");

	sprintf(buf2,"%slenj",temp->len_append == NULL ? "" : temp->len_append);
	basic_add_StructElement(sh,buf2,"int");

	sprintf(buf2,"%smaxlenj",temp->len_append == NULL ? "" : temp->len_append);
	basic_add_StructElement(sh,buf2,"int");
	***/

      }

      if( ishidden == TRUE ) {
	temp->ishidden = TRUE;
      }

      if( isfunc == TRUE ) { 
	temp->isfunc = TRUE;
	runner = name_from_func_type(temp->name);
	sprintf(buffer,"%s %s",temp->element_type,temp->name);
	temp->element_type = stringalloc(buffer);
	temp->name = runner;
      }

      if( strchr(temp->element_type,'*') != NULL)
	temp->isapointer = TRUE;
    }

  return TRUE;
}

# line 1530 "wisec.dy"
void show_StructElement(StructElement * se,FILE * ofp)
{
  fprintf(ofp,"Name [%s] Type [%s]\n",se->name,se->element_type);
} 

# line 1535 "wisec.dy"
void show_StructHolder(StructHolder   * sh,FILE * ofp)
{
  register int i;

  for(i=0;i<sh->len;i++)
    show_StructElement(sh->el[i],ofp);
}

# line 1543 "wisec.dy"
StructElement * basic_add_StructElement(StructHolder * sh,char * name,char * element)
{
  StructElement * se;

  se = StructElement_alloc();

  /*  fprintf(stderr,"allocating %s %s\n",name,element);
   */

  se->name = stringalloc(name);
  se->element_type = stringalloc(element);

  add_StructHolder(sh,se);

  return se;
}


# line 1556 "wisec.c"
/* Function:  hard_link_StructElement(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StructElement *]
 *
 * Return [UNKN ]  Undocumented return value [StructElement *]
 *
 */
StructElement * hard_link_StructElement(StructElement * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a StructElement object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  StructElement_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StructElement *]
 *
 */
StructElement * StructElement_alloc(void) 
{
    StructElement * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(StructElement *) ckalloc (sizeof(StructElement))) == NULL)  {  
      warn("StructElement_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->element_type = NULL;    
    out->comment = NULL; 
    out->def = NULL; 
    out->isapointer = FALSE; 
    out->islinked = FALSE;   
    out->islist = FALSE; 
    out->ismatrix = FALSE;   
    out->len_append = NULL;  
    out->isfunc = FALSE; 
    out->islocal = FALSE;    
    out->ishidden = FALSE;   


    return out;  
}    


/* Function:  free_StructElement(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StructElement *]
 *
 * Return [UNKN ]  Undocumented return value [StructElement *]
 *
 */
StructElement * free_StructElement(StructElement * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a StructElement obj. Should be trappable"); 
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
    if( obj->element_type != NULL)   
      ckfree(obj->element_type);     
    if( obj->comment != NULL)    
      ckfree(obj->comment);  
    if( obj->def != NULL)    
      ckfree(obj->def);  
    if( obj->len_append != NULL) 
      ckfree(obj->len_append);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_StructHolder(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_StructHolder
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [StructElement **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_StructHolder(StructElement ** list,int i,int j)  
{
    StructElement * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_StructHolder(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_StructHolder which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [StructElement **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_StructHolder(StructElement ** list,int left,int right,int (*comp)(StructElement * ,StructElement * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_StructHolder(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_StructHolder (list,++last,i);   
      }  
    swap_StructHolder (list,left,last);  
    qsort_StructHolder(list,left,last-1,comp);   
    qsort_StructHolder(list,last+1,right,comp);  
}    


/* Function:  sort_StructHolder(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_StructHolder
 *
 *
 * Arg:         obj [UNKN ] Object containing list [StructHolder *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_StructHolder(StructHolder * obj,int (*comp)(StructElement *, StructElement *)) 
{
    qsort_StructHolder(obj->el,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_StructHolder(obj,len)
 *
 * Descrip:    Really an internal function for add_StructHolder
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [StructHolder *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_StructHolder(StructHolder * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_StructHolder called with no need");   
      return TRUE;   
      }  


    if( (obj->el = (StructElement ** ) ckrealloc (obj->el,sizeof(StructElement *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_StructHolder, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_StructHolder(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [StructHolder *]
 * Arg:        add [OWNER] Object to add to the list [StructElement *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_StructHolder(StructHolder * obj,StructElement * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_StructHolder(obj,obj->len + StructHolderLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->el[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_StructHolder(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [StructHolder *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_StructHolder(StructHolder * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->el[i] != NULL)    {  
        free_StructElement(obj->el[i]);  
        obj->el[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  StructHolder_alloc_std(void)
 *
 * Descrip:    Equivalent to StructHolder_alloc_len(StructHolderLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StructHolder *]
 *
 */
StructHolder * StructHolder_alloc_std(void) 
{
    return StructHolder_alloc_len(StructHolderLISTLENGTH);   
}    


/* Function:  StructHolder_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [StructHolder *]
 *
 */
StructHolder * StructHolder_alloc_len(int len) 
{
    StructHolder * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = StructHolder_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->el = (StructElement ** ) ckcalloc (len,sizeof(StructElement *))) == NULL)   {  
      warn("Warning, ckcalloc failed in StructHolder_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_StructHolder(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StructHolder *]
 *
 * Return [UNKN ]  Undocumented return value [StructHolder *]
 *
 */
StructHolder * hard_link_StructHolder(StructHolder * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a StructHolder object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  StructHolder_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StructHolder *]
 *
 */
StructHolder * StructHolder_alloc(void) 
{
    StructHolder * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(StructHolder *) ckalloc (sizeof(StructHolder))) == NULL)    {  
      warn("StructHolder_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->el = NULL;  
    out->len = out->maxlen = 0;  
    out->name = NULL;    
    out->placing_function = NULL;    
    out->oi = NULL;  


    return out;  
}    


/* Function:  free_StructHolder(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StructHolder *]
 *
 * Return [UNKN ]  Undocumented return value [StructHolder *]
 *
 */
StructHolder * free_StructHolder(StructHolder * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a StructHolder obj. Should be trappable");  
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
    if( obj->el != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->el[i] != NULL)  
          free_StructElement(obj->el[i]);    
        }  
      ckfree(obj->el);   
      }  
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->placing_function != NULL)   
      ckfree(obj->placing_function);     
    if( obj->oi != NULL) 
      free_ObjectInfo(obj->oi);  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

#ifdef _cplusplus
extern "C" {
#endif
#include "api.h"

  
# line 42 "api.dy"
boolean write_pfdoc_dynAPI(dynAPI * api,char * package_name,FILE * ofp)
{
  int i,j;

  for(i=0;i<api->len;i++) {
    auto APIObject * obj;
    obj = api->obj[i];
    fprintf(ofp,":object %s\n",obj->name);
    for(j=0;j<obj->len;j++) {
      fprintf(ofp,"\t:func %s\n",obj->member[j]->name);
    } 
    fprintf(ofp,"\t:des %s\n",obj->destructor->name);
    fprintf(ofp,"\n\n");
  }


  for(j=0;j<api->non_len;j++) {
      fprintf(ofp,":func %s\n",api->non_obj[j]->name);
  }

  fprintf(ofp,"\n\n");

  for(i=0;i<api->len;i++) {
    auto APIObject * obj;
    obj = api->obj[i];
    if( obj->info != NULL ) { 
      fprintf(ofp,"\n:define %s\n",obj->name);
      fprintf(ofp,":desc\n");
      dump_Ftext(obj->info->ft,ofp);
      fprintf(ofp,"!desc\n!defined\n\n");
    }

    for(j=0;j<obj->len;j++) {
      write_pfdoc_func_def(obj->member[j]->fi,ofp);
      fprintf(ofp,"\n\n");
    } 
    write_pfdoc_func_def(obj->destructor->fi,ofp);
    fprintf(ofp,"\n\n");
  }

  for(i=0;i<api->non_len;i++) {
    write_pfdoc_func_def(api->non_obj[i]->fi,ofp);
    fprintf(ofp,"\n\n");
  }

  return TRUE;
}

# line 90 "api.dy"
boolean write_pfdoc_func_def(FuncInfo * fi,FILE * ofp)
{
  int i;

  fprintf(ofp,":define %s\n",fi->name);
  fprintf(ofp,":desc\n");
  dump_Ftext(fi->ft,ofp);
  fprintf(ofp,"!desc\n");
  fprintf(ofp,":arg\n");
  for(i=0;i<fi->len;i++) {
    fprintf(ofp,"%s [%s] %s\n",fi->arg[i]->name,fi->arg[i]->type,fi->arg[i]->desc);
  }
  fprintf(ofp,"!arg\n");
  fprintf(ofp,":return\n");
  fprintf(ofp,"%s [%s] %s\n",fi->ret->name,fi->ret->type,fi->ret->desc);
  fprintf(ofp,"!return\n");
  fprintf(ofp,"!define");
  
  return TRUE;
}


# line 112 "api.dy"
boolean write_latex_dynAPI(dynAPI * api,char * module,char * package,FILE * ofp)
{
  int i,j;

  fprintf(ofp,"\\section{%s}\n",module);
  fprintf(ofp,"\\label{module_%s}\n",module);
  if( api->len == 0 ) {
    fprintf(ofp,"This module only contains factory methods\n\n");
  } else {
    fprintf(ofp,"This module contains the following objects\n\n\\begin{itemize}\n");
    for(i=0;i<api->len;i++)
      fprintf(ofp,"\\item \\ref{object_%s} %s\n\n",api->obj[i]->name,api->obj[i]->name);
    if( api->non_len > 0 ) {
      fprintf(ofp,"\\item This module also contains some factory methods\n");
    }

    fprintf(ofp,"\\end{itemize}\n");
  }

  if( api->non_len > 0) {
    fprintf(ofp,"\\subsection{%s factory methods}\n",module);
    for(i=0;i<api->non_len;i++) {
      fprintf(ofp,"\\subsubsection{%s}\n",api->non_obj[i]->name);
      write_latex_APIfunc(api->non_obj[i],package,FALSE,NULL,ofp);
    }
    fprintf(ofp,"\n\n");
  }

  for(i=0;i<api->len;i++) {
    fprintf(ofp,"\\subsection{Object %s}\n\n",api->obj[i]->name);
    fprintf(ofp,"\\label{object_%s}\n\n",api->obj[i]->name);
    fprintf(ofp,"The %s object has the following fields. To see how to access them refer to \\ref{accessing_fields}\n",api->obj[i]->name);
    fprintf(ofp,"\\begin{description}\n");
    for(j=0;j<api->obj[i]->sh->len;j++) {
      fprintf(ofp,"\\item{%s} Type [%s : %s] %s\n\n",api->obj[i]->sh->el[j]->name,api->obj[i]->sh->el[j]->element_type,api->obj[i]->sh->el[j]->islist == TRUE ? "List" : "Scalar",api->obj[i]->sh->el[j]->comment == NULL ? "No documentation" : api->obj[i]->sh->el[j]->comment);
    }
    fprintf(ofp,"\\end{description}\n");

    if( api->obj[i]->info == NULL ) {
      fprintf(ofp,"No documentation for %s\n\n",api->obj[i]->name);
    } else {
      latex_Ftext(api->obj[i]->info->ft,ofp);
      fprintf(ofp,"\n\n");
    }
    
    fprintf(ofp,"Member functions of %s\n\n",api->obj[i]->name);
    for(j=0;j<api->obj[i]->len;j++) {
      if( api->obj[i]->member[j]->fi->is_hand_written == TRUE) {
	fprintf(ofp,"\\subsubsection{%s}\n\n",api->obj[i]->member[j]->name);
	write_latex_APIfunc(api->obj[i]->member[j],package,TRUE,api->obj[i]->name,ofp);
      }
    }
  }

  return TRUE;
}


# line 170 "api.dy"
boolean write_latex_APIfunc(APIfunc * f,char * package,boolean isobj,char * objectname,FILE * ofp)
{
  FuncInfo * fi;
  int i;

  fi = f->fi;
  if( fi == NULL ) {
    warn("Could not make latex non API function due to no info!");
    return FALSE;
  }
  
  fprintf(ofp,"\\begin{description}\n");
  fprintf(ofp,"\\item[External C] {\\tt %s_%s (%s",package,fi->name,fi->len == 0 ? "void" : fi->arg[0]->name);
  if( f->only_C == TRUE ) {
    fprintf(ofp,"(This function is only available in the C api)");
  } else {
    for(i=1;i<fi->len;i++)
      fprintf(ofp,",%s",fi->arg[i]->name);
    fprintf(ofp,")}\n");
    if( isobj == TRUE ) {
      fprintf(ofp,"\\item[Perl] {\\tt &%s::%s::%s (%s",package,objectname,fi->simple == NULL ? fi->name : fi->simple,fi->len == 0 ? "" : fi->arg[0]->name); 
    } else {
      fprintf(ofp,"\\item[Perl] {\\tt &%s::%s (%s",package,fi->simple == NULL ? fi->name : fi->simple,fi->len == 0 ? "" : fi->arg[0]->name);
    }
    
    
    for(i=1;i<fi->len;i++) {
      if( fi->arg[i]->should_NULL == TRUE ) 
	continue; /* skip it out */
      fprintf(ofp,",%s",fi->arg[i]->name);
    }
    fprintf(ofp,")}\n\n");

    if( isobj == TRUE ) {
      if( fi->len == 0 ) {
	warn("Trying to indicate that %s is an object method when it doesn't have an argument",fi->name);
      } else {
	fprintf(ofp,"\\item[Perl-OOP call] {\\tt $obj->%s(%s",fi->simple == NULL ? fi->name : fi->simple,fi->len == 1 ? "" : fi->arg[1]->name);
	for(i=2;i<fi->len;i++) {
	  if( fi->arg[i]->should_NULL == TRUE ) 
	    continue; /* skip it out */
	  fprintf(ofp,",%s",fi->arg[i]->name);
	}
	fprintf(ofp,")}\n\n");
      }
    }
  } /* end of is an non C based function */
  fprintf(ofp,"\\end{description}\n");

  fprintf(ofp,"Arguments\n\\begin{description}\n");

  for(i=0;i<fi->len;i++) {
    auto ArgInfo * ai;
    ai = fi->arg[i];
    if( ai->should_NULL == TRUE ) 
      fprintf(ofp,"\\item[%s] \\em{only for C api} [%s] %s [%s]\n",ai->name,ArgType_to_string(ai->argtype),ai->desc,CKS(ai->type));
    else {  
      fprintf(ofp,"\\item[%s] [%s] %s [%s]\n",ai->name,ArgType_to_string(ai->argtype),ai->desc,CKS(ai->type));
    }
  }
  if( strcmp(fi->ret->type,"void") == 0 ) {
    fprintf(ofp,"\\item[returns] Nothing - no return value\n");
  } else {
    fprintf(ofp,"\\item[returns] [%s] %s [%s]\n",ArgType_to_string(fi->ret->argtype),fi->ret->desc,CKS(fi->ret->type));
  }

  fprintf(ofp,"\\end{description}\n");

  latex_Ftext(fi->ft,ofp);

  return TRUE;
}
  

# line 244 "api.dy"
boolean write_pod_dynAPI(dynAPI * api,char * module,char * package,FILE * ofp)
{
  int i,j;

  fprintf(ofp,"=head1 NAME\n\n%s module - part of the %s package\n\n",module,package);
  if( api->len == 0 ) {
    fprintf(ofp,"=head1 SYNOPSIS\n\nThis module contains helper functions for the %s package\n\n",package);
  } else {
    fprintf(ofp,"=head1 SYNOPSIS\n\nThis module contains the following objects\n\n=over\n\n");
    for(i=0;i<api->len;i++)
      fprintf(ofp,"=item %s\n\n",api->obj[i]->name);
    fprintf(ofp,"\n=back\n\n");
  }
  fprintf(ofp,"=head1 DESCRIPTION\n\n");
  for(i=0;i<api->len;i++) {
    fprintf(ofp,"=head2 Object %s\n\n=over\n\n",api->obj[i]->name);
    for(j=0;j<api->obj[i]->sh->len;j++) {
      fprintf(ofp,"=item %s\n\n Type [%s] %s %s\n\n",api->obj[i]->sh->el[j]->name,api->obj[i]->sh->el[j]->element_type,api->obj[i]->sh->el[j]->islist == TRUE ? "List" : "Scalar",api->obj[i]->sh->el[j]->comment == NULL ? "No documentation" : api->obj[i]->sh->el[j]->comment);
    }
    fprintf(ofp,"\n\n=back\n\n");
    if( api->obj[i]->info == NULL ) {
      fprintf(ofp,"No documentation for %s\n\n",api->obj[i]->name);
    } else {
      dump_Ftext(api->obj[i]->info->ft,ofp);
      fprintf(ofp,"\n\n");
    }
    fprintf(ofp,"=head2 Member functions of %s\n\n",api->obj[i]->name);
    fprintf(ofp,"=over\n\n");
    for(j=0;j<api->obj[i]->len;j++) {
      fprintf(ofp,"=item %s\n\n",api->obj[i]->member[j]->fi->simple == NULL ? api->obj[i]->member[j]->name : api->obj[i]->member[j]->fi->simple);
      write_pod_obj_APIfunc(api->obj[i]->member[j],package,api->obj[i]->name,ofp);
    }
    fprintf(ofp,"=back\n\n");
  }
  if( api->non_len > 0) {
    fprintf(ofp,"=over\n\n");
    for(i=0;i<api->non_len;i++) {
      fprintf(ofp,"=item %s\n\n",api->non_obj[i]->name);
      write_pod_non_APIfunc(api->non_obj[i],package,ofp);
    }
    fprintf(ofp,"=back\n\n");
  }

  return TRUE;
}
  
# line 290 "api.dy"
void write_pod_obj_APIfunc(APIfunc * f,char * package,char * obj,FILE * ofp)
{
  FuncInfo * fi;
  int i;

  fi = f->fi;

  fprintf(ofp,"&%s::%s::%s(%s",package,obj,fi->simple == NULL ? fi->name : fi->simple ,fi->len == 0 ? "void" : fi->arg[0]->name);
  for(i=1;i<fi->len;i++)
    fprintf(ofp,",%s",fi->arg[i]->name);
  fprintf(ofp,")\n\n");
  
  dump_Ftext_pre("  ",fi->ft,ofp);
  
  fprintf(ofp,"\n\n");
  for(i=0;i<fi->len;i++) {
    auto ArgInfo * ai;
    ai = fi->arg[i];
    if( ai->should_NULL == TRUE ) {
      continue;
    }

    fprintf(ofp,"  Argument %-12s [%s] %s [%s]\n",ai->name,ArgType_to_string(ai->argtype),ai->desc,CKS(ai->type));
  }
  fprintf(ofp,"  Return [%s] %s [%s]\n",ArgType_to_string(fi->ret->argtype),fi->ret->desc,CKS(fi->ret->type));

  fprintf(ofp,"\n\n");
}
  
# line 319 "api.dy"
void write_pod_non_APIfunc(APIfunc * f,char * package,FILE * ofp)
{
  FuncInfo * fi;
  int i;

  fi = f->fi;

  fprintf(ofp,"&%s::%s(%s",package,fi->simple == NULL ? fi->name : fi->simple ,fi->len == 0 ? "void" : fi->arg[0]->name);
  for(i=1;i<fi->len;i++)
    fprintf(ofp,",%s",fi->arg[i]->name);
  fprintf(ofp,")\n\n");
  
  dump_Ftext_pre("  ",fi->ft,ofp);
  
  fprintf(ofp,"\n\n");
  for(i=0;i<fi->len;i++) {
    auto ArgInfo * ai;
    ai = fi->arg[i];
    fprintf(ofp,"  Argument %-12s [%s] %s [%s]\n",ai->name,ArgType_to_string(ai->argtype),ai->desc,CKS(ai->type));
  }
  fprintf(ofp,"  Return [%s] %s [%s]\n",ArgType_to_string(fi->ret->argtype),fi->ret->desc,CKS(fi->ret->type));

  fprintf(ofp,"\n\n");
}
    
  
  
  
  
  
# line 349 "api.dy"
boolean write_type_C_dynAPI(dynAPI * api,char * package_name,FILE * ofp)
{
  int i;
  
  for(i=0;i<api->len;i++) {
    auto APIObject * obj;
    obj = api->obj[i];
  
    fprintf(ofp,"typedef struct %s%s %s%s;\n\n",package_name,obj->name,package_name,obj->name);
  }

  return TRUE;
}

# line 363 "api.dy"
boolean write_C_dynAPI(dynAPI * api,char * package_name,FILE * ofp)
{
  int i,j;

  for(i=0;i<api->len;i++) {
    auto APIObject * obj;
    obj = api->obj[i];
    fprintf(ofp,"\n\n/* Functions that create, manipulate or act on %s\n *\n",obj->name);
    for(j=0;j<obj->len;j++) {
      fprintf(ofp," * %s%s\n",package_name,obj->member[j]->name);
    }
    fprintf(ofp," * %s%s [destructor]\n",package_name,obj->destructor->name);
    fprintf(ofp," *\n */\n\n");
  }
  if( api->non_len > 0 ) { 
    fprintf(ofp,"\n\n/* Helper functions in the module\n *\n");
    for(i=0;i<api->non_len;i++) {
      fprintf(ofp," * %s%s\n",package_name,api->non_obj[i]->name);
    }
    fprintf(ofp," *\n\n");
  }

  for(i=0;i<api->len;i++) {
    auto APIObject * obj;
    obj = api->obj[i];
    fprintf(ofp,"/* API for object %s */\n",obj->name);

    for(j=0;j<obj->len;j++) {
      write_C_APIfunc(obj->member[j],package_name,ofp);
    }
    fprintf(ofp,"/* This is the destructor function, ie, call this to free object*/\n");
    write_C_APIfunc(obj->destructor,package_name,ofp);
  }
  if( api->non_len > 0 ) 
    fprintf(ofp,"\n\n/* These functions are not associated with an object */\n");
  for(i=0;i<api->non_len;i++) {
    write_C_APIfunc(api->non_obj[i],package_name,ofp);
  }

  return TRUE;
}

# line 405 "api.dy"
boolean write_C_APIfunc(APIfunc * api,char * package_name,FILE * ofp)
{
  int i;
  FuncInfo * fi;

  fi = api->fi;

  fprintf(ofp,"/* Function:  %s%s(%s",package_name,CKS(fi->name),fi->len == 0 ? "void" : fi->arg[0]->name);
  for(i=1;i<fi->len;i++)
    fprintf(ofp,",%s",fi->arg[i]->name);
  fprintf(ofp,")\n *\n");

  show_eddystyle_Ftext(fi->ft,"Descrip:",15,ofp,"No Description");
  fprintf(ofp," *\n");

  for(i=0;i<fi->len;i++)
      fprintf(ofp," * %*s%-12s %s [%s%s]\n",-12,"Arg:",fi->arg[i]->name,fi->arg[i]->desc,is_basic_type_API(fi->arg[i]->type) == TRUE ? "" : package_name,fi->arg[i]->type);
  fprintf(ofp," *\n");
  fprintf(ofp," * Returns %s [%s%s]\n",fi->ret->desc,is_basic_type_API(fi->ret->type) == TRUE ? "" : package_name,fi->ret->type);

  fprintf(ofp," *\n */\n");

  fprintf(ofp,"%s%s ",is_basic_type_API(fi->ret->type) == TRUE ? "" : package_name,fi->ret->type);
  fprintf(ofp,"%s%s(",package_name,fi->name);
  for(i=0;i<fi->len;i++) {
    if( fi->arg[i]->argtype == ARGTYPE_P2FUNC ) {
      fprintf(ofp,"%c%s",i == 0 ? ' ' : ',',fi->arg[i]->func_decl);
    } else {
      fprintf(ofp,"%c%s%s %s",i == 0 ? ' ' : ',',is_basic_type_API(fi->arg[i]->type) == TRUE ? "" : package_name,fi->arg[i]->type,fi->arg[i]->name);
    }
  }

  fprintf(ofp,");\n\n");

  return TRUE;
}

# line 442 "api.dy"
boolean is_membasic_type_API(char * type)
{
  char * temp;

  if( is_basic_type_API(type) == FALSE) 
    return FALSE;

  if( strstartcmp(type,"char") == 0 ) {
    if( (temp=strchr(type,'*')) != NULL ) {
      if( strchr(++temp,'*') != NULL ) {
	warn("Can't cope with char **'s or above!");
	return TRUE;
      } else {
	return FALSE; /* char *'s are not membasic */
      }
    } 
  }

  return TRUE; /* default */
}


  
# line 465 "api.dy"
boolean is_basic_type_API(char * type)
{
  if( strstartcmp(type,"const") == 0 ) 
    return TRUE; /*** oops! ***/

  if( strstartcmp(type,"char") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"FILE") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"double") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"float") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"short") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"long") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"int") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"void") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"boolean") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"Score") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"Probability") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"Bits") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"base") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"codon") == 0 ) 
    return TRUE;
  if( strstartcmp(type,"aa") == 0 ) 
    return TRUE;
  return FALSE;
}



# line 505 "api.dy"
dynAPI * read_dynAPI_line(char * line,FILE * ifp)
{
  char buffer[MAXLINE];
  dynAPI * out;
  APIfunc * func;
  APIObject * obj;
  

  if( strstartcmp(line,"api") != 0 ) {
    warn("In building dynAPI, passed in a non api [%s]",line);
    return NULL;
  }

  out = dynAPI_alloc_std();

  while( get_watched_line(buffer,MAXLINE,ifp) != NULL ) {
    if ( buffer[0] == '#' ) {
      continue;
    }
    if( strstartcmp(buffer,"endapi") == 0 ) {
      break;
    }
    if( strstartcmp(buffer,"end") == 0 ) {
      warn("Got an end line but not endapi in api section. Don't like it [%s]",buffer);
      break;
    }

    if( strstartcmp(buffer,"func") == 0 ) {
      func = APIfunc_from_buffer(buffer);
      if( func == NULL )
	continue;
      add_non_dynAPI(out,func);
    } else if ( strstartcmp(buffer,"object") == 0 ) {
      obj = APIObject_from_line(buffer,MAXLINE,ifp);
      if( obj == NULL ) 
	continue;
      add_dynAPI(out,obj);
    } else {
      warn("In reading an api specification, [%s] not understood\n",buffer);
    }
  }

  return out;

}

# line 551 "api.dy"
boolean write_perl_XS_accessor_functions(APIObject * obj,char * package,FILE * ofp)
{
  int j;
  StructHolder * sh;
  char strippack[64];
  char * listappend;
  char * runner;

  /* HACK coming up! */

  strcpy(strippack,package);
  if( strippack[strlen(strippack)-1] == '_') 
    strippack[strlen(strippack)-1] = '\0';

  sh  = obj->sh;


  for(j=0;j<sh->len;j++) {
    auto StructElement * temp;
      
    temp = sh->el[j];
    if( temp->ishidden == TRUE ) 
      continue; /* skip it! */

    if( temp->islist == TRUE ) {
      listappend=CKN(temp->len_append);
      runner = depointer_element(temp->element_type);
      runner = depointer_element(runner);
      while ( !isalpha(runner[strlen(runner)-1]) ) {
	runner[strlen(runner)-1] = '\0';
      }

      fprintf(ofp,"void\neach_%s(obj)\n",temp->name);
      fprintf(ofp,"\t%s%s * obj\n",package,sh->name);
      fprintf(ofp,"\tPPCODE:\n\tint i=0;\n\tint len;\n\tSV* temp;\n\tlen = %s_length_%s_%s(obj);\n\tfor(i=0;i<len;i++){\n",strippack,temp->name,obj->name);
      fprintf(ofp,"\t  temp = sv_newmortal();\n");
      fprintf(ofp,"\t  sv_setref_pv(temp, \"%s::%s\", (void*) (%shard_link_%s(%s_access_%s_%s(obj,i))));\n",strippack,runner,package,runner,strippack,temp->name,obj->name);
      fprintf(ofp,"\t  XPUSHs(temp);\n");
      fprintf(ofp,"\t  }\n");
      fprintf(ofp,"\tXSRETURN(len);\n\n");
    }
  }

  return TRUE;
}
      

# line 598 "api.dy"
boolean write_dynAPI_accessor_functions(DYNFILE * dfp,dynAPI * api)
{
  int i,j;
  FuncInfo * fi;
  ArgInfo * ai;
  APIfunc * af;
  boolean islist;
  
  char * runner;
  char * listappend;

  /* Ngggg! We have to check if there is a list or not. this was bad API design
     in wisec.dy sometime ago. (like - years ago ) Yuk! yuk! */




  for(i=0;i<api->len;i++) {
    auto APIObject * obj;
    auto StructHolder * sh;


    obj = api->obj[i];
    sh  = obj->sh;

    islist = FALSE;
    for(j=0;j<sh->len;j++) {
      if( sh->el[j]->islist == TRUE ) {
	islist = TRUE;
      }
    }

    /* promote hard link into API */
    af = APIfunc_alloc();
    af->name = stringallocf("hard_link_%s",sh->name);
    add_APIObject(obj,af);

    /* promote alloc into API */

    /* yuk. stupid is list problem! */

    af = APIfunc_alloc();
    af->name = stringallocf("%s_alloc%s",sh->name,islist == TRUE ? "_std" : "");
    af->only_C = TRUE;
    add_APIObject(obj,af);

    for(j=0;j<sh->len;j++) {
      auto StructElement * temp;
      
      temp = sh->el[j];
      if( temp->ishidden == TRUE ) 
	continue; /* skip it! */

      if( temp->islist == TRUE ) {
	listappend=CKN(temp->len_append);
	runner = depointer_element(temp->element_type);


	/* make a "access element" function */

	
	/* build the function documentation. hide it from the C header file */
	fi = FuncInfo_named_from_varstr(FI_INTERNAL,"access_%s_%s",temp->name,sh->name);
	fi->simple = stringalloc(temp->name);
	add_line_to_Ftext(fi->ft,"Access members stored in the %s list",temp->name);
	add_line_to_Ftext(fi->ft,"For use principly by API functions",temp->name);
	ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"obj");
	ai->desc = stringalloc("Object holding the list");
	ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"i");
	ai->desc = stringalloc("Position in the list");

	ai = ArgInfo_alloc();
	ai->name = stringalloc("return");
	ai->desc = stringalloc("Element of the list");
	ai->argtype = ARGTYPE_STATIC;
	fi->ret  = ai;

	/* promote it into API */
	af = APIfunc_alloc();
	af->name = stringalloc(fi->name);
	add_APIObject(obj,af);

	/* build the function */
	
	
	start_function_FuncInfo(fi,dfp,"%s access_%s_%s(%s * obj,int i)",runner,temp->name,sh->name,sh->name);
	expr(dfp,"if( obj == NULL) ");
	startbrace(dfp);
	warn_expr(dfp,"In accessor function %s for object %s, got a NULL object",temp->name,sh->name);
	expr(dfp,"return %s",def_from_element(runner));
	closebrace(dfp);
	
	expr(dfp,"if( obj->%slen <= i )",listappend);
	startbrace(dfp);
	expr(dfp,"warn(\"In accessor function %s for object %s, index %%%%d is greater than list length %%%%d\",i,obj->len);",temp->name,sh->name);
	expr(dfp,"return %s",def_from_element(runner));
	closebrace(dfp);

	expr(dfp,"return obj->%s[i]",temp->name);
	close_function(dfp);
	add_break(dfp);

	/* make a "length" function */

	/* build the function documentation. hide it from the C header file */
	fi = FuncInfo_named_from_varstr(FI_INTERNAL,"length_%s_%s",temp->name,sh->name);
	fi->simple=stringallocf("length_%s",temp->name);
	add_line_to_Ftext(fi->ft,"discover the length of the list",temp->name);
	add_line_to_Ftext(fi->ft,"For use principly by API functions",temp->name);
	ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"obj");
	ai->desc = stringalloc("Object holding the list");

	ai = ArgInfo_alloc();
	ai->name = stringalloc("return");
	ai->desc = stringalloc("length of the list");
	fi->ret = ai;

	/* promote it into API */
	af = APIfunc_alloc();
	af->name = stringalloc(fi->name);
	add_APIObject(obj,af);

	/* build the function */
	
	start_function_FuncInfo(fi,dfp,"int length_%s_%s(%s * obj)",temp->name,sh->name,sh->name);
	expr(dfp,"if( obj == NULL) ");
	startbrace(dfp);
	warn_expr(dfp,"In length function %s for object %s, got a NULL object",temp->name,sh->name);
	expr(dfp,"return -1");
	closebrace(dfp);
	

	expr(dfp,"return obj->%slen",listappend);
	close_function(dfp);
	add_break(dfp);
	
	/* promote flush and add to api */

	af = APIfunc_alloc();
	af->name = stringallocf("flush_%s%s",listappend,sh->name);
	add_APIObject(obj,af);

	af = APIfunc_alloc();
	af->name = stringallocf("add_%s%s",listappend,sh->name);
	add_APIObject(obj,af);


      } else if ( temp->ismatrix == TRUE ) {
	warn("Cannot make matrix accessor functions yet!");
      } else if ( temp->isfunc == TRUE ) {
	warn("Cannot make pointer to functions accessor functions yet (if at all!)");
      } else {

	if ( strchr(temp->name,'[') != NULL  ) {
	  /* should be a compiler option */
	  /* warn("Cannot build accessor for %s, as is an array",temp->name); */
	  continue;
	}
	
	/* simple replace guy */

	/* build the function documentation. hide it from the C header file */
	fi = FuncInfo_named_from_varstr(FI_INTERNAL,"replace_%s_%s",temp->name,sh->name);
	fi->simple = stringallocf("set_%s",temp->name);
	add_line_to_Ftext(fi->ft,"Replace member variable %s",temp->name);
	add_line_to_Ftext(fi->ft,"For use principly by API functions",temp->name);
	ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"obj");
	ai->desc = stringalloc("Object holding the variable");

	ai =  ArgInfo_in_FuncInfo_from_varstr(fi,temp->name);
	ai->argtype = ARGTYPE_OWNER;
	ai->desc = stringalloc("New value of the variable");

	ai = ArgInfo_alloc();
	ai->name = stringalloc("return");
	ai->desc = stringallocf("member variable %s",temp->name);

	ai->argtype = ARGTYPE_STATIC;
	fi->ret =ai;

	/* promote it into API */
	af = APIfunc_alloc();
	af->name = stringalloc(fi->name);
	add_APIObject(obj,af);

	/* build the function */
	
	start_function_FuncInfo(fi,dfp,"boolean replace_%s_%s(%s * obj,%s %s)",temp->name,sh->name,sh->name,temp->element_type,temp->name);
	expr(dfp,"if( obj == NULL) ");
	startbrace(dfp);
	warn_expr(dfp,"In replacement function %s for object %s, got a NULL object",temp->name,sh->name);
	expr(dfp,"return FALSE");
	closebrace(dfp);
	
	expr(dfp,"obj->%s = %s;",temp->name,temp->name);
	expr(dfp,"return TRUE;");
	close_function(dfp);
	add_break(dfp);

	/* simple accessor */

	/* build the function documentation. hide it from the C header file */
	fi = FuncInfo_named_from_varstr(FI_INTERNAL,"access_%s_%s",temp->name,sh->name);
	fi->simple = stringalloc(temp->name);
	add_line_to_Ftext(fi->ft,"Access member variable %s",temp->name);
	add_line_to_Ftext(fi->ft,"For use principly by API functions",temp->name);
	ai =  ArgInfo_in_FuncInfo_from_varstr(fi,"obj");
	ai->desc = stringalloc("Object holding the variable");

	ai = ArgInfo_alloc();
	ai->name = stringalloc("return");
	ai->desc = stringallocf("member variable %s",temp->name);

	ai->argtype = ARGTYPE_STATIC;
	fi->ret =ai;

	/* promote it into API */
	af = APIfunc_alloc();
	af->name = stringalloc(fi->name);
	add_APIObject(obj,af);

	/* build the function */
	
	start_function_FuncInfo(fi,dfp,"%s access_%s_%s(%s * obj)",temp->element_type,temp->name,sh->name,sh->name);
	expr(dfp,"if( obj == NULL) ");
	startbrace(dfp);
	warn_expr(dfp,"In accessor function %s for object %s, got a NULL object",temp->name,sh->name);
	expr(dfp,"return %s",def_from_element(temp->element_type));
	closebrace(dfp);
	
	expr(dfp,"return obj->%s",temp->name);
	close_function(dfp);
	add_break(dfp);
      }
    }
  }

  return TRUE;
}


# line 839 "api.dy"
boolean write_XS_typemap_dynAPI(dynAPI * api,char * package_name,FILE * ofp)
{
  int i;
  char strip_package[64];

  /* HACK coming up! */

  strcpy(strip_package,package_name);
  if( strip_package[strlen(strip_package)-1] == '_') 
    strip_package[strlen(strip_package)-1] = '\0';

  for(i=0;i<api->len;i++) {
    fprintf(ofp,"\nTYPEMAP\n%s%s *    T_%s_%s\n",package_name,api->obj[i]->name,strip_package,api->obj[i]->name);
    fprintf(ofp,"\nINPUT\nT_%s_%s\n\t$var = ($type) (SvROK($arg) == 0 ? NULL : (%s_%s *) SvIV((SV*)SvRV($arg)))\n",strip_package,api->obj[i]->name,strip_package,api->obj[i]->name);
    
    fprintf(ofp,"\nOUTPUT\nT_%s_%s\n\tsv_setref_pv($arg, \"%s::%s\", (void*) $var);\n",strip_package,api->obj[i]->name,strip_package,api->obj[i]->name);
  }
    
  return TRUE;
}

# line 860 "api.dy"
boolean write_XS_dynAPI(dynAPI * api,char * package_name,FILE * ofp)
{
  int i;
  char strip_package[64];

  /* HACK coming up! */

  strcpy(strip_package,package_name);
  if( strip_package[strlen(strip_package)-1] == '_') 
    strip_package[strlen(strip_package)-1] = '\0';

  for(i=0;i<api->len;i++) {
    write_XS_header_APIObject(api,package_name,api->obj[i],ofp);
    write_perl_XS_accessor_functions(api->obj[i],package_name,ofp);
  }

  fprintf(ofp,"\n\nMODULE = %s PACKAGE = %s\n\n",strip_package,strip_package);
  for(i=0;i<api->non_len;i++) {
    write_XS_APIfunc(api->non_obj[i],package_name,ofp);
  }
  
  return TRUE;
}

# line 884 "api.dy"
boolean write_XS_header_APIObject(dynAPI * api,char * package_name,APIObject * obj,FILE * ofp)
{
  int i;
  boolean islist;
  char strip_package[64];

  /* HACK coming up! */

  strcpy(strip_package,package_name);
  if( strip_package[strlen(strip_package)-1] == '_') 
    strip_package[strlen(strip_package)-1] = '\0';
  

  fprintf(ofp,"\n\nMODULE = %s PACKAGE = %s::%s\n\n",strip_package,strip_package,obj->name);

  for(i=0;i<obj->len;i++) {
    write_XS_APIfunc(obj->member[i],package_name,ofp);
  }

  /* now to write a simple constructor */

  /* Ngggg! We have to check if there is a list or not. this was bad API design
     in wisec.dy sometime ago. (like - years ago ) Yuk! yuk! */

  islist = FALSE;
  for(i=0;i<obj->sh->len;i++) {
    if( obj->sh->el[i]->islist == TRUE ) {
      islist = TRUE;
    }
  }

  fprintf(ofp,"\n%s%s *\nnew(class)\n\tchar * class\n\tPPCODE:\n\t%s%s * out;\n\tout = %s%s_alloc%s();\n\tST(0) = sv_newmortal();\n\tsv_setref_pv(ST(0),class,(void*)out);\n\tXSRETURN(1);\n",package_name,obj->name,package_name,obj->name,package_name,obj->name,islist == TRUE ? "_std" : "");
  /* now to write the destructor */

  fprintf(ofp,"\nvoid\nDESTROY(obj)\n\t%s%s * obj\n\tCODE:\n\t%s%s(obj);\n\n",package_name,obj->name,package_name,obj->destructor->name);
  
  return TRUE;
}

# line 923 "api.dy"
boolean write_XS_APIfunc(APIfunc * fu,char * package_name,FILE * ofp)
{
  int j;
  boolean is_basic;
  

 
  /* return statement */
  is_basic = is_basic_type_API(fu->fi->ret->type);
  if( is_basic == TRUE ) 
    fprintf(ofp,"%s\n",fu->fi->ret->type);
  else
    fprintf(ofp,"%s%s\n",(is_basic == TRUE ? "" : package_name),fu->fi->ret->type);
  
  
  /* if we have a simple specification use that - 
     the package will protect us from nasty name clashes!
     */
  if( fu->fi->simple != NULL ) {
    fprintf(ofp,"%s(",fu->fi->simple);
  } else {
    /* function name - we do need the package stuff...*/
    fprintf(ofp,"%s(",fu->name);
  }
  
  for(j=0;j<fu->fi->len;j++) {
    if( fu->fi->arg[j]->should_NULL == TRUE) {
      continue; /* don't put it into the Perl prototype */
    }
    fprintf(ofp,"%s%s",j == 0 ? "" : ",",fu->fi->arg[j]->name);
  }
  fprintf(ofp,")\n");
  
  /* each arguments */
  for(j=0;j<fu->fi->len;j++) {

    if( fu->fi->arg[j]->should_NULL == TRUE) {
      continue; /* don't put it into the Perl prototype */
    }

    is_basic = is_basic_type_API(fu->fi->arg[j]->type);
    if( is_basic == TRUE) 
      fprintf(ofp,"\t%s %s\n",fu->fi->arg[j]->type,fu->fi->arg[j]->name);
    else
      fprintf(ofp,"\t%s%s %s\n",(is_basic == TRUE ? "" : package_name),fu->fi->arg[j]->type,fu->fi->arg[j]->name);
  }
  
  
  /* ok - now lets do the glue */
  if( strcmp(fu->fi->ret->type,"void") != 0 ) {
    if( fu->fi->ret->argtype == ARGTYPE_STATIC && is_basic_type_API(fu->fi->ret->type) == FALSE) {
      /* we need to put in the memory handler... */
      /*
       * Horrible hack. We assumme the the C type MyType * means dynamite type MyType.
       *
       * YUK!!!! Also - v.v.v.v.v bad hard coded Wise2
       *
       */
      fprintf(ofp,"\tINIT:\nWise2_%s temp;\n\tCODE:\n\ttemp = Wise2_hard_link_%s(%s%s(",fu->fi->ret->type,c2dyn_type(fu->fi->ret->type),package_name,fu->name);
    }
    else if(fu->fi->ret->argtype == ARGTYPE_STATIC &&  strcmp(fu->fi->ret->type,"char *") == 0 ) {
	fprintf(ofp,"\tINIT:\n\t%s temp;\n\tCODE:\n\ttemp = Wise2_stringalloc(%s%s(",fu->fi->ret->type,package_name,fu->name);
      }
    else { /* a basic type - don't need to handle, despite STATIC linkage */
      fprintf(ofp,"\tCODE:\n\tRETVAL = %s%s(",package_name,fu->name);
    } 
  } /* end of is not void */ 
  else { /* if it is void */
    fprintf(ofp,"\tCODE:\n\t%s%s(",package_name,fu->name);
  }
  
  for(j=0;j<fu->fi->len;j++) {
    if( fu->fi->arg[j]->argtype == ARGTYPE_OWNER && is_membasic_type_API(fu->fi->arg[j]->type) == FALSE && fu->fi->arg[j]->should_NULL == FALSE ) {
      if( (strcmp(fu->fi->arg[j]->type,"char*") == 0)  || (strcmp(fu->fi->arg[j]->type,"char *") == 0) ) {
	fprintf(ofp,"%sWise2_stringalloc(%s)",j == 0 ? "" : ",",fu->fi->arg[j]->name);
      } else {
	fprintf(ofp,"%sWise2_hard_link_%s(%s)",j == 0 ? "" : ",",c2dyn_type(fu->fi->arg[j]->type),fu->fi->arg[j]->name); 
      }
    } else {
      fprintf(ofp,"%s%s",j == 0 ? "" : ",",fu->fi->arg[j]->should_NULL == TRUE ? "NULL" : fu->fi->arg[j]->name);
    }
  }

  if( fu->fi->ret->argtype == ARGTYPE_STATIC && (is_basic_type_API(fu->fi->ret->type) == FALSE || strcmp(fu->fi->ret->type,"char *") == 0) ) {
    fprintf(ofp,"));\n\tRETVAL = temp;\n\tOUTPUT:\n\tRETVAL\n\n");
  } else {

    if( strcmp(fu->fi->ret->type,"void") != 0 ) {
      fprintf(ofp,");\n\tOUTPUT:\n\tRETVAL\n\n");
    }  else {
      fprintf(ofp,");\n\n");
    }
  }
  
  fprintf(ofp,"\n\n");

  return TRUE;
}
  

static char statbuf[128];

/* Function:  c2dyn_type(c)
 *
 * Descrip:    internal thing that maps C * types to dynamite types
 *
 *
 * Arg:        c [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 1029 "api.dy"
char * c2dyn_type(char * c)
{
  strcpy(statbuf,c);

  c = statbuf + strlen(statbuf)-1;
  for(;c != statbuf && *c != '*';c--)
    ;
  c--;
  for(;c != statbuf && isspace(*c);c--)
    ;
  c++;
  *c = '\0';
  return statbuf;
}
  
  

# line 1046 "api.dy"
boolean reconcile_dynAPI_with_FuncInfo(dynAPI * api,DYNFILE * dfp)
{
  int i;
  int j;
  boolean ret = TRUE;
  FuncInfo * fi;

  for(i=0;i<api->len;i++) {
    for(j=0;j<api->obj[i]->len;j++) {
      if( (fi= FuncInfo_from_name_DYNFILE(dfp,api->obj[i]->member[j]->name)) == NULL ) {
	warn("Cannot find any documentation for %s\n",api->obj[i]->member[j]->name);
	ret = FALSE;
      } else {
	api->obj[i]->member[j]->fi = fi;
      }
    }
    if(api->obj[i]->destructor != NULL) {
      if( (fi= FuncInfo_from_name_DYNFILE(dfp,api->obj[i]->destructor->name)) == NULL ) {
	warn("Cannot find any documentation for %s (destructor!)\n",api->obj[i]->destructor->name);
	ret = FALSE;
      } else {
	api->obj[i]->destructor->fi = fi;
      }
    }
  }

  for(j=0;j<api->non_len;j++) {
      if( (fi= FuncInfo_from_name_DYNFILE(dfp,api->non_obj[j]->name)) == NULL ) {
	warn("Cannot find any documentation for %s\n",api->non_obj[j]->name);
	ret = FALSE;
      } else {
	api->non_obj[j]->fi = fi;
      }
  }
  
  return ret;
}

# line 1084 "api.dy"
APIObject * APIObject_from_line(char * line,int maxline,FILE * ifp)
{
  APIObject * out;
  APIfunc * func;
  char * runner;

  if( strstartcmp(line,"object") != 0 ) {
    warn("In building APIObject, passed in a non object line [%s]",line);
    return NULL;
  }

  (void) strtok(line,spacestr);

  if( (runner=strtok(NULL,spacestr)) == NULL ) {
    warn("In building APIObject, no name for object",line);
    /* should skip to endobject ? */
    return NULL;
  }

  out = APIObject_alloc_std();
  out->name = stringalloc(runner);

  while( get_watched_line(line,maxline,ifp) != NULL ) {
    if( line[0] == '#' ) 
      continue;

    if( strstartcmp(line,"endobject") == 0 ) {
      break;
    }

    if( strstartcmp(line,"end") == 0 ) {
      warn("Got an end line [%s] but not an end object line. Don't like!");
      break;
    }

    if( strstartcmp(line,"func") == 0 ) {
      func = APIfunc_from_buffer(line);
      if( func == NULL )
	continue;
      add_APIObject(out,func);
    } else if ( strstartcmp(line,"des") == 0 ) {
      func = APIfunc_from_buffer(line);
      if( func == NULL )
	continue;
      out->destructor=  func;
    } else {
      warn("Did not understand [%s] as an APIObject line",line);
    }

  }


  return out;

}




# line 1143 "api.dy"
APIfunc * APIfunc_from_buffer(char * line)
{
  APIfunc * out;

  if( line == NULL ) {
    warn("Passed in a NULL line to APIfunc_from_buffer");
    return NULL;
  }
  if( strstartcmp(line,"func") != 0 && strstartcmp(line,"des") != 0 ) {
    warn("Passed a non function line to APIfunc_from_buffer");
    return NULL;
  }
  (void)strtok(line,spacestr);

  if( (line = strtok(NULL,spacestr)) == NULL ) {
    warn("line to APIfunc_from_buffer was empty. Yikes");
    return NULL;
  }

  out = APIfunc_alloc();
  if( out == NULL)
    return NULL;
  
  out->name = stringalloc(line);
  
  return out;
}



# line 1164 "api.c"
/* Function:  hard_link_APIfunc(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [APIfunc *]
 *
 * Return [UNKN ]  Undocumented return value [APIfunc *]
 *
 */
APIfunc * hard_link_APIfunc(APIfunc * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a APIfunc object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  APIfunc_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [APIfunc *]
 *
 */
APIfunc * APIfunc_alloc(void) 
{
    APIfunc * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(APIfunc *) ckalloc (sizeof(APIfunc))) == NULL)  {  
      warn("APIfunc_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->fi = NULL;  
    out->only_C = FALSE; 


    return out;  
}    


/* Function:  free_APIfunc(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [APIfunc *]
 *
 * Return [UNKN ]  Undocumented return value [APIfunc *]
 *
 */
APIfunc * free_APIfunc(APIfunc * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a APIfunc obj. Should be trappable");   
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
    if( obj->fi != NULL) 
      free_FuncInfo(obj->fi);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_APIObject(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_APIObject
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [APIfunc **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_APIObject(APIfunc ** list,int i,int j)  
{
    APIfunc * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_APIObject(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_APIObject which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [APIfunc **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_APIObject(APIfunc ** list,int left,int right,int (*comp)(APIfunc * ,APIfunc * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_APIObject(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_APIObject (list,++last,i);  
      }  
    swap_APIObject (list,left,last); 
    qsort_APIObject(list,left,last-1,comp);  
    qsort_APIObject(list,last+1,right,comp); 
}    


/* Function:  sort_APIObject(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_APIObject
 *
 *
 * Arg:         obj [UNKN ] Object containing list [APIObject *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_APIObject(APIObject * obj,int (*comp)(APIfunc *, APIfunc *)) 
{
    qsort_APIObject(obj->member,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_APIObject(obj,len)
 *
 * Descrip:    Really an internal function for add_APIObject
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [APIObject *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_APIObject(APIObject * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_APIObject called with no need");  
      return TRUE;   
      }  


    if( (obj->member = (APIfunc ** ) ckrealloc (obj->member,sizeof(APIfunc *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_APIObject, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_APIObject(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [APIObject *]
 * Arg:        add [OWNER] Object to add to the list [APIfunc *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_APIObject(APIObject * obj,APIfunc * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_APIObject(obj,obj->len + APIObjectLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->member[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_APIObject(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [APIObject *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_APIObject(APIObject * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->member[i] != NULL)    {  
        free_APIfunc(obj->member[i]);    
        obj->member[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  APIObject_alloc_std(void)
 *
 * Descrip:    Equivalent to APIObject_alloc_len(APIObjectLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [APIObject *]
 *
 */
APIObject * APIObject_alloc_std(void) 
{
    return APIObject_alloc_len(APIObjectLISTLENGTH); 
}    


/* Function:  APIObject_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [APIObject *]
 *
 */
APIObject * APIObject_alloc_len(int len) 
{
    APIObject * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = APIObject_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->member = (APIfunc ** ) ckcalloc (len,sizeof(APIfunc *))) == NULL)   {  
      warn("Warning, ckcalloc failed in APIObject_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_APIObject(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [APIObject *]
 *
 * Return [UNKN ]  Undocumented return value [APIObject *]
 *
 */
APIObject * hard_link_APIObject(APIObject * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a APIObject object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  APIObject_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [APIObject *]
 *
 */
APIObject * APIObject_alloc(void) 
{
    APIObject * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(APIObject *) ckalloc (sizeof(APIObject))) == NULL)  {  
      warn("APIObject_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->destructor = NULL;  
    out->hard_linker = NULL; 
    out->member = NULL;  
    out->len = out->maxlen = 0;  
    out->info = NULL;    
    out->sh = NULL;  


    return out;  
}    


/* Function:  free_APIObject(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [APIObject *]
 *
 * Return [UNKN ]  Undocumented return value [APIObject *]
 *
 */
APIObject * free_APIObject(APIObject * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a APIObject obj. Should be trappable"); 
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
    if( obj->destructor != NULL) 
      free_APIfunc(obj->destructor);     
    if( obj->hard_linker != NULL)    
      free_APIfunc(obj->hard_linker);    
    if( obj->member != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->member[i] != NULL)  
          free_APIfunc(obj->member[i]);  
        }  
      ckfree(obj->member);   
      }  
    if( obj->info != NULL)   
      free_ObjectInfo(obj->info);    
    if( obj->sh != NULL) 
      free_StructHolder(obj->sh);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_non_dynAPI(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_non_dynAPI
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [APIfunc **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_non_dynAPI(APIfunc ** list,int i,int j)  
{
    APIfunc * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_non_dynAPI(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_non_dynAPI which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [APIfunc **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_non_dynAPI(APIfunc ** list,int left,int right,int (*comp)(APIfunc * ,APIfunc * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_non_dynAPI(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_non_dynAPI (list,++last,i); 
      }  
    swap_non_dynAPI (list,left,last);    
    qsort_non_dynAPI(list,left,last-1,comp); 
    qsort_non_dynAPI(list,last+1,right,comp);    
}    


/* Function:  sort_non_dynAPI(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_non_dynAPI
 *
 *
 * Arg:         obj [UNKN ] Object containing list [dynAPI *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_non_dynAPI(dynAPI * obj,int (*comp)(APIfunc *, APIfunc *)) 
{
    qsort_non_dynAPI(obj->non_obj,0,obj->non_len-1,comp);    
    return;  
}    


/* Function:  expand_non_dynAPI(obj,len)
 *
 * Descrip:    Really an internal function for add_non_dynAPI
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [dynAPI *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_non_dynAPI(dynAPI * obj,int len) 
{


    if( obj->non_maxlen > obj->non_len )     {  
      warn("expand_dynAPInon_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->non_obj = (APIfunc ** ) ckrealloc (obj->non_obj,sizeof(APIfunc *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_dynAPI, returning FALSE");   
      return FALSE;  
      }  
    obj->non_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_non_dynAPI(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [dynAPI *]
 * Arg:        add [OWNER] Object to add to the list [APIfunc *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_non_dynAPI(dynAPI * obj,APIfunc * add) 
{
    if( obj->non_len >= obj->non_maxlen) {  
      if( expand_non_dynAPI(obj,obj->non_len + dynAPILISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->non_obj[obj->non_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_non_dynAPI(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [dynAPI *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_non_dynAPI(dynAPI * obj) 
{
    int i;   


    for(i=0;i<obj->non_len;i++)  { /*for i over list length*/ 
      if( obj->non_obj[i] != NULL)   {  
        free_APIfunc(obj->non_obj[i]);   
        obj->non_obj[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->non_len = 0;    
    return i;    
}    


/* Function:  swap_dynAPI(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_dynAPI
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [APIObject **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_dynAPI(APIObject ** list,int i,int j)  
{
    APIObject * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_dynAPI(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_dynAPI which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [APIObject **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_dynAPI(APIObject ** list,int left,int right,int (*comp)(APIObject * ,APIObject * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_dynAPI(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_dynAPI (list,++last,i); 
      }  
    swap_dynAPI (list,left,last);    
    qsort_dynAPI(list,left,last-1,comp); 
    qsort_dynAPI(list,last+1,right,comp);    
}    


/* Function:  sort_dynAPI(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_dynAPI
 *
 *
 * Arg:         obj [UNKN ] Object containing list [dynAPI *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_dynAPI(dynAPI * obj,int (*comp)(APIObject *, APIObject *)) 
{
    qsort_dynAPI(obj->obj,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_dynAPI(obj,len)
 *
 * Descrip:    Really an internal function for add_dynAPI
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [dynAPI *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_dynAPI(dynAPI * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_dynAPI called with no need"); 
      return TRUE;   
      }  


    if( (obj->obj = (APIObject ** ) ckrealloc (obj->obj,sizeof(APIObject *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_dynAPI, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_dynAPI(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [dynAPI *]
 * Arg:        add [OWNER] Object to add to the list [APIObject *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_dynAPI(dynAPI * obj,APIObject * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_dynAPI(obj,obj->len + dynAPILISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->obj[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_dynAPI(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [dynAPI *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_dynAPI(dynAPI * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->obj[i] != NULL)   {  
        free_APIObject(obj->obj[i]); 
        obj->obj[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  dynAPI_alloc_std(void)
 *
 * Descrip:    Equivalent to dynAPI_alloc_len(dynAPILISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [dynAPI *]
 *
 */
dynAPI * dynAPI_alloc_std(void) 
{
    return dynAPI_alloc_len(dynAPILISTLENGTH);   
}    


/* Function:  dynAPI_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [dynAPI *]
 *
 */
dynAPI * dynAPI_alloc_len(int len) 
{
    dynAPI * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = dynAPI_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->non_obj = (APIfunc ** ) ckcalloc (len,sizeof(APIfunc *))) == NULL)  {  
      warn("Warning, ckcalloc failed in dynAPI_alloc_len");  
      return NULL;   
      }  
    out->non_len = 0;    
    out->non_maxlen = len;   


    if((out->obj = (APIObject ** ) ckcalloc (len,sizeof(APIObject *))) == NULL)  {  
      warn("Warning, ckcalloc failed in dynAPI_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_dynAPI(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [dynAPI *]
 *
 * Return [UNKN ]  Undocumented return value [dynAPI *]
 *
 */
dynAPI * hard_link_dynAPI(dynAPI * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a dynAPI object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  dynAPI_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [dynAPI *]
 *
 */
dynAPI * dynAPI_alloc(void) 
{
    dynAPI * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(dynAPI *) ckalloc (sizeof(dynAPI))) == NULL)    {  
      warn("dynAPI_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->non_obj = NULL; 
    out->non_len = out->non_maxlen = 0;  
    out->obj = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_dynAPI(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [dynAPI *]
 *
 * Return [UNKN ]  Undocumented return value [dynAPI *]
 *
 */
dynAPI * free_dynAPI(dynAPI * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a dynAPI obj. Should be trappable");    
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
    if( obj->non_obj != NULL)    {  
      for(i=0;i<obj->non_len;i++)    {  
        if( obj->non_obj[i] != NULL) 
          free_APIfunc(obj->non_obj[i]); 
        }  
      ckfree(obj->non_obj);  
      }  
    if( obj->obj != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->obj[i] != NULL) 
          free_APIObject(obj->obj[i]);   
        }  
      ckfree(obj->obj);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif

/*  Last edited: Apr 29 10:16 1997 (birney) */



#include "module.h"
#include "dyna2.h"
#include "telegraph.h"
#include "dynafunc.h"
#include "modulefunc.h"
#include "api.h"

/* remove this #define if you don't want to link to compugen conversion code */
#define COMPUGEN

#ifdef COMPUGEN
#include "compugen.h"
#endif

char * program_name = "Dynamite Compiler";
char * description  = "This takes a Dynamite source file to make .c and .h files";
char * version      = "1.1";
char * author       = "Ewan Birney";
char * email        = "birney@sanger.ac.uk";

FILE * tele_file =NULL;

/*boolean do_compugen = FALSE;*/

typedef int FailureType;

#define   FailureType_dyc_All 01
#define   FailureType_dyc_UF  02

DPImplementation * dpi;


void show_usage(FILE * ofp)
{
  fprintf(ofp,"%s %s\nWritten by %s <%s>\nUSAGE: dyc [options] <dynamite source files>\n",program_name,version,author,email);
  fprintf(ofp,"OPTIONS:\n");
  fprintf(ofp,"  -u -h Show this message \n\n");
/*
  fprintf(ofp,"  -i  produce information on STDOUT of Dynamite types\n");
  fprintf(ofp,"  -m  file containing method information\n"); 
  fprintf(ofp,"  -U  Do not load configured method information, only -m options\n");
*/
  fprintf(ofp,"  -I  <configdir> Overrides WISECONFIGDIR\n"); /* arve */
  fprintf(ofp,"  -n  package name for file\n");
  fprintf(ofp,"  -P  [no] Protection level of the code (run-time warning checks)\n\n  Advanced dy compiler options for Dynamite development\n");
  fprintf(ofp,"  -l  hard link option           \n");
  fprintf(ofp,"  -F  Fail on incomplete modules \n");
  fprintf(ofp,"  -D  Warn on undocumented funcs \n");
  fprintf(ofp,"API layer options\n");
  fprintf(ofp,"  -perl Generate perl xs and stubs from api\n");
  fprintf(ofp,"  -latex Generate latex section\n");
  fprintf(ofp,"  -tele <filename> generate telegraph pseudo-grammar\n");
  show_help_DPImplementation(ofp);

}


/***
  Error system: to be registered
  in wiseerror.c
  
  Will give line and position of
  error'd file.
  
  ***/

char * parsed_file = NULL;
static char pbuff[1024];

char * parse_line_error(void)
{
  sprintf(pbuff,"%s:%d:",parsed_file,get_watched_linecount()); /* arve */
  return pbuff;
}




/*** 
  parse function parses the function section of the file,
  returning NULL on error. input is the pointer to the .dy
  stream, dfp is the DYNFILE (found in dynfile.dy) output structure.

  ModuleFunction * is going to have a list of all the
  Modules which have constructor or deconstructor functions
  like  Module_alloc()
  and   free_Module()

  if they are in place, then we will write out different
  cons/decons in the function files 

****/




ModuleFunctionList * parse_function_section(FILE * input,DYNFILE * dfp,FailureType * ft,boolean addnumbers)
{
  ModuleFunctionList * out;
  boolean outoffunction=TRUE;
  char * name;
  boolean next_line_alloc = FALSE;
  boolean next_line_free  = FALSE;
  char buffer[MAXLINE];
  char * tempb;
  ModuleFunction * mf;
  FuncInfo * fi = NULL;


  out = ModuleFunctionList_alloc_std();

  while( get_watched_line(buffer,MAXLINE,input) != NULL) {
	chop_newline(buffer);

      if( strstartcmp(buffer,"%}") == 0)
	break;

      if( strstartcmp(buffer,"%func") == 0 ) {
    /*	fprintf(stderr,"Getting info line!\n"); */
	fi =  read_FuncInfo_line(buffer,input);
	fi->is_hand_written = TRUE;
	if( fi == NULL ) {
		warn("Unable to read function info line at %s. Probably going to give a parse error",buffer);
		}
	else add_DYNFILE(dfp,fi);
	continue;
      }
	


      if( outoffunction == TRUE ) {
	
	if( buffer[0] == '{' ) {	    
	  outoffunction = FALSE;
	  fputs_func_DYNFILE(buffer,dfp);
	  continue;
	}
	
	/*	fprintf(stderr,"Got here [%s]\n",buffer); */
	
	if( buffer[0] == '\0' || isspace(buffer[0]) || strstartcmp(buffer,"static") == 0 || buffer[0] == '#' 
	    || buffer[0] == '/' || buffer[0] == '*' ) {
	  next_line_alloc = next_line_free = FALSE;
	  fputs_func_DYNFILE(buffer,dfp);
	  continue;
	}



	/*** if the line starts with a ! we want to interpret it ***/
	/*** as a dynamite specific tag ***/
	
	
	if( strstartcmp(buffer,"!constructor") == 0 ) {
	  next_line_alloc = TRUE;
	  continue;
	}

	if( strstartcmp(buffer,"!deconstructor") == 0 ) { 
	  next_line_free = TRUE;
	  continue;
	}
	  
	  
	if ( buffer[0] == '!' ) {
	  next_line_alloc = next_line_free = FALSE;
	  continue;
	}


	/*fprintf(stderr,"Got here again [%s]\n",buffer);*/

	/*** ok this should be a function! ****/
	  
	
	if( next_line_alloc == TRUE || next_line_free == TRUE ) {
	  tempb = stringalloc(buffer);
	  name = parse_and_get_module_name_from_func(tempb,next_line_alloc);
	  if(name == NULL) {
	    warn("Function [%s] specified as cons/decons, but no module name parsed...",buffer);
	  }
	  else {
	    if( (mf=get_ModuleFunction_from_name(out,name)) == NULL ) 
	      mf = new_ModuleFunction(out,name);
	    
	    if( next_line_alloc == TRUE ) {
	      if( mf->has_cons == TRUE ) 
		warn("Module %s already has a constructor... double overload!",name);
	      else mf->has_cons = TRUE;
	    }
	    else {
	      if( mf->has_decons == TRUE ) 
		warn("Module %s already has a deconstructor... double overload",name);
	      else mf->has_decons = TRUE;
	    }
	    ckfree(name);
	  }
	  ckfree(tempb);
	} /*** end of if cons/decons ***/
	
	
	next_line_alloc = next_line_free = FALSE;
	
	
	/*** reconcile with fi if any ***/
	
	if( fi != NULL ) {
	  reconcile_FuncInfo_with_funcstr(fi,buffer);
	  dfp->funcpos += show_eddystyle_FuncInfo(fi,dfp->func);
	  /** already added to DYNFILE **/
	} else { 
	  /*  fprintf(stderr,"Using [%s] as an unknown function\n",buffer); */
	  fi = unknown_user_FuncInfo(buffer);
	  /*  fprintf(stderr,"Unknown function %s...\n",fi->name);*/
	  add_DYNFILE(dfp,fi);
	}


	/*** put #line compiler information ***/
	
	/*** 
	  using global parsed file which 
	  is also used for errors... sort of yukky.
	  ***/
	if( addnumbers == TRUE) {
	  fprintf(dfp->func,"# line %d \"%s\"\n",get_watched_linecount(),parsed_file);
	  dfp->funcpos++;
	}
	
	/*** put line ***/
	
	fi->line_in_c = dfp->funcpos;
	fi = NULL;
	fputs_func_DYNFILE(buffer,dfp);
	
	
	/*** put away function now ***/
	
	
	
	/*** header information will be put away in dfp call ***/
	
      } else {
	
	/*** actually inside a function ***/
	
	if( buffer[0] == '}' )
	  outoffunction = TRUE;
	fputs_func_DYNFILE(buffer,dfp);
      }
  }

  return out;
}


/***
  basic top-level loop over the file sourcename.

  mts = method scope file

  write_html_now does not currently work!

  ***/

typedef struct {
  char * c_extension_name;
  char * t_extension_name;
  char * pfdoc_ext;
  char * xs_ext;
  char * typemap_ext;
  char * pod_ext;
  char * latex_ext;
  boolean make_perl;
  boolean make_latex;
  boolean all_callable;
} APIpara;

boolean should_make_API(APIpara * api)
{
  if( api->c_extension_name != NULL ) {
    if( api->t_extension_name == NULL ) {
      warn("You can't specifiy a c extension name without a type extension name as well!");
      return FALSE;
    }
    
    return TRUE;
  } 

  if( api->make_perl == TRUE ) {
    return TRUE;
  }

  return FALSE;
}

boolean make_API(DYNFILE * dfp,Dynamite * dyn,APIpara * api_para)
{
  int i;
  dynAPI * api = NULL;
  char buffer[512];
  char strpack[64];
  FILE * ofp;

  strcpy(strpack,dfp->package_name);
  if( strpack[strlen(strpack)-1] == '_') {
    strpack[strlen(strpack)-1] = '\0';
  }

  if( api_para->all_callable == TRUE) {
    api = promote_all_callable_dynAPI(dyn,dfp);
  } else {
    for(i=0;i<dyn->len;i++) {
      if( dyn->list[i]->type == DYNAMITETYPEAPI ) {
	if( api == NULL ) 
	  api = (dynAPI*)dyn->list[i]->data;
	else {
	  warn("You have got more than one API defintion for this module. Currently not handled\n");
	  return FALSE;
	}
	
      }
    }
    
    
    if( api == NULL ) {
      return TRUE;
    }

    hard_link_dynAPI(api);
    crosslink_to_StructHolder_dynAPI(api,dyn);

  }

  if( api_para->make_perl == TRUE ) {
    write_dynAPI_accessor_functions(dfp,api);
  }


  if( prepare_dynAPI(api,dfp,dyn) == FALSE ) {
    warn("Could not prepare API");
    return FALSE;
  }

  sprintf(buffer,"%s%s",dfp->sourceroot,api_para->c_extension_name);

  ofp = openfile(buffer,"w");
  if( ofp == NULL ){
    warn("Could not open [%s] as a API file",buffer);
    return FALSE;
  }

  write_C_dynAPI(api,dfp->package_name,ofp);

  fclose(ofp);

  sprintf(buffer,"%s%s",dfp->sourceroot,api_para->t_extension_name);

  ofp = openfile(buffer,"w");
  if( ofp == NULL ){
    warn("Could not open [%s] as a API file",buffer);
    return FALSE;
  }

  write_type_C_dynAPI(api,dfp->package_name,ofp);

  fclose(ofp);


  if( api_para->make_latex == TRUE ) {
    if( api_para->latex_ext == NULL ) {
      api_para->latex_ext = "tex";
    }
    sprintf(buffer,"%s.%s",dfp->sourceroot,api_para->latex_ext);

    ofp = openfile(buffer,"w");
    if( ofp == NULL ){
      warn("Could not open [%s] as a API file",buffer);
      return FALSE;
    }
    write_latex_dynAPI(api,dfp->sourceroot,strpack,ofp);
    fclose(ofp);
  }

  if( api_para->pfdoc_ext != NULL ) {
    sprintf(buffer,"%s%s",dfp->sourceroot,api_para->pfdoc_ext);

    ofp = openfile(buffer,"w");
    if( ofp == NULL ){
      warn("Could not open [%s] as a API file",buffer);
      return FALSE;
    }

    fprintf(ofp,"\n:Module %s\n\n",dfp->sourceroot);
    write_pfdoc_dynAPI(api,dfp->package_name,ofp);
    fprintf(ofp,"\n\n!Module\n\n");
    
    fclose(ofp);
  }
  if( api_para->make_perl == TRUE ) {
    sprintf(buffer,"%s%s",dfp->sourceroot,api_para->pod_ext);
    ofp = openfile(buffer,"w");
    if( ofp == NULL ){
      warn("Could not open [%s] as a API file",buffer);
      return FALSE;
    }
    write_pod_dynAPI(api,dfp->sourceroot,strpack,ofp);
    fclose(ofp);

    sprintf(buffer,"%s%s",dfp->sourceroot,api_para->xs_ext);

    ofp = openfile(buffer,"w");
    if( ofp == NULL ){
      warn("Could not open [%s] as a API file",buffer);
      return FALSE;
    }


    write_XS_dynAPI(api,dfp->package_name,ofp);
    fclose(ofp);

    sprintf(buffer,"%s%s",dfp->sourceroot,api_para->typemap_ext);

    ofp = openfile(buffer,"w");
    if( ofp == NULL ){
      warn("Could not open [%s] as a API file",buffer);
      return FALSE;
    }


    write_XS_typemap_dynAPI(api,dfp->package_name,ofp);
    fclose(ofp);
  }


  /*  free_dynAPI(api); */
    
  return TRUE;
}
  

boolean do_a_file(char * sourcename,MethodTypeSet * mts,boolean write_html_now,int c_protection_level,boolean should_hard_link,boolean should_warn_undocs,boolean addnumbers,char * package,APIpara * api_para,FailureType * ft)
{
  FILE * ifp = NULL;  
  DYNFILE * dfp = NULL;
  char * outputname = NULL;
  char * runner;
  char * rt;
  char c;
  char buffer[MAXLINE];
  Dynamite * dyn;
  ModuleFunctionList * mfl = NULL;
  boolean end_of_file;
  ParseError pfail;
  int no;
  int i;
#ifdef COMPUGEN
  OneModel * om;
#endif



  /* Read first block */

  /*	read_plain_scan_and_replace("test.kar"); */

  /*** open source ***/

  ifp = openfile(sourcename,"r");
  if( ifp == NULL  ){
    warn("Could not open [%s] as Dynamite source",sourcename);
    return FALSE;
  }
  open_watch_file(ifp);
  parsed_file = stringalloc(sourcename);
  push_errormsg_stack_call(parse_line_error);
  

  /*** open DYNFILE ***/
  runner = strrchr(sourcename,'.');
  if( runner == NULL ) {
    warn("Source %s has no . character. Not critical, but worrying...",sourcename);
    runner = sourcename + strlen(sourcename);
  }

  /*** hard coding debug for the moment ***/


  c = *runner;
  *runner = '\0';

  if( (dfp = open_std_no_html_dynfile(sourcename)) == NULL) {
    warn("Could not open DYNFILE system for output %s\n",outputname);
    goto error;
  }

  dfp->code_debug_level = c_protection_level;
  dfp->should_hard_link = should_hard_link;
  if( package != NULL ) {
    dfp->package_name = stringalloc(package);
  }

  *runner = c;


  /* Skip to %{ */

  while( (rt=get_watched_line(buffer,MAXLINE,ifp)) != NULL) {
      if(strstartcmp(buffer,"%{") == 0)
	break;
    }

  if( rt == NULL) {
    log_full_error(WARNING,0,"Exhausted dynamite source [%s] without reading any blocks - consult man/help pages!",sourcename);
    goto error;
  }
	
  /* 
   * Header region
   * Dump until next block %} 
   *
   */

  while( (rt=get_watched_line(buffer,MAXLINE,ifp)) != NULL) {
      if(strstartcmp(buffer,"%}") == 0)
	break;
      fputs(buffer,dfp->head);
    }


  if( rt == NULL) {
    warn("had only header information in dynamite file, no close header %}");
    goto error;
  }

  /* Load in Dynamite */


  dyn = read_Dynamite(ifp,"%{",&end_of_file,mts);
	
  if( dyn == NULL ) {
      log_full_error(WARNING,0,"Unable to read Dynamite models. Parse error.");
      goto error;
    }




  /* probably should check it now! */

  pop_errormsg_stack();
	
  pfail = (PERR_ARG_MISTYPE | PERR_ARG_NUM_MIS | PERR_SYNTAX);

  if( prepare_Dynamite(dyn,mts,dpi->dycw,pfail) == FALSE ) {
      log_full_error(FATAL,0,"A Dynamite blueprint fails semantic checks. Please refer to previous errors for the precise problem");
      goto error;
    }

  if( tele_file != NULL ) {
    for(i=0;i<dyn->len;i++) {
      if( dyn->list[i]->type == DYNAMITETYPEGENERICMATRIX ) {
	write_Telegraph_grammar((GenericMatrix*)dyn->list[i]->data,tele_file);
      }
    }
  }


  /*** writes out user functions ***/

  
  if( end_of_file == FALSE) {
    push_errormsg_stack_call(parse_line_error);
    mfl = parse_function_section(ifp,dfp,ft,addnumbers);
    pop_errormsg_stack();
  }
  else {
    write_Dynamite_minimal_func(dfp);
    mfl = NULL;
  }


  /*** write out the header ***/

  write_Dynamite_header(dyn,dpi,dfp);

  /*** writes out dynamite functions ***/

  /*** should write line number now ***/

  fprintf(dfp->func,"# line %d \"%s.c\"\n",dfp->funcpos,dfp->sourceroot);
  dfp->funcpos++;
  write_Dynamite_func(dyn,mfl,dpi,dfp);

#ifdef COMPUGEN

  if( dyn->len > 0 && dyn->list[0]->type == DYNAMITETYPEGENERICMATRIX && dpi->doone == TRUE) {
      
      om = OneModel_from_GenericMatrix((GenericMatrix *) dyn->list[0]->data);
      if( om == NULL ) {
	warn("Cannot built Compugen port!");
      } else {
	write_cugen_funcs(om,mts,dfp);
      }
  }
#endif

  positionify_DYNFILE(dfp);


  /*
     check errors... could fail here

     */

  if( have_got_unknown_func_DYNFILE(dfp,&no,should_warn_undocs) == TRUE) {
    *ft = (*ft | FailureType_dyc_UF);
    warn("You have %d undocumented functions",no);
  }


  if( should_make_API(api_para) == TRUE ) 
    make_API(dfp,dyn,api_para);

  place_declarations_DYNFILE(dfp,FO_CUI_POSITION);

  /***
    close Dynfile

    ***/

  if(dyn != NULL)
    free_Dynamite(dyn);
  if( mfl != NULL ) 
    free_ModuleFunctionList(mfl);

  close_dynfile(dfp);
  close_watch_file();
  fclose(ifp);
  return TRUE;

  error :
  if( ifp != NULL) {
    fclose(ifp);
    close_watch_file();
  }
  if( dfp != NULL) 
    close_dynfile(dfp);
  return FALSE;  
	
}



int main(int argc,char * argv[]) 
{
  FailureType fail = 0 ;
  FailureType should_fail_on = 0;
  int i;
  boolean doinfo = FALSE;
  boolean noaddnumbers = FALSE;
  MethodTypeSet * mts;
  MethodTypeSet * cp;
  boolean no_config_mts = FALSE;
  int prot_level = 0;
  int should_hard_link = 0;
  boolean should_warn_undoc = FALSE;
  char * prot_str;
  char * runner;
  char *config_dir=NULL;
  char buffer[64]; /** really for removing files **/
  char * telegraph;
  
  APIpara api;

  char * pack;
	
  /** we no longer read in configs **/

  mts = standard_dynamite_MethodTypeSet();


  if( strip_out_boolean_argument(&argc,argv,"h") == TRUE 
      || strip_out_boolean_argument(&argc,argv,"u") == TRUE /* arve */
      || argc == 1 ) {
	show_usage(stdout);
	exit(1);
	}

  noaddnumbers   = strip_out_boolean_argument(&argc,argv,"m");
  doinfo        = strip_out_boolean_argument(&argc,argv,"i");
  no_config_mts = strip_out_boolean_argument(&argc,argv,"U");
  should_hard_link = strip_out_boolean_argument(&argc,argv,"l");
  prot_str      = strip_out_assigned_argument(&argc,argv,"P");
  should_warn_undoc = strip_out_boolean_argument(&argc,argv,"D");
  telegraph    = strip_out_assigned_argument(&argc,argv,"tele");

  pack = strip_out_assigned_argument(&argc,argv,"n");

  api.xs_ext = NULL;
  api.typemap_ext = NULL;
  api.pod_ext = NULL;

  api.c_extension_name = strip_out_assigned_argument(&argc,argv,"a");
  api.t_extension_name = strip_out_assigned_argument(&argc,argv,"b");
  api.pfdoc_ext = strip_out_assigned_argument(&argc,argv,"p");
  api.xs_ext = strip_out_assigned_argument(&argc,argv,"x");
  api.typemap_ext = strip_out_assigned_argument(&argc,argv,"tym");
  api.all_callable = strip_out_boolean_argument(&argc,argv,"c");
  api.make_perl = strip_out_boolean_argument(&argc,argv,"perl");
  api.latex_ext = strip_out_assigned_argument(&argc,argv,"exttex");
  api.make_latex = strip_out_boolean_argument(&argc,argv,"latex");


  if( api.make_perl == TRUE) {
   
    if( api.xs_ext == NULL ) {
      api.xs_ext = ".xs";
    }
    if( api.typemap_ext == NULL ) {
      api.typemap_ext = ".typemap";
    }
    if( api.pod_ext == NULL ) {
      api.pod_ext = ".pod";
    }
  }


  if( strip_out_boolean_argument(&argc,argv,"F") == TRUE) {
    should_fail_on = FailureType_dyc_All;
  }

  /* do DPImplementation */

  dpi = new_DPImplementation_from_argstr(&argc,argv);

  if( prot_str != NULL ) {
      if( is_integer_string(prot_str,&prot_level) == FALSE ) {
	  warn("Protection level %s is no integer!");
	  prot_level = 0;
      }
  }

				/* Override/set WISECONFIGDIR on the cmdline. (arve) */
  config_dir = strip_out_assigned_argument(&argc, argv, "I");
  if (config_dir != NULL) {
      set_config_dir(config_dir); 
  }

  if( read_into_MethodTypeSet_filename(mts,"methods") == FALSE){
    warn("You have no config file called 'methods'. This is bad news for dynamite matrices. I will attempt to compile, but you cannot use logical types. 'methods' should be either in the current directory, the $WISECONFIGDIR or your $WISEPERSONALDIR");
  }



  /*** ok,loop over and do it ***/

  if( argc < 1 ) {
    warn("You must have at least one dynamite source file to compile!");
    show_usage(stdout);
    exit(1);
  }

  if( telegraph != NULL ) {
    tele_file= fopen(telegraph,"w");
  }

  
  for(i=1;i<argc;i++) {
    if( mts != NULL) 
      cp = copy_MethodTypeSet(mts); /* actually very cheap */
    if( do_a_file(argv[i],mts,FALSE,prot_level,should_hard_link,should_warn_undoc,noaddnumbers == TRUE ? FALSE : TRUE,pack,&api,&fail) == FALSE ) {
      
      fatal("Terminated dyc one %d argument %s",i,argv[i]);
    }
    if( (should_fail_on == 01 && fail != 0) || (fail & should_fail_on) ) {

      /*** remove files which fail ****/

      /*** ugh this should be done better ***/

      for(runner=argv[i]+strlen(argv[i]) - 1;runner > argv[i] && *runner != '.';runner--)
	;
      if( runner != argv[i] ) {
	*runner = '\0';
	sprintf(buffer,"%s.c",argv[i]);
	if( remove_file(buffer) == FALSE ) {
	  warn("Could not remove file %s from filesystem",buffer);
	}
	sprintf(buffer,"%s.h",argv[i]);
	if( remove_file(buffer) == FALSE ) {
	  warn("Could not remove file %s from filesystem",buffer);
	}
      }
      /* else - well - something bad has happened */
      
      fatal("Failed on file %s due to user defined fails",argv[i]);
    }

    if( mts != NULL ) {
      free_MethodTypeSet(mts);
      mts = cp;
    }
  }

  free_MethodTypeSet(mts);

  return 0;
}  



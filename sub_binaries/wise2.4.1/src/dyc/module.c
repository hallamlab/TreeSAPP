/*  Last edited: Apr 29 10:12 1997 (birney) */
 
/********************************************************/
/* Wisetools version 3 file                             */
/*                                                      */
/* this file is copyright (c) Ewan Birney 1996          */
/* this file is part of the wisetools sequence analysis */
/* package                                              */
/********************************************************/

/********************************************************/
/* Dynamite structure - a wrapper for the variety of    */
/* of different Dynamite inputs possible                */
/*                                                      */
/*                                                      */
/********************************************************/

/***** RCS Info *****************************************/
/*
   $Id: module.c,v 1.1.1.1 2001/06/18 13:59:56 birney Exp $

   $Log: module.c,v $
   Revision 1.1.1.1  2001/06/18 13:59:56  birney
   moved wise2 to ensembl cvs repository

   Revision 1.5  1999/05/03 14:52:58  birney
   added in alot of fixes to dyanmite.

   added in def boolean's in commandline

   added in seqalign/pwmdna stuff and linked it to hmmer bridge.

   Revision 1.4  1998/12/15 10:04:25  birney
   parsing '#' as comments now in dynamite files

 * Revision 1.3  1998/12/03  23:01:28  birney
 * getting the latex from teh api docs
 *
 * Revision 1.2  1998/11/16  15:24:59  birney
 * started the threads .... ;)
 *
 * Revision 1.1.1.1  1998/08/28  09:30:57  birney
 * Wise2
 *
 * Revision 1.8  1998/07/27  11:39:36  birney
 * before compugen merge
 *
 * Revision 1.7  1998/01/09  15:12:45  birney
 * Added api handling, namespace protection and pdoc production. Lots of bugs
 * out from the main matrix construction including start point errors on
 * packaln
 *
 * Revision 1.6  1997/07/24  15:40:48  birney
 * added MethodTypeSet support for on-the-fly methods
 *
 * Revision 1.5  1997/06/10  13:55:37  birney
 * Lots of things: moved to teh dynfile system of keeping files around, addedfunction information, started to remodelled the matrix part
 *
 * Revision 1.4  1996/12/06  11:17:23  birney
 * added prepare display
 *
 * Revision 1.3  1996/10/14  10:38:37  birney
 * added display handling
 *
 * Revision 1.2  1996/10/04  16:41:42  birney
 * added cons/decons control for basic structures
 *
 * Revision 1.1  1996/03/04  12:15:33  birney
 * Initial revision
 *

*/
/********************************************************/

#include "module.h"

#define WRITECHEADER

#include "dyna2.h"
#include "wisec.h"


boolean prepare_dynAPI(dynAPI * api,DYNFILE * dfp,Dynamite * dyn)
{
  int i,j;
  boolean ret = TRUE;
  StructHolder * sh;

  for(i=0;i<api->len;i++) {
    if(api->obj[i]->destructor == NULL ){
      warn("Dynamite object %s does not have a deconstructor. Yikes!",api->obj[i]->name);
      ret = FALSE;
    }
  }

  if ( reconcile_dynAPI_with_FuncInfo(api,dfp) == FALSE)
    ret = FALSE;

  for(i=0;i<api->len;i++) {
    for(j=0;j<dyn->len;j++) {
      if( dyn->list[j]->type == DYNAMITETYPESTRUCT ) {
	sh = (StructHolder *)dyn->list[j]->data;
	if( strcmp(sh->name,api->obj[i]->name) == 0 ) {
	  if( sh->oi != NULL )
	    api->obj[i]->info = hard_link_ObjectInfo(sh->oi);
	  break;
	}
      } 
    }
    if( j == dyn->len ) {
      warn("Dynamite api object %s has no actual physical object!\n",api->obj[i]->name);
      ret = FALSE;
    }
  }
    

  return ret;
}

boolean crosslink_to_StructHolder_dynAPI(dynAPI * api,Dynamite * dyn)
{
  int i,j;
  StructHolder * sh;

  for(j=0;j<api->len;j++) {

    for(i=0;i<dyn->len;i++) {
      if( dyn->list[i]->type != DYNAMITETYPESTRUCT ) 
	continue;

      sh = (StructHolder *)dyn->list[i]->data;

      if( strcmp(api->obj[j]->name,sh->name) == 0 ) {
	break;
      }
    }
	
    if( i == dyn->len ) {
      warn("Api object %s is not supported by a dynamite object. Weird!",api->obj[j]->name);
      continue;
    }

    api->obj[j]->sh = hard_link_StructHolder(sh);
  }

  return TRUE; /* should error check! Nggggg!*/
}  


dynAPI * promote_all_callable_dynAPI(Dynamite * dyn,DYNFILE * dfp)
{
  int i,j;
  StructHolder * sh;
  dynAPI * out;
  char buffer[512];
  char * first;
  APIObject * obj;
  APIfunc * func;
  
  out = dynAPI_alloc_std();
    
  /*** ok  - promote each struct to an object ***/


  for(i=0;i<dyn->len;i++)
    if( dyn->list[i]->type == DYNAMITETYPESTRUCT ) {
      sh = (StructHolder *)dyn->list[i]->data;
      obj = APIObject_alloc_std();
      obj->name = stringalloc(sh->name);
      sprintf(buffer,"free_%s",sh->name);
      obj->destructor = APIfunc_alloc();
      obj->destructor->name = stringalloc(buffer);
      obj->sh = hard_link_StructHolder(sh);
      add_dynAPI(out,obj);
    }
  
  for(i=0;i<dfp->len;i++) {
    if( dfp->info[i]->functype != FI_CALLABLE ) {
      continue;
    }

    func = APIfunc_alloc();
    func->name = stringalloc(dfp->info[i]->name);
    if( dfp->info[i]->len > 0 ) {
      
      first = dfp->info[i]->arg[0]->type;
      for(j=0;j<out->len;j++) {
	if( strstartcmp(first,out->obj[j]->name) == 0 ) {
	  add_APIObject(out->obj[j],func);
	  break;
	}
      } 
      if( j >= out->len ) {
	/* no object */
	add_non_dynAPI(out,func);
      }
    } else {
      add_non_dynAPI(out,func);
    }
  }

  return out;
  
}




void write_Dynamite_func(Dynamite * dyn,ModuleFunctionList * mfl,DPImplementation * dpi,DYNFILE * dfp)
{
  register int i;




  for(i=0;i<dyn->len;i++) {
    switch ( dyn->list[i]->type ) {
    case DYNAMITETYPEGENERICMATRIX :
      write_GenericMatrix_func(dfp,(GenericMatrix *)dyn->list[i]->data,dpi);
      break;
    case DYNAMITETYPESTRUCT     :
      write_StructHolder_function(dfp,(StructHolder *)dyn->list[i]->data,mfl);
      break;
    case DYNAMITETYPEDISPLAY :
      write_Aln2Display_function(dfp,(Aln2Display *) dyn->list[i]->data);
      break;
    case DYNAMITETYPEFRIEND :
      break;
    case DYNAMITETYPEAPI :
      break; 
    default	: log_full_error(WARNING,0,"Cannot not currently write dynamite type %d\n",dyn->list[i]->type);
    }
  }

}

void write_Dynamite_header(Dynamite * dyn,DPImplementation * dpi,DYNFILE * dfp)
{
  register int i;


  for(i=0;i<dyn->len;i++) {
    switch ( dyn->list[i]->type ) {
    case DYNAMITETYPEGENERICMATRIX :
      write_GenericMatrix_header  (dfp,(GenericMatrix *)dyn->list[i]->data,dpi);
      break;
    case DYNAMITETYPESTRUCT     :
      write_StructHolder_header  (dfp,(StructHolder *)dyn->list[i]->data);
      break;
    case DYNAMITETYPEDISPLAY : 
      write_Aln2Display_header (dfp,(Aln2Display *) dyn->list[i]->data);
      break;
    case DYNAMITETYPEFRONTEND   :
      break;
    case DYNAMITETYPEFRIEND :
      write_Friend_header(dfp,(Friend *) dyn->list[i]->data);
      break;
    case DYNAMITETYPEAPI :
      break; /*** API's produce *different* header files ***/
    default	: log_full_error(WARNING,0,"Cannot not currently write dynamite type %d\n",dyn->list[i]->type);
    }
  }
  
}

boolean  prepare_Dynamite(Dynamite * dyn,MethodTypeSet * mts,DycWarning * dycw,ParseError fail)
{
  register int i;
  boolean ret = TRUE;
  boolean temp;

  for(i=0;i<dyn->len;i++) {
    switch( dyn->list[i]->type ) {
    case DYNAMITETYPEGENERICMATRIX :
      temp = prepare_matrix( (GenericMatrix *) dyn->list[i]->data,mts,dycw,fail);
      if( temp == FALSE )
	ret = FALSE;
      break;
    case DYNAMITETYPEDISPLAY :
      if( prepare_Aln2Display( (Aln2Display *) dyn->list[i]->data) == FALSE)
	ret = FALSE;
      break;
    case DYNAMITETYPEAPI :
      break; /* maybe something in here */
    default :
      break;
    }
  }			
  
  return ret;
}

Dynamite * read_Dynamite(FILE * ifp,char * fend,boolean * end_of_file,MethodTypeSet * mts)
{
  char buffer[MAXLINE];
  char * runner;
  GenericMatrix * gm;
  StructHolder * sh;
  Aln2Display  * dis;
  Friend * fr;
  dynAPI * api;
  Dynamite * out;
  boolean ret = TRUE;

  *end_of_file = FALSE;

  out = Dynamite_alloc_std();

  while( get_watched_line(buffer,MAXLINE,ifp) != NULL) {
    if( strstr(buffer,fend) != NULL)
      break;

    chop_newline(buffer);

    if( buffer[0] == '#' ) {
      continue;
    }

    if( strstartcmp(buffer,"struct") == 0) {
      runner=strtok(buffer,spacestr);
      runner=strtok(NULL,spacestr);
      if( runner == NULL) {
	log_full_error(WARNING,0,"Picked up a line called struct but no name - can't read struct!");
	ret = FALSE;
	continue;
      }
      sh = StructHolder_alloc_std();
      sh->name=stringalloc(runner);

      if(read_StructHolder_elements(sh,ifp) == FALSE) {
	log_full_error(WARNING,0,"Unable to read structure --- problemo!");
	sh = free_StructHolder(sh);
	ret = FALSE;

	continue;
      }
      add_Dynamite(out,wrap_StructHolder(sh));
    }
    else if ( strstartcmp(buffer,"matrix") == 0) {
      gm = read_GenericMatrix_line(buffer,ifp);
      
      if( gm == NULL ) {
	ret = FALSE;
	warn("Could not read GenericMatrix in line [%s]",buffer);
	continue;
      }
      
      add_Dynamite(out,wrap_GenericMatrix(gm));
    }
    else if ( strstartcmp(buffer,"display") == 0 ) {
      dis = read_Aln2Display_line(buffer,ifp);
      
      if( dis == NULL ) {
	ret = FALSE;
	warn("Could not read Display in line [%s]",buffer);
	continue;
      }
      
      add_Dynamite(out,wrap_Aln2Display(dis));
    }
    else if ( strstartcmp(buffer,"friend") == 0 ) {
      fr = read_Friend_line(buffer,ifp);
      if( fr == NULL ) {
	warn("Could not read friend line");
	continue;
      }
      add_Dynamite(out,wrap_Friend(fr));
    }
    else if ( strstartcmp(buffer,"api") == 0 ) {
      api = read_dynAPI_line(buffer,ifp);
      if( api == NULL  ){
	warn("Could not read API lines");
	continue;
      }
      add_Dynamite(out,wrap_dynAPI(api));
    }

    else if ( strstartcmp(buffer,"method") == 0 ) {
      add_me_MethodTypeSet(mts,read_Method_line(buffer,ifp));
    }
    else if ( strstartcmp(buffer,"type") == 0 ) {
      add_ty_MethodTypeSet(mts,read_Type_line(buffer,ifp));
    }
    else if( only_whitespace(buffer,spacestr) == TRUE ) {
      continue;
    }
  
    else {
      fatal("[Dynamite Level] Did not understand line [%s]. Probably a run-away parsing error, so failing now",buffer);
    }
  }

  if( feof(ifp) || ferror(ifp) ) {
    *end_of_file = TRUE;
  }
  else *end_of_file = FALSE;

  if( ret == FALSE )
    log_full_error(FATAL,0,"Could not read Dynamite instances");  

  return out;
  
}


DynamiteHolder * wrap_Aln2Display(Aln2Display * dis)
{
  DynamiteHolder * out;
  
  out = DynamiteHolder_alloc();
  
  out->data = dis;
  out->type = DYNAMITETYPEDISPLAY;

  return out;
}

DynamiteHolder * wrap_Friend(Friend * fr)
{
  DynamiteHolder * out;
  
  out = DynamiteHolder_alloc();
  
  out->data = fr;
  out->type = DYNAMITETYPEFRIEND;

  return out;
}

DynamiteHolder * wrap_dynAPI(dynAPI * api)
{
  DynamiteHolder * out;
  
  out = DynamiteHolder_alloc();
  
  out->data = api;
  out->type = DYNAMITETYPEAPI;

  return out;
}


DynamiteHolder * wrap_GenericMatrix(GenericMatrix * wm)
{
  DynamiteHolder * out;
  
  out = DynamiteHolder_alloc();
  
  out->data = wm;
  out->type = DYNAMITETYPEGENERICMATRIX;
  return out;
}

DynamiteHolder * wrap_StructHolder(StructHolder * wm)
{
  DynamiteHolder * out;
  
  out = DynamiteHolder_alloc();
  
  out->data = wm;
  out->type = DYNAMITETYPESTRUCT;

  return out;
}


/* ************************************************* */ 
/*  Memory allocation functions written by mite on   */ 
/*          Wed Nov 15 20:14:23 1995                 */ 
/*                                                   */ 
/* ************************************************* */ 


/* Simple alloc function for DynamiteHolder */ 
DynamiteHolder * DynamiteHolder_alloc(void) 
{
  DynamiteHolder * out; /* out is exported at end of function */ 
  
  
  /* call ckalloc and see if NULL */ 
  if((out=(DynamiteHolder *) ckalloc (sizeof(DynamiteHolder))) == NULL)   {
    log_full_error(WARNING,0,"DynamiteHolder_alloc failed ");  
    return NULL; /* calling function should respond! */ 
  }  
  out->data = NULL;  
  out->type = 0;  
  
  
  return out;  
}  


/* Free function for DynamiteHolder */ 
DynamiteHolder * free_DynamiteHolder(DynamiteHolder * obj) 
{


  if( obj == NULL)  
    return NULL;  
  
  
  ckfree(obj);  
  return NULL;  
}  


/* swap function for qsort function */ 
void swap_Dynamite(DynamiteHolder ** list,int i,int j)  
{
  register DynamiteHolder * temp;  
  temp=list[i];  
  list[i]=list[j];  
  list[j]=temp;  
}  


/* qsort - lifted from K&R - for Dynamite */ 
void qsort_Dynamite(DynamiteHolder ** list,int left,int right,int (*comp)(DynamiteHolder * ,DynamiteHolder * )) 
{
  int i,last;  
  if( left >= right )  
    return;  
  

  swap_Dynamite(list,left,(left+right)/2);  
  last = left;  
  for ( i=left+1; i <= right;i++)   {
    if( (*comp)(list[i],list[left]) < 0)  
      swap_Dynamite (list,++last,i);  
		}  
  swap_Dynamite (list,left,last);  
  qsort_Dynamite(list,left,last-1,comp);  
  qsort_Dynamite(list,last+1,right,comp);  
}  


/* sort function to be called  */ 
void sort_Dynamite(Dynamite * obj,int (*comp)(DynamiteHolder *, DynamiteHolder *)) 
{
  qsort_Dynamite(obj->list,0,obj->len-1,comp);  
  return;  
}  


/* will expand function if necessary */ 
boolean add_Dynamite(Dynamite * obj,DynamiteHolder * add) 
{
  if( obj->len >= obj->maxlen)   {
    if( expand_Dynamite(obj,obj->len + DynamiteLISTLENGTH) == FALSE)  
      return FALSE;  
  }  
  
  
  obj->list[obj->len++]=add;  
  return TRUE;  
}  


/* Expander function for Dynamite */ 
boolean expand_Dynamite(Dynamite * obj,int len) 
{
  
  
  if( obj->maxlen > len )   {
    log_full_error(PEDANTIC,0,"expand_Dynamite called with no need");  
    return TRUE;  
  }  
  
  
  if( (obj->list = (DynamiteHolder ** ) ckrealloc (obj->list,sizeof(DynamiteHolder *)*len)) == NULL)   {
    log_full_error(WARNING,0,"ckrealloc failed for expand_Dynamite, returning FALSE");  
    return FALSE;  
  }  
  obj->maxlen = len;  
  return TRUE;  
}  


Dynamite * Dynamite_alloc_std(void) 
{
  return Dynamite_alloc_len(DynamiteLISTLENGTH);  
}  


/* Alloc length function for Dynamite */ 
/* This function automatically allocates memory for list components */ 
Dynamite * Dynamite_alloc_len(int len) 
{
  Dynamite * out; /* out is exported at the end of function */ 


  /* Call alloc function: return NULL if NULL */ 
  /* Warning message alread in alloc function */ 
  if((out = Dynamite_alloc()) == NULL)  
    return NULL;  
  
  
  /* Calling ckcalloc for list elements */ 
  if((out->list = (DynamiteHolder ** ) ckcalloc (len,sizeof(DynamiteHolder *))) == NULL)   {
    log_full_error(WARNING,0,"Warning, ckcalloc failed in Dynamite_alloc_len");  
		return NULL;  
  }  
  out->len = 0;  
  out->maxlen = len;  
  
  
  return out;  
}  


/* Simple alloc function for Dynamite */ 
Dynamite * Dynamite_alloc(void) 
{
  Dynamite * out; /* out is exported at end of function */ 

  
  /* call ckalloc and see if NULL */ 
  if((out=(Dynamite *) ckalloc (sizeof(Dynamite))) == NULL)   {
    log_full_error(WARNING,0,"Dynamite_alloc failed ");  
    return NULL; /* calling function should respond! */ 
  }  
  out->list = NULL;  
  out->len = out->maxlen = 0;  
  out->name = NULL;  
  

  return out;  
}  


/* Free function for Dynamite */ 
Dynamite * free_Dynamite(Dynamite * obj) 
{
  register int i;  
  
  
  if( obj == NULL)  
    return NULL;  
  
  
  if( obj->list != NULL)   {
    for(i=0;i<obj->len;i++)   {
      if( obj->list[i] != NULL)  
	free_DynamiteHolder(obj->list[i]);  
    }  
  }  
  if( obj->name != NULL)  
    ckfree(obj->name);  
  
  
  ckfree(obj);  
  return NULL;  
}  

